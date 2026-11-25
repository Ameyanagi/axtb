# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run official xtb optimizations for the prepared optbench XYZ starts.

Example: python reshiki/optbench/xtb_opt.py --starts starts --potential gfn2 \
    --level normal --jobs 4 --out results/xtb-gfn2.jsonl

Each start runs once in its own temporary directory. Results, logs, and copied
geometries live beside --out. The parser targets official xtb 6.8.0 ANCopt.
The v6.8.0 upstream optimizer.f90, geoopt_driver.f90, scf_module.F90,
scc_core.f90, and type/timer.f90 define the corresponding output formats.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote

FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
CYCLE = re.compile(r"^\s*\.+\s+CYCLE\s+(\d+)\s+\.+\s*$", re.MULTILINE | re.IGNORECASE)
CONVERGED = re.compile(r"GEOMETRY OPTIMIZATION CONVERGED AFTER\s+(\d+)\s+ITERATIONS", re.IGNORECASE)
MAX_STEPS = re.compile(r"FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN\s+(\d+)\s+ITERATIONS", re.IGNORECASE)
CYCLE_END = re.compile(
    r"GEOMETRY OPTIMIZATION CONVERGED AFTER|FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN|"
    r"\bANCopt \(total\)|\bFinal Singlepoint\b|\b(?:abnormal|normal) termination\b|\[ERROR\]",
    re.IGNORECASE,
)
CYCLE_ENERGY = rf"^\s*\*\s*total energy\s*:\s*({FLOAT})\s+Eh\b"
CYCLE_GRADIENT = rf"^\s*gradient norm\s*:\s*({FLOAT})\s+Eh/"
SCF_HEADER = re.compile(r"^\s*iter\s+E\s+dE\s+RMSdq\s+gap\s+omega\s+full diag\s*$", re.IGNORECASE)
SCF_ROW = re.compile(rf"^\s*(\d+)\s+({FLOAT})\s*({FLOAT})\s*({FLOAT})\s*({FLOAT})\s*({FLOAT})\s+[TF]\s*$")
KINDS = ("etkdg", "mmff", "given", "rough")


@dataclass(frozen=True)
class Start:
    path: Path
    set: str
    case: str
    kind: str


def _float(value: str) -> float | None:
    number = float(value.replace("D", "E").replace("d", "e"))
    return number if math.isfinite(number) else None


def _last_float(pattern: str, text: str) -> float | None:
    matches = re.findall(pattern, text, re.MULTILINE | re.IGNORECASE)
    return _float(matches[-1]) if matches else None


def _cycle_trace(text: str, cycles: Sequence[re.Match[str]]) -> list[list]:
    """Read only complete energy/gradient pairs inside their CYCLE blocks."""
    trace = []
    for index, cycle in enumerate(cycles):
        end = cycles[index + 1].start() if index + 1 < len(cycles) else len(text)
        block = CYCLE_END.split(text[cycle.end() : end], maxsplit=1)[0]
        energy = _last_float(CYCLE_ENERGY, block)
        gradient_norm = _last_float(CYCLE_GRADIENT, block)
        if energy is not None and gradient_norm is not None:
            # xtb prints neither per-cycle wall times nor a maximum gradient.
            trace.append([int(cycle[1]), 0, energy, None, "normal", gradient_norm])
    return trace


def _scf_tables(text: str) -> list[int | None]:
    """Count each SCC table once, including the quieter tables inside ANCopt."""
    tables: list[int | None] = []
    in_table = False
    for line in text.splitlines():
        if SCF_HEADER.match(line):
            tables.append(None)
            in_table = True
        elif in_table:
            row = SCF_ROW.match(line)
            if row:
                tables[-1] = int(row[1])
            elif line.strip() and tables[-1] is not None:
                in_table = False
    return tables


def _timer(text: str, label: str, *, total: bool = False) -> float | None:
    if total:
        match = re.search(
            rf"^\s*{re.escape(label)}\s+(\d+)\s+d,\s*(\d+)\s+h,\s*(\d+)\s+min,\s*({FLOAT})\s+sec\s*$",
            text,
            re.MULTILINE | re.IGNORECASE,
        )
        if match:
            seconds = _float(match[4])
            if seconds is not None:
                return int(match[1]) * 86400 + int(match[2]) * 3600 + int(match[3]) * 60 + seconds
        return None
    match = re.search(
        rf"^\s*{re.escape(label)}\s+\.\.\.\s*(\d+)\s+min,\s*({FLOAT})\s+sec\b",
        text,
        re.MULTILINE | re.IGNORECASE,
    )
    if match:
        seconds = _float(match[2])
        if seconds is not None:
            return int(match[1]) * 60 + seconds
    return None


def parse_output(text: str, *, potential: str, returncode: int = 0) -> dict:
    """Parse printed measurements without treating gradient norms as maxima.

    ``steps`` uses the explicit convergence or failure report when available.
    ``evaluations`` and ``cycle_blocks`` count printed CYCLE blocks, including
    incomplete blocks. Initial and final singlepoints remain in the independent
    GFN2 SCC iteration total. Trace rows use the cycle energy and gradient norm;
    unavailable timestamps are zero and unavailable maximum gradients are null.
    """
    cycles = list(CYCLE.finditer(text))
    converged = CONVERGED.findall(text)
    exhausted = MAX_STEPS.findall(text)
    reported_iterations = int((converged or exhausted)[-1]) if converged or exhausted else None
    scf_tables = _scf_tables(text) if potential == "gfn2" else []
    normal = bool(re.search(r"\bnormal termination\b", text, re.IGNORECASE))
    scf_failed = bool(
        re.search(
            r"SCF not converged|SCC.{0,30}not converged|convergence criteria cannot be satisfied|"
            r"Self consistent charge iterator terminated",
            text,
            re.IGNORECASE,
        )
    )
    if scf_failed:
        status = "scf_failed"
    elif exhausted:
        status = "max_steps"
    elif returncode or re.search(r"abnormal termination|\[ERROR\]", text, re.IGNORECASE):
        status = "error"
    elif converged and normal:
        status = "converged"
    else:
        status = "error"

    complete_scf = normal and len(scf_tables) == len(cycles) + 2 and all(n is not None for n in scf_tables)
    ancopt_seconds = _timer(text, "ANCopt (total)", total=True)
    potential_seconds = _timer(text, "single point calculation")
    model_hessian_seconds = _timer(text, "model hessian")
    optimizer_seconds = None
    if ancopt_seconds is not None and potential_seconds is not None:
        # Independent millisecond rounding of printed timers can give -0.001.
        overhead = ancopt_seconds - potential_seconds
        optimizer_seconds = max(0.0, overhead) if overhead >= -0.002 else None
    energy = _last_float(rf"^\s*\|\s*TOTAL ENERGY\s+({FLOAT})\s+Eh\b", text)
    if energy is None:
        energy = _last_float(CYCLE_ENERGY, text)
    gradient_norm = _last_float(rf"^\s*\|\s*GRADIENT NORM\s+({FLOAT})\s+Eh/", text)
    if gradient_norm is None:
        gradient_norm = _last_float(CYCLE_GRADIENT, text)
    return {
        "status": status,
        "steps": reported_iterations if reported_iterations is not None else len(cycles) or None,
        "evaluations": len(cycles),
        "hessian_evaluations": None,
        "potential_seconds": potential_seconds,
        "hessian_seconds": None,
        "optimizer_seconds": optimizer_seconds,
        "scf_iterations_total": sum(scf_tables) if complete_scf else None,
        "energy": energy,
        "max_gradient": None,
        "trace": _cycle_trace(text, cycles),
        "cycle_blocks": len(cycles),
        "reported_iterations": reported_iterations,
        "gradient_norm": gradient_norm,
        "ancopt_seconds": ancopt_seconds,
        "model_hessian_seconds": model_hessian_seconds,
    }


def read_xyz(path: Path) -> tuple[list[str], dict[str, str]]:
    """Validate a single XYZ frame and return its ordered element symbols."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: incomplete XYZ header")
    try:
        count = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"{path}: invalid XYZ atom count") from exc
    if count < 1 or len(lines) < count + 2 or any(line.strip() for line in lines[count + 2 :]):
        raise ValueError(f"{path}: expected one XYZ frame of {count} atoms")
    symbols = []
    for line in lines[2 : count + 2]:
        fields = line.split()
        if len(fields) != 4 or not re.fullmatch(r"[A-Za-z]{1,2}", fields[0]):
            raise ValueError(f"{path}: invalid XYZ atom row {line!r}")
        try:
            valid = all(_float(value) is not None for value in fields[1:])
        except ValueError:
            valid = False
        if not valid:
            raise ValueError(f"{path}: invalid or nonfinite XYZ coordinates")
        symbols.append(fields[0].capitalize())
    metadata = dict(token.split("=", 1) for token in lines[1].split() if "=" in token)
    return symbols, metadata


def discover_starts(
    directory: Path,
    *,
    sets: Sequence[str] | None = None,
    cases: Sequence[str] | None = None,
    kinds: Sequence[str] | None = None,
) -> list[Start]:
    starts = []
    for path in sorted(directory.glob("*/*.xyz")):
        name, separator, kind = path.stem.rpartition(".")
        if not separator or not name or kind not in KINDS:
            continue
        if sets and path.parent.name not in sets:
            continue
        if cases and name not in cases and f"{path.parent.name}/{name}" not in cases:
            continue
        if kinds and kind not in kinds:
            continue
        starts.append(Start(path.resolve(), path.parent.name, name, kind))
    return starts


def _relative_artifact(start: Start, potential: str, config: str, category: str, suffix: str) -> Path:
    return Path(
        category, potential, config, quote(start.set, safe="-_"), f"{quote(start.case, safe='-_')}.{start.kind}{suffix}"
    )


def run_start(
    start: Start, *, xtb: str, potential: str, level: str, output_parent: Path, timeout: float | None
) -> dict:
    config = f"xtb-opt-{level}"
    result = {
        "case": start.case,
        "set": start.set,
        "start": start.kind,
        "potential": potential,
        "config": config,
        "atoms": None,
        **parse_output("", potential=potential),
        "wall_seconds": None,
        "final_xyz": None,
        "returncode": None,
    }
    try:
        symbols, metadata = read_xyz(start.path)
        result["atoms"] = len(symbols)
        charge, uhf = int(metadata["charge"]), int(metadata["uhf"])
        if uhf < 0:
            raise ValueError("uhf must be nonnegative")
    except (OSError, ValueError, KeyError) as exc:
        result.update(status="input_error", error=str(exc))
        return result

    command = [xtb, "input.xyz", *(["--gfn", "2"] if potential == "gfn2" else ["--gfnff"])]
    command += ["--opt", level, "--chrg", str(charge), "--uhf", str(uhf)]
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OMP_THREAD_LIMIT="1",
        OMP_STACKSIZE="1G",
        MKL_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )
    log_relative = _relative_artifact(start, potential, config, "logs", ".log")
    log_path = output_parent / log_relative
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="optbench-xtb-") as scratch:
            work = Path(scratch)
            # Match the reference runner: avoid an inherited user .xtbrc or
            # parameter search path changing the optimizer or potential.
            env.pop("XTBPATH", None)
            env["XTBHOME"] = str(work.resolve())
            shutil.copyfile(start.path, work / "input.xyz")
            started = time.perf_counter()
            try:
                completed = subprocess.run(
                    command,
                    cwd=work,
                    env=env,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=timeout,
                    check=False,
                )
                result["wall_seconds"] = time.perf_counter() - started
                text = completed.stdout
                result.update(parse_output(text, potential=potential, returncode=completed.returncode))
                result["returncode"] = completed.returncode
                if result["status"] in ("error", "scf_failed"):
                    result["error"] = f"xtb exited {completed.returncode}. See {log_relative.as_posix()}"
            except subprocess.TimeoutExpired as exc:
                result["wall_seconds"] = time.perf_counter() - started
                text = exc.stdout or ""
                if isinstance(text, bytes):
                    text = text.decode("utf-8", errors="replace")
                result.update(
                    parse_output(text, potential=potential), status="timeout", error=f"exceeded {timeout:g} seconds"
                )
            log_path.write_text(text, encoding="utf-8")
            result["log"] = log_relative.as_posix()
            final = work / "xtbopt.xyz"
            if not final.is_file() and result["status"] in ("max_steps", "scf_failed"):
                final = work / "xtblast.xyz"
            if final.is_file():
                final_symbols, _ = read_xyz(final)
                if final_symbols != symbols:
                    raise ValueError(f"{final.name} atom identities/order differ from the input")
                relative = _relative_artifact(start, potential, config, "geoms", ".xyz")
                destination = output_parent / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(final, destination)
                result["final_xyz"] = relative.as_posix()
            elif result["status"] == "converged":
                raise ValueError("convergence reported without xtbopt.xyz")
            if (work / "NOT_CONVERGED").exists() and result["status"] == "converged":
                result.update(status="error", error="xtb wrote NOT_CONVERGED despite its convergence report")
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        result.update(status="error", error=str(exc))
    return result


def _comma_list(value: str) -> list[str]:
    items = [item.strip() for item in value.split(",")]
    if not items or any(not item for item in items):
        raise argparse.ArgumentTypeError("provide a nonempty comma-separated list")
    return list(dict.fromkeys(items))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--starts", type=Path, required=True, help="prepared <set>/<case>.<kind>.xyz directory")
    parser.add_argument("--potential", choices=("gfn2", "gfnff"), required=True)
    parser.add_argument("--xtb", default=os.environ.get("AXTB_XTB680"), help="xtb executable (default: $AXTB_XTB680)")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True, help="JSONL file to replace on each invocation")
    parser.add_argument("--level", choices=("normal", "tight"), default="normal")
    parser.add_argument("--sets", type=_comma_list, help="exact comma-separated set names")
    parser.add_argument("--cases", type=_comma_list, help="comma-separated case names or set/name pairs")
    parser.add_argument("--starts-kind", type=_comma_list, help=f"comma-separated start kinds: {', '.join(KINDS)}")
    parser.add_argument("--timeout", type=float, help="optional wall seconds allowed per job")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.timeout is not None and (not math.isfinite(args.timeout) or args.timeout <= 0):
        parser.error("--timeout must be finite and positive")
    if not args.starts.is_dir():
        parser.error(f"--starts is not a directory: {args.starts}")
    if args.starts_kind and set(args.starts_kind) - set(KINDS):
        parser.error(f"--starts-kind must use {', '.join(KINDS)}")
    if not args.xtb:
        parser.error("provide --xtb or set AXTB_XTB680")
    xtb = shutil.which(os.path.expanduser(args.xtb))
    if xtb is None:
        parser.error(f"xtb is not an executable file or command: {args.xtb}")
    xtb = str(Path(xtb).absolute())
    starts = discover_starts(args.starts, sets=args.sets, cases=args.cases, kinds=args.starts_kind)
    if not starts:
        parser.error("no XYZ starts matched the requested filters")
    output = args.out.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with (
        output.open("w", encoding="utf-8") as stream,
        concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool,
    ):
        # map preserves the deterministic, sorted start order in the JSONL.
        results = pool.map(
            lambda start: run_start(
                start,
                xtb=xtb,
                potential=args.potential,
                level=args.level,
                output_parent=output.parent,
                timeout=args.timeout,
            ),
            starts,
        )
        for result in results:
            stream.write(json.dumps(result, allow_nan=False, separators=(",", ":")) + "\n")
            stream.flush()
            failures += result["status"] != "converged"
            print(f"{result['set']}/{result['case']}.{result['start']}: {result['status']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
