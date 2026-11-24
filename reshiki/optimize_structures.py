# SPDX-License-Identifier: LGPL-3.0-or-later
"""Prepare method-specific reference geometries for Reshiki Hessian checks.

Run ``python reshiki/optimize_structures.py`` after configuring ``AXTB_XTB680``
and ``AXTB_GXTB2``. Only the reference binaries needed by the selected methods
are required; tblite and axtb are not used to optimize these structures.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

from cases import METHODS, ROOT, load_cases, optimized_path, read_structure
from references import ReferenceError, optimize_structure, reference_programs


def _comma_list(value: str) -> list[str]:
    items = [item.strip() for item in value.split(",")]
    if not items or any(not item for item in items):
        raise argparse.ArgumentTypeError("provide a nonempty comma-separated list")
    return list(dict.fromkeys(items))


def _source_name(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.name


def _write_geometry(path: Path, result: dict, case: dict, method: str, reference: dict) -> None:
    program = "g-xTB 2.0.1, xtb" if method == "gxtb" else "xtb"
    comment = f"charge={case['charge']} uhf={case['uhf']} optimized with {program} {reference['version']} --opt tight"
    lines = [str(len(result["symbols"])), comment]
    lines.extend(
        f"{symbol:2s} {x:.17g} {y:.17g} {z:.17g}" for symbol, (x, y, z) in zip(result["symbols"], result["positions"])
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix="opt-") as directory:
        temporary = Path(directory) / path.name
        temporary.write_text("\n".join([*lines, ""]), encoding="utf-8")
        temporary.replace(path)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", help="case ID substring, matched without case sensitivity")
    parser.add_argument("--methods", type=_comma_list, default=list(METHODS), help="comma-separated methods")
    parser.add_argument("--overwrite", action="store_true", help="replace existing optimized geometries")
    parser.add_argument(
        "--timeout", type=float, default=1800.0, help="seconds allowed per optimization (default: 1800)"
    )
    args = parser.parse_args(argv)
    if set(args.methods) - set(METHODS):
        parser.error(f"--methods must use {', '.join(METHODS)}")
    if not 0 < args.timeout < float("inf"):
        parser.error("--timeout must be finite and positive")
    try:
        cases = load_cases(args.only)
        jobs = [
            (case, method, optimized_path(case, method))
            for case in cases
            if case["hessian"]
            for method in args.methods
            if case["methods"].get(method) == "match"
        ]
        if not jobs:
            print("No Hessian geometries are selected.")
            return 0
        required = tuple(dict.fromkeys("gxtb" if method == "gxtb" else "xtb" for _, method, _ in jobs))
        programs, metadata = reference_programs(required=required)
    except (ReferenceError, ValueError, OSError) as exc:
        print(f"Geometry preparation failed: {exc}", file=sys.stderr)
        return 1
    failures = 0
    for case, method, destination in jobs:
        if destination.exists() and not args.overwrite:
            print(f"{case['id']} {method}: exists ({_source_name(destination)}); use --overwrite to replace")
            continue
        name = "gxtb" if method == "gxtb" else "xtb"
        try:
            structure = read_structure(case["path"])
            result = optimize_structure(
                programs[name], method, structure, case["charge"], case["uhf"], timeout=args.timeout
            )
            _write_geometry(destination, result, case, method, metadata[name])
        except (ReferenceError, ValueError, OSError) as exc:
            failures += 1
            print(f"{case['id']} {method}: FAILED: {exc}", file=sys.stderr)
            continue
        print(f"{case['id']} {method}: wrote {_source_name(destination)} ({result['energy']:.12f} Eh)")
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())
