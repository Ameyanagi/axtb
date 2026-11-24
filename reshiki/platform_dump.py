# SPDX-License-Identifier: LGPL-3.0-or-later
"""Collect a portable, single-threaded axtb CLI dump for cross-platform comparison.

Run ``python reshiki/platform_dump.py --out linux-x64.json`` on the baseline,
then use ``platform_compare.py`` to compare another host's dump. ``--only water``
selects matching case IDs for a smoke run. No reference program or Python binding
is used. Expected unsupported inputs remain in the dump as structured errors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from axtb_cli import AxtbError, axtb_binary, axtb_metadata, calculation_payload, host_metadata, request
from cases import METHODS, ROOT, hessian_path, load_cases, read_structure

SCHEMA_VERSION = 1
RESULT_FIELDS = ("energy", "gradient", "charges", "iterations")
IDENTITY_FIELDS = ("method", "conventions", "numerics", "eigensolver")
VIBRATIONAL_FIELDS = ("frequencies", "linear", "n_imaginary", "vibrational_conventions", "hessian_method")
OPERATIONAL_ERRORS = {"operational_error", "process", "parse", "timeout", "eigensolver"}
CONVENTIONS = "axtb"


def rustc_version() -> str | None:
    """Return the available compiler version, without making Rust a prerequisite."""
    try:
        proc = subprocess.run(["rustc", "--version"], capture_output=True, text=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def display_path(path: Path) -> str:
    """Use checkout-relative paths so geometry identities travel between hosts."""
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def input_fingerprint(payload: dict[str, Any], run: str) -> str:
    """Hash the effective request identity and the exact coordinates sent to the CLI."""
    identity = {name: payload[name] for name in ("method", "conventions", "charge", "uhf", "numbers", "positions")}
    identity["run"] = run
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def collect_record(binary: str, case: dict[str, Any], method: str, hessian: bool = False) -> dict[str, Any]:
    """Run one CLI request, retaining diagnostics and the quantities being compared."""
    path = hessian_path(case, method) if hessian else Path(case["path"])
    record: dict[str, Any] = {
        "case": case["id"],
        "method": method,
        "run": "hessian-analytic" if hessian else "gradient",
        "geometry": display_path(path),
        "expected": case["methods"][method],
        "input_fingerprint": None,
    }
    try:
        structure = read_structure(path)
        settings: dict[str, Any] = {"gradient": True, "conventions": CONVENTIONS}
        if hessian:
            settings.update(hessian=True, hessian_method="analytic", thermo=False)
        payload = calculation_payload(structure, case["charge"], case["uhf"], method, **settings)
        record["input_fingerprint"] = input_fingerprint(payload, record["run"])
        response = request(binary, payload)
    except AxtbError as exc:
        record.update(error_kind=exc.kind, error=exc.message)
        if exc.response is not None:
            record["error_response"] = exc.response
    except (OSError, ValueError, KeyError, RuntimeError, subprocess.SubprocessError) as exc:
        record.update(error_kind="operational_error", error=str(exc), error_type=type(exc).__name__)
    else:
        result = {name: response.get(name) for name in RESULT_FIELDS}
        result.update({name: response[name] for name in IDENTITY_FIELDS if name in response})
        if hessian:
            # Default axtb analysis diagonalizes the rigid-body complement, so every
            # reported frequency is vibrational (including negative/imaginary modes).
            # Do not drop near-zero values or assume that only positive modes count.
            result.update({name: response[name] for name in VIBRATIONAL_FIELDS if name in response})
        record["result"] = result
        issues = []
        if "eigensolver" in response and response["eigensolver"] != "faer":
            issues.append(f"expected faer eigensolver, received {response['eigensolver']!r}")
        for name in ("energy", "gradient", "charges"):
            if result[name] is None:
                issues.append(f"missing {name} in CLI response")
        if hessian and "frequencies" not in result:
            issues.append("missing frequencies in analytic Hessian response")
        if issues:
            record["issues"] = issues
    return record


def make_dump(binary: str, cases: list[dict[str, Any]]) -> dict[str, Any]:
    """Collect all method/case pairs, including expected failures, in stable order."""
    metadata = dict(host_metadata())
    metadata.update(
        axtb=axtb_metadata(binary),
        rustc_version=rustc_version(),
        generated_at=datetime.now(timezone.utc).isoformat(),
    )
    records = []
    for case in cases:
        for method in METHODS:
            for hessian in (False, True) if case["hessian"] else (False,):
                record = collect_record(binary, case, method, hessian)
                records.append(record)
                status = record.get("error_kind", "ok")
                print(f"{case['id']} {method} {record['run']}: {status}", file=sys.stderr, flush=True)
    return {
        "schema_version": SCHEMA_VERSION,
        "metadata": metadata,
        "selections": [
            {
                "case": case["id"],
                "methods": list(METHODS),
                "runs": ["gradient", "hessian-analytic"] if case["hessian"] else ["gradient"],
            }
            for case in cases
        ],
        "settings": {
            "methods": list(METHODS),
            "threads": 1,
            "conventions": CONVENTIONS,
            "input_fingerprint": "sha256 of canonical JSON request identity, including exact bohr coordinates",
            "AXTB_LAPACK": None,
            "AXTB_LAPACK_BLAS": None,
            "frequency_modes": "vibrational only; all default-convention frequencies retained",
            "units": {"energy": "Eh", "gradient": "Eh/bohr", "frequencies": "cm^-1", "charges": "e"},
        },
        "records": records,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="destination JSON file")
    parser.add_argument("--only", metavar="SUBSTR", help="select case IDs containing this substring")
    args = parser.parse_args(argv)
    try:
        cases = load_cases(only=args.only)
        if not cases:
            parser.error("no cases matched --only")
        dump = make_dump(axtb_binary(), cases)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(dump, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"platform dump failed: {exc}", file=sys.stderr)
        return 1
    records = dump["records"]
    errors = sum("error_kind" in record for record in records)
    invalid = sum(bool(record.get("issues")) or record.get("error_kind") in OPERATIONAL_ERRORS for record in records)
    print(f"Wrote {len(records)} runs ({errors} errors, {invalid} invalid/operational runs) to {args.out}")
    return int(invalid > 0)


if __name__ == "__main__":
    raise SystemExit(main())
