# SPDX-License-Identifier: LGPL-3.0-or-later
"""Run ReShiki's documented reference gates on Linux x64."""

import argparse
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from axtb_cli import AxtbError, axtb_binary, axtb_metadata, host_metadata, request, run_axtb
from cases import (
    ATOMIC_NUMBERS,
    BOHR_2018,
    CUTOFF_BOHR,
    EXPECTED_BOTH_FAIL,
    METHODS,
    REQUIRED_QUANTITIES,
    ROOT,
    hessian_path,
    load_cases,
    read_structure,
    tolerances,
    vibrational_mode_count,
)
from references import ReferenceRefusal, parse_molden, reference_programs, run_tblite, run_xtb

STATUSES = ("pass", "fail", "expected_error", "expected_difference", "both_failed", "error")
HESSIAN_CHECKS = ("hessian_analytic", "hessian_same", "hessian_numerical")
CONVENTIONS = {"gfn2": "xtb", "gfnff": "xtb", "gxtb": "gxtb"}
ANALYSIS = {**CONVENTIONS, "gxtb": "xtb-6.7.1"}
VALUE_SETTINGS = {"gfn2": {"accuracy": 1e-6}, "gxtb": {"accuracy": 1e-4, "max_iter": 250}}
MOLDEN_ACCURACY = {"gfn2": 1e-4, "gxtb": 1.0}
HESSIAN_ACCURACY = {"gfn2": 1e-6}


def attempt(function, *args, **kwargs):
    try:
        return function(*args, **kwargs), None
    except (OSError, RuntimeError, ValueError, KeyError, IndexError, TypeError, subprocess.SubprocessError) as exc:
        return None, exc


def error_details(exc):
    details = {"kind": getattr(exc, "kind", type(exc).__name__), "message": str(exc)}
    if isinstance(exc, AxtbError) and exc.response is not None:
        details["response"] = exc.response
    return details


def max_difference(left, right, shape=None):
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != right.shape or (shape is not None and left.shape != shape):
        raise ValueError(f"incompatible shapes {left.shape} and {right.shape} (expected {shape})")
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError("non-finite values in response")
    return float(np.max(np.abs(left - right), initial=0))


def orbital_difference(left, right, method, uhf):
    left, right = dict(left), dict(right)
    if uhf == 0 or method == "gfn2":
        # Either program can write a restricted set twice.
        for one, other in ((left, right), (right, left)):
            if "beta" in other and one.keys() == {"alpha"}:
                one["beta"] = one["alpha"]
    elif method == "gxtb" and ("beta" not in left or "beta" not in right):
        raise ValueError("open-shell g-xTB requires beta orbitals in both Molden files")
    if left.keys() != right.keys() or "alpha" not in left:
        raise ValueError(f"different Molden spin channels: {list(left)} vs {list(right)}")
    worst = 0.0
    for spin in left:
        a, b = np.sort(left[spin]), np.sort(right[spin])
        if len(a) < len(b):
            raise ValueError(f"{spin}: axtb has {len(a)} orbitals, reference has {len(b)}")
        # The reference Molden writer may omit orbitals above 2 Eh.
        extra = a[len(b) :]
        if not (np.isfinite(extra) & (extra > 2.0)).all():
            raise ValueError(f"{spin}: the reference may omit only orbitals above 2 Eh")
        worst = max(worst, max_difference(a[: len(b)], b))
    return worst


def frequency_difference(left, right, positions, method, check):
    a, b = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if any(freq.ndim != 1 or not np.isfinite(freq).all() for freq in (a, b)):
        raise ValueError("expected finite frequency vectors")
    expected = vibrational_mode_count(positions)
    if len(a) != expected:
        raise ValueError(f"axtb has {len(a)} vibrational modes, expected {expected}")
    if method != "gfnff" and len(b) < expected:
        # xtb's g98.out omits modes at or below 0.1 cm⁻¹ in magnitude (src/freq/io.f90).
        keep = np.abs(a) > 0.1
        if np.count_nonzero(keep) == len(b):
            a = a[keep]
    if len(a) != len(b):
        raise ValueError(f"reference has {len(b)} vibrational modes, incompatible with {expected} axtb modes")
    a, b = np.sort(a), np.sort(b)
    dropped = expected - len(a)
    if (method, check) == ("gfn2", "hessian_analytic"):
        keep = np.abs(b) > 50  # Match the reference matrix's check_analytic_xtb selection.
        a, b = a[keep], b[keep]
    return max_difference(a, b), dropped


def measure(left, right, row, nat, uhf):
    shapes = {"energy": (), "gradient": (nat, 3), "charges": (nat,), "dipole": (3,)}
    measured = row["deviations"]
    for quantity in row["tolerances"]:
        a, b = left.get(quantity), right.get(quantity)
        if a is None or b is None:
            raise ValueError(f"missing {quantity} from {'axtb' if a is None else row['reference']}")
        if quantity == "orbital_energies":
            measured[quantity] = orbital_difference(a, b, row["method"], uhf)
        elif quantity == "frequencies":
            measured[quantity], row["axtb_modes_dropped"] = frequency_difference(
                a, b, right["positions"], row["method"], row["check"]
            )
        else:
            measured[quantity] = max_difference(a, b, shapes[quantity])


def compare(case, method, check, reference, path, nat, left, right=None, *, long=False):
    lhs, left_error = left
    rhs, right_error = right if right is not None else (None, None)
    row = {
        "case": case["id"],
        "method": method,
        "check": check,
        "geometry": str(path.relative_to(ROOT)),
        "expected": case["methods"][method],
        "reference": reference,
        "status": "error",
        "deviations": {},
        "tolerances": tolerances(case, method, check, long),
        "iterations": {
            "axtb": (lhs or {}).get("iterations"),
            "reference": (rhs or {}).get("iterations"),
            "status": "not_gated",
        },
        "eigensolver": {"axtb": (lhs or {}).get("eigensolver"), "reference": (rhs or {}).get("eigensolver")},
    }
    errors = {side: error_details(exc) for side, exc in (("axtb", left_error), ("reference", right_error)) if exc}
    if errors:
        row["errors"] = errors
    expected_both_fail = EXPECTED_BOTH_FAIL.get((case["id"], method, check))

    if row["expected"] == "both_error":
        if (
            isinstance(left_error, AxtbError)
            and left_error.kind == "unsupported_element"
            and (left_error.response or {}).get("z") == 92
            and isinstance(right_error, ReferenceRefusal)
            and 92 in right_error.unsupported_elements
        ):
            row.update(status="expected_error", note="Both programs refused uranium (Z=92)")
        else:
            row["note"] = "Expected both programs to refuse uranium (Z=92)"
            succeeded = [side for side in ("axtb", "reference") if side not in errors]
            if succeeded:
                row["note"] += ". Calculation succeeded for " + ", ".join(succeeded)
    elif (
        expected_both_fail
        and isinstance(left_error, AxtbError)
        and left_error.kind == "scf_not_converged"
        and isinstance(right_error, ReferenceRefusal)
        and re.search(expected_both_fail["diagnostic"], str(right_error), re.MULTILINE)
    ):
        row.update(status="both_failed", note=f"Both SCFs failed to converge ({expected_both_fail['source']})")
    elif method == "gfnff" and check == "molden" and right_error is None:
        if isinstance(left_error, AxtbError) and left_error.kind == "not_applicable":
            row.update(status="pass", note="axtb refused GFN-FF Molden with not_applicable")
        elif left_error is None:
            row.update(status="fail", note="axtb accepted GFN-FF Molden")
    elif not errors:
        try:
            measure(lhs, rhs, row, nat, case["uhf"])
        except (TypeError, ValueError, KeyError) as exc:
            row["errors"] = {"comparison": error_details(exc)}
        else:
            measured = row["deviations"]
            # tests/gate.py: the documented limits have two significant digits, and a value passes
            # when it is at most its limit at two digits
            rounded = {q: float(f"{val:.1e}") for q, val in measured.items()}
            row["deviations_at_source_precision"] = rounded
            failed = [
                q for q, spec in row["tolerances"].items() if spec["limit"] is not None and rounded[q] > spec["limit"]
            ]
            row["status"] = "fail" if failed else "expected_difference" if long else "pass"
            if failed:
                row["note"] = "Tolerance exceeded for " + ", ".join(failed)
            elif long:
                row["note"] = "Documented real-space cutoff difference at the source's two-significant-figure precision"
    return row


def axtb_molden(binary, structure, case, method):
    settings = {"accuracy": MOLDEN_ACCURACY[method]} if method in MOLDEN_ACCURACY else {}
    result = run_axtb(
        binary, structure, case["charge"], case["uhf"], method, conventions=CONVENTIONS[method], molden=True, **settings
    )
    if method != "gfnff":
        result["orbital_energies"] = parse_molden(result["molden"])
    return result


def reference_molden(binary, structure, case, method):
    result = run_xtb(binary, method, structure, case["charge"], case["uhf"], "molden")
    result["orbital_energies"] = result["orbitals"]
    return result


def analyze_reference(binary, structure, reference, method):
    symbols, _ = structure
    return request(
        binary,
        {
            "numbers": [ATOMIC_NUMBERS[s] for s in symbols],
            "positions": np.asarray(reference["positions"]).tolist(),
            "analyze_hessian": np.asarray(reference["hessian"]).tolist(),
            "analysis_conventions": ANALYSIS[method],
            "thermo": False,
        },
    )


def check_case(job):
    case, selected_methods, binaries = job
    binary = binaries["axtb"]
    path = case["path"]
    structure = read_structure(path)
    nat = len(structure[0])
    charge, uhf = case["charge"], case["uhf"]
    rows = []
    for method in selected_methods:
        reference = "gxtb" if method == "gxtb" else "xtb"
        ref_binary = binaries[reference]
        values_settings = VALUE_SETTINGS.get(method, {})
        ax_values = attempt(
            run_axtb,
            binary,
            structure,
            charge,
            uhf,
            method,
            gradient=True,
            conventions=CONVENTIONS[method],
            **values_settings,
        )
        ref_values = attempt(run_xtb, ref_binary, method, structure, charge, uhf, "values")
        rows.append(compare(case, method, f"{reference}_values", reference, path, nat, ax_values, ref_values))
        if case["methods"][method] == "both_error":
            continue

        if method == "gfn2":
            ax_tblite = attempt(
                run_axtb, binary, structure, charge, uhf, method, gradient=True, conventions="tblite", **values_settings
            )
            ref_tblite = attempt(run_tblite, binaries["tblite"], structure, charge, uhf)
            rows.append(compare(case, method, "tblite_values", "tblite", path, nat, ax_tblite, ref_tblite))
            # The upstream matrix passes the same bohr coordinates to both axtb modes.
            bohr = structure[1] / BOHR_2018
            ax_default = attempt(
                run_axtb, binary, structure, charge, uhf, method, positions_bohr=bohr, gradient=True, **values_settings
            )
            extent = np.linalg.norm(bohr[:, None, :] - bohr[None, :, :], axis=2).max()
            rows.append(
                compare(
                    case,
                    method,
                    "default_vs_xtb",
                    "axtb[xtb]",
                    path,
                    nat,
                    ax_default,
                    ax_values,
                    long=bool(extent > CUTOFF_BOHR),
                )
            )

        ax_molden = attempt(axtb_molden, binary, structure, case, method)
        ref_molden = attempt(reference_molden, ref_binary, structure, case, method) if method != "gfnff" else None
        rows.append(
            compare(case, method, "molden", reference if ref_molden else None, path, nat, ax_molden, ref_molden)
        )

        if case["hessian"]:
            optimized = hessian_path(case, method)
            geometry, geometry_error = attempt(read_structure, optimized)
            if geometry_error:
                rows.extend(
                    compare(case, method, check, reference, optimized, nat, (None, geometry_error))
                    for check in HESSIAN_CHECKS
                )
                continue
            ref_hessian = attempt(run_xtb, ref_binary, method, geometry, charge, uhf, "hessian")
            hess_settings = {"hessian_accuracy": HESSIAN_ACCURACY[method]} if method in HESSIAN_ACCURACY else {}
            ax_analytic = attempt(
                run_axtb,
                binary,
                geometry,
                charge,
                uhf,
                method,
                hessian=True,
                hessian_method="analytic",
                thermo=False,
                conventions=CONVENTIONS[method],
                analysis_conventions=ANALYSIS[method],
                **hess_settings,
            )
            ax_same = (
                attempt(analyze_reference, binary, geometry, ref_hessian[0], method) if not ref_hessian[1] else None
            )
            rows.append(compare(case, method, "hessian_analytic", reference, optimized, nat, ax_analytic, ref_hessian))
            rows.append(
                compare(case, method, "hessian_same", reference, optimized, nat, ax_same or (None, None), ref_hessian)
            )
            if method == "gfnff":
                # xtb's own procedure (central differences, step 0.005 bohr) in axtb: separates
                # xtb's finite-difference error from the analytic Hessian's
                ax_numerical = attempt(
                    run_axtb,
                    binary,
                    geometry,
                    charge,
                    uhf,
                    method,
                    hessian=True,
                    hessian_method="numerical",
                    thermo=False,
                    conventions=CONVENTIONS[method],
                    analysis_conventions=ANALYSIS[method],
                )
                rows.append(
                    compare(case, method, "hessian_numerical", reference, optimized, nat, ax_numerical, ref_hessian)
                )
    return rows


def totals(ncases, rows):
    counts = Counter(row["status"] for row in rows)
    return f"{ncases} cases, {len(rows)} checks: " + ", ".join(f"{counts[status]} {status}" for status in STATUSES)


def cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ").strip()


def numbers(values):
    return ", ".join(f"{key} {value:.6g}" if value is not None else f"{key} —" for key, value in values.items()) or "—"


def summary(report):
    setup = report["metadata"]
    host = setup["host"]
    rows = report["checks"]
    lines = [
        "# ReShiki reference validation",
        "",
        f"Run date (UTC) `{setup['date']}`",
        f"Host `{host['platform']}` / `{host['machine']}` / `{host['cpu']}`",
        "",
    ]
    if report.get("setup_error"):
        lines += [f"Setup failed. {report['setup_error']}", ""]
    else:
        lines += [
            f"axtb checkout revision `{setup['git_revision']}`",
            "",
            "| Program | Version | Binary SHA256 |",
            "|---|---|---|",
        ]
        lines += [
            f"| {name} | {cell(info['version'])} | `{info.get('sha256', '—')}` |"
            for name, info in setup["programs"].items()
        ]
        lines += ["", "Settings", ""] + [f"- **{key}** {value}" for key, value in report["settings"].items()]

    lines += [
        "",
        "## Checks",
        "",
        "| Method | Check | Cases | Pass | Fail | Other | Worst deviation by quantity | Tolerance by quantity |",
        "|---|---|---:|---:|---:|---:|---|---|",
    ]
    groups = defaultdict(list)
    for row in rows:
        groups[row["method"], row["check"]].append(row)
    for (method, check), group in groups.items():
        counts = Counter(row["status"] for row in group)
        worst, limits = {}, defaultdict(set)
        for row in group:
            for key, deviation in row["deviations"].items():
                if deviation is not None:
                    worst[key] = max(worst.get(key, 0.0), deviation)
            for key, gate in row["tolerances"].items():
                limits[key].add(gate["limit"])
        limit_text = (
            ", ".join(
                f"{key} " + ("not gated" if None in values else " / ".join(f"{v:.6g}" for v in sorted(values)))
                for key, values in limits.items()
            )
            or "—"
        )
        lines.append(
            f"| {method} | {check} | {len(group)} | {counts['pass']} | {counts['fail']} | "
            f"{len(group) - counts['pass'] - counts['fail']} | {numbers(worst)} | {limit_text} |"
        )

    lines += [
        "",
        "## Non-pass cases",
        "",
        "| Case | Method | Check | Status | Deviation | Details |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row["status"] != "pass":
            details = [row.get("note", "")]
            details.extend(
                f"{side} {err['kind']}: {err['message'][:240]}" for side, err in row.get("errors", {}).items()
            )
            lines.append(
                f"| {row['case']} | {row['method']} | {row['check']} | {row['status']} | "
                f"{numbers(row['deviations'])} | {cell(' / '.join(filter(None, details)))} |"
            )

    lines += [
        "",
        "## Every check",
        "",
        "Deviation is max |axtb − reference| per quantity; tolerance is the documented limit (— not gated).",
        "",
        "| Case | Method | Check | Reference | Deviation | Tolerance | Iterations axtb / ref | Status |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in sorted(rows, key=lambda r: (r["case"].casefold(), r["method"], r["check"])):
        limits = {q: gate["limit"] for q, gate in row["tolerances"].items()}
        it = row["iterations"]
        deviation = numbers(row["deviations"])
        if dropped := row.get("axtb_modes_dropped"):
            deviation += f" ({dropped} axtb mode(s) omitted from g98.out)"
        lines.append(
            f"| {row['case']} | {row['method']} | {row['check']} | {row['reference'] or '—'} | "
            f"{deviation} | {numbers(limits)} | {it['axtb'] or '—'} / {it['reference'] or '—'} | "
            f"{row['status']} |"
        )
    return "\n".join(lines + ["", totals(report["case_count"], rows), ""])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", help="case ID substring, case insensitive")
    parser.add_argument("--methods", default=",".join(METHODS), help="comma-separated subset of " + ",".join(METHODS))
    parser.add_argument("--jobs", type=int, default=1, help="parallel case processes")
    parser.add_argument("--out", type=Path, default=ROOT / "reshiki/results")
    parser.add_argument("--summarize", type=Path, help="only rewrite summary.md from this results.json")
    args = parser.parse_args()
    if args.summarize:
        args.out.mkdir(parents=True, exist_ok=True)
        (args.out / "summary.md").write_text(summary(json.loads(args.summarize.read_text())))
        return 0
    methods = tuple(dict.fromkeys(method.strip() for method in args.methods.split(",")))
    if not set(methods) <= set(METHODS) or args.jobs < 1:
        parser.error("choose valid --methods and a positive --jobs")
    cases = load_cases(args.only)
    if not cases:
        parser.error(f"no case IDs match {args.only!r}")
    metadata = {"date": datetime.now(timezone.utc).isoformat(), "host": host_metadata()}
    report = {"metadata": metadata, "case_count": len(cases), "checks": []}
    try:
        for method, check in REQUIRED_QUANTITIES:
            if method in methods:
                tolerances(cases[0], method, check)
        binaries, programs = reference_programs()
        binaries["axtb"] = axtb_binary()
        metadata["programs"] = {"axtb": axtb_metadata(binaries["axtb"]), **programs}
        for info in metadata["programs"].values():  # no home directories in published results
            info["path"] = Path(info["path"]).name
        metadata["git_revision"] = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as exc:
        report["setup_error"] = str(exc)
        print(f"Setup failed. {exc}", file=sys.stderr)
    else:
        report["settings"] = {
            "Selection": f"methods {', '.join(methods)}, case substring {args.only!r}, jobs {args.jobs}",
            "Threads": "one per calculation, axtb SCF eigensolver faer (AXTB_LAPACK and AXTB_LAPACK_BLAS unset)",
            "xtb and g-xTB environment": "XTBPATH unset and XTBHOME set to each run's scratch directory",
            "GFN2": "values axtb and tblite 1e-6, xtb 1e-4. Molden both 1e-4. Analytic Hessian 1e-6, xtb 0.001.",
            "g-xTB": "values both 1e-4 and 250 iterations maximum. Molden default accuracy. Reference Hessian 1e-4.",
            "GFN-FF": "xtb 6.8.0 conventions in axtb. Reference gradient and Hessian use xtb defaults.",
            "Long molecules": "cutoff energy and gradient differences compared at their published precision",
            "Hessian geometries": "method-specific optimized XYZ when present, otherwise original XYZ",
            "GFN2 analytic frequencies": "matrix comparison uses modes with reference |ν| > 50 cm⁻¹",
            "Units": "energy and orbital energies Eh, gradient Eh/bohr, charges e, dipole e·bohr, frequency cm⁻¹",
        }
        jobs = [(case, methods, binaries) for case in cases]
        if args.jobs == 1:
            report["checks"] = [row for job in jobs for row in check_case(job)]
        else:
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                report["checks"] = [row for batch in pool.map(check_case, jobs) for row in batch]

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "results.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    (args.out / "summary.md").write_text(summary(report))
    print(totals(len(cases), report["checks"]))
    return int("setup_error" in report or any(row["status"] in ("fail", "error") for row in report["checks"]))


if __name__ == "__main__":
    raise SystemExit(main())
