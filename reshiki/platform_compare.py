# SPDX-License-Identifier: LGPL-3.0-or-later
"""Compare another platform's axtb dump against a Linux x64 baseline.

Usage: ``python reshiki/platform_compare.py BASE.json OTHER.json [--md report.md]``.
Numerical gates are |dE| <= 1e-8 Eh, max|dG| <= 1e-7 Eh/bohr, and
max sorted |dnu| <= 0.1 cm^-1. Charges and iteration counts are reported without
gates. Input fingerprints must agree. Declared expected error pairs can have
different diagnostic text across operating systems.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from cases import EXPECTED_ERRORS

# The ReShiki lead chose these cross-platform tolerances.
ENERGY_LIMIT = 1e-8
GRADIENT_LIMIT = 1e-7
FREQUENCY_LIMIT = 0.1
Key = tuple[str, str, str]


def record_key(record: dict[str, Any]) -> Key:
    return tuple(str(record[name]) for name in ("case", "method", "run"))


def index_records(dump: dict[str, Any]) -> dict[Key, dict[str, Any]]:
    """Reject duplicate identities instead of silently losing a run."""
    if dump.get("schema_version") != 1:
        raise ValueError(f"unsupported schema_version: {dump.get('schema_version')!r}")
    records = dump.get("records")
    if not isinstance(records, list):
        raise TypeError("dump must contain a records list")
    indexed = {}
    for record in records:
        key = record_key(record)
        if key in indexed:
            raise ValueError(f"duplicate record key: {key!r}")
        if ("result" in record) == ("error_kind" in record):
            raise ValueError(f"record must contain exactly one of result or error_kind: {key!r}")
        if "result" in record and not isinstance(record["result"], dict):
            raise ValueError(f"result must be an object: {key!r}")
        if "error_kind" in record and (
            not isinstance(record["error_kind"], str) or not record["error_kind"] or "error" not in record
        ):
            raise ValueError(f"error record has no stable kind or diagnostic: {key!r}")
        indexed[key] = record
    return indexed


def selected_keys(dump: dict[str, Any]) -> set[Key]:
    """Expand explicitly scheduled case/method/run identities, including lost records."""
    selections = dump.get("selections", [])
    if not isinstance(selections, list):
        raise TypeError("dump selections must be a list")
    keys = set()
    for selection in selections:
        if not isinstance(selection, dict) or not isinstance(selection.get("case"), str) or not selection["case"]:
            raise ValueError("each dump selection must name a case")
        for field in ("methods", "runs"):
            values = selection.get(field)
            if not isinstance(values, list) or not values or any(not isinstance(v, str) or not v for v in values):
                raise ValueError(f"selection {selection['case']!r} must contain a non-empty {field} list")
        keys.update((selection["case"], method, run) for method in selection["methods"] for run in selection["runs"])
    return keys


def matching_fingerprints(
    left: dict[str, Any] | None, right: dict[str, Any] | None, base_label: str, other_label: str, problems: list[str]
) -> bool:
    """Require valid SHA-256 fingerprints and explain how to replace legacy dumps."""
    fingerprints = []
    for label, record in ((base_label, left), (other_label, right)):
        value = record.get("input_fingerprint") if record is not None else None
        if record is not None:
            if "input_fingerprint" not in record:
                problems.append(
                    f"missing input fingerprint in {label} (regenerate this legacy dump with platform_dump.py)"
                )
            elif value is None:
                problems.append(f"unavailable input fingerprint in {label} (CLI request was not constructed)")
            elif not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                problems.append(f"invalid input fingerprint in {label} (expected a SHA-256 hex digest)")
                value = None
        fingerprints.append(value)
    if all(fingerprints) and fingerprints[0] != fingerprints[1]:
        problems.append("input fingerprints differ")
    return bool(all(fingerprints) and fingerprints[0] == fingerprints[1])


def frequency_count(values: Any) -> int | None:
    """Return the number of modes only for a one-dimensional frequency vector."""
    try:
        array = np.asarray(values)
    except (TypeError, ValueError):
        return None
    return int(array.size) if array.ndim == 1 else None


def max_difference(base: Any, other: Any, *, sort: bool = False, ndim: int | None = None) -> float:
    """Largest absolute difference, or infinity for malformed/non-finite data."""
    if base is None or other is None:
        return math.inf
    try:
        left, right = np.asarray(base), np.asarray(other)
    except (TypeError, ValueError):
        return math.inf
    if left.dtype.kind not in "iuf" or right.dtype.kind not in "iuf":
        return math.inf
    left, right = left.astype(float), right.astype(float)
    if left.shape != right.shape or (ndim is not None and left.ndim != ndim):
        return math.inf
    if not np.isfinite(left).all() or not np.isfinite(right).all():
        return math.inf
    if sort:
        left, right = np.sort(left), np.sort(right)
    if left.size == 0:
        return math.inf
    with np.errstate(over="ignore", invalid="ignore"):
        return float(np.max(np.abs(right - left)))


def metric_text(value: float | None) -> str:
    if value is None:
        return "—"
    return f"{value:.3e}" if math.isfinite(value) else "invalid"


def markdown_text(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ").replace("\r", " ")


def host_description(dump: dict[str, Any]) -> str:
    metadata = dump.get("metadata", {})
    return " / ".join(str(metadata.get(name, "unknown")) for name in ("platform", "machine", "cpu"))


def compare_dumps(
    base: dict[str, Any], other: dict[str, Any], base_label: str = "BASE", other_label: str = "OTHER"
) -> tuple[str, int]:
    """Return a Markdown report and the number of comparison violations."""
    baseline, candidate = index_records(base), index_records(other)
    scheduled = selected_keys(base) | selected_keys(other)
    stats: dict[tuple[str, str], dict[str, Any]] = defaultdict(
        lambda: {
            "runs": 0,
            "compared": 0,
            "errors": 0,
            "failed": 0,
            "iterations": 0,
            "iterations_base": 0,
            "iterations_other": 0,
            "iteration_pairs": 0,
            "energy": None,
            "gradient": None,
            "frequencies": None,
            "charges": None,
        }
    )
    violations: list[str] = []
    per_run: list[str] = []
    iteration_changes: list[str] = []
    error_changes: list[str] = []
    for key in sorted(baseline.keys() | candidate.keys() | scheduled):
        label = "/".join(key)
        group = stats[key[1:]]
        group["runs"] += 1
        left, right = baseline.get(key), candidate.get(key)
        problems = []
        deltas = {}
        fingerprints_match = matching_fingerprints(left, right, base_label, other_label, problems)
        if left is None or right is None:
            missing = [label for label, record in ((base_label, left), (other_label, right)) if record is None]
            problems.append(
                "missing from " + " and ".join(missing) + (" (selected in dump header)" if key in scheduled else "")
            )
        left_kind, right_kind = (record.get("error_kind") if record else None for record in (left, right))
        if left_kind is not None or right_kind is not None:
            expected_kind = EXPECTED_ERRORS.get(key[:2])
            if expected_kind is not None and left_kind == right_kind == expected_kind:
                group["errors"] += 1
                if left["error"] != right["error"]:
                    error_changes.append(f"{label}: diagnostic text differs (kind {left['error_kind']!r}, ungated)")
            else:
                for side, kind in ((base_label, left_kind), (other_label, right_kind)):
                    if kind is not None:
                        problems.append(f"{side} error {kind!r} is not a matching declared expected error")
        elif left is not None and right is not None and fingerprints_match:
            group["compared"] += 1
            left_result, right_result = left["result"], right["result"]
            deltas = {
                "energy": max_difference(left_result.get("energy"), right_result.get("energy"), ndim=0),
                "gradient": max_difference(left_result.get("gradient"), right_result.get("gradient"), ndim=2),
                "charges": max_difference(left_result.get("charges"), right_result.get("charges"), ndim=1),
            }
            try:
                gradient = np.asarray(left_result.get("gradient"))
            except (TypeError, ValueError):
                gradient = np.empty((0, 0))
            if gradient.ndim != 2 or gradient.shape[1] != 3 or gradient.shape[0] == 0:
                deltas["gradient"] = math.inf
            if key[2] == "hessian-analytic":
                left_count = frequency_count(left_result.get("frequencies"))
                right_count = frequency_count(right_result.get("frequencies"))
                if left_count != right_count:
                    problems.append(f"frequency mode counts differ: {left_count!r} -> {right_count!r}")
                deltas["frequencies"] = max_difference(
                    left_result.get("frequencies"), right_result.get("frequencies"), sort=True, ndim=1
                )
            for field, delta in deltas.items():
                old_max = group[field]
                group[field] = delta if old_max is None else max(old_max, delta)
            for field, limit in (
                ("energy", ENERGY_LIMIT),
                ("gradient", GRADIENT_LIMIT),
                ("frequencies", FREQUENCY_LIMIT),
            ):
                if field in deltas and deltas[field] > limit:
                    problems.append(f"{field} delta {metric_text(deltas[field])} exceeds {limit:g}")
            left_iterations, right_iterations = left_result.get("iterations"), right_result.get("iterations")
            if isinstance(left_iterations, int) and isinstance(right_iterations, int):
                group["iterations_base"] += left_iterations
                group["iterations_other"] += right_iterations
                group["iteration_pairs"] += 1
            if left_iterations != right_iterations:
                group["iterations"] += 1
                iteration_changes.append(f"{label}: {left_iterations!r} -> {right_iterations!r}")
        per_run.append(
            f"| {markdown_text(label)} | "
            + " | ".join(metric_text(deltas.get(f)) for f in ("energy", "gradient", "frequencies"))
            + f" | {'fail' if problems else ('expected error' if not deltas else 'pass')} |"
        )
        if problems:
            group["failed"] += 1
            violations.append(f"{label} — " + ", ".join(problems))
    if not any(group["compared"] for group in stats.values()):
        violations.append(
            "No numeric pairs were compared. At least one successful pair with matching input fingerprints is required"
        )
    failed = len(violations)
    lines = [
        f"# axtb platform comparison — {'FAIL' if failed else 'PASS'}",
        "",
        f"BASE (Linux x64 reference): `{markdown_text(base_label)}` — {markdown_text(host_description(base))}",
        f"OTHER: `{markdown_text(other_label)}` — {markdown_text(host_description(other))}",
        "",
        (
            "Limits are 1e-8 Eh for |dE|, 1e-7 Eh/bohr for max|dG|, and 0.1 cm⁻¹ for max sorted |dν|. "
            "Charges and iteration counts are measured without gates. Input fingerprints must match, selected runs "
            "must be present, and at least one numeric pair is required. Only declared expected error pairs are accepted."
        ),
        "",
        (
            "| Method / run | Runs | Numeric pairs | Matching errors | Max ΔE (Eh) | Max ΔG (Eh/bohr) "
            "| Max Δν (cm⁻¹) | Max Δq (e, ungated) | SCF iterations BASE / OTHER | Iteration changes | Failed runs |"
        ),
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    totals = {
        field: sum(group[field] for group in stats.values())
        for field in ("runs", "compared", "errors", "iterations", "iterations_base", "iterations_other", "failed")
    }
    maxima = {
        field: max((group[field] for group in stats.values() if group[field] is not None), default=None)
        for field in ("energy", "gradient", "frequencies", "charges")
    }
    for (method, run), group in sorted(stats.items()):
        iterations = f"{group['iterations_base']} / {group['iterations_other']}" if group["iteration_pairs"] else "—"
        lines.append(
            f"| {markdown_text(method)} / {markdown_text(run)} | {group['runs']} | {group['compared']} "
            f"| {group['errors']} | {metric_text(group['energy'])} | {metric_text(group['gradient'])} "
            f"| {metric_text(group['frequencies'])} | {metric_text(group['charges'])} "
            f"| {iterations} | {group['iterations']} | {group['failed']} |"
        )
    lines.extend(
        [
            (
                f"| **Total** | {totals['runs']} | {totals['compared']} | {totals['errors']} "
                f"| {metric_text(maxima['energy'])} | {metric_text(maxima['gradient'])} "
                f"| {metric_text(maxima['frequencies'])} | {metric_text(maxima['charges'])} "
                f"| {totals['iterations_base']} / {totals['iterations_other']} "
                f"| {totals['iterations']} | {totals['failed']} |"
            ),
            "",
            (f"Violations **{failed}**. Changed diagnostics for expected error pairs: {len(error_changes)} (ungated)."),
        ]
    )
    for title, details in (
        ("Violations", violations),
        ("Iteration changes (ungated)", iteration_changes),
        ("Diagnostic text changes (ungated)", error_changes),
    ):
        if details:
            lines.extend(["", f"**{title}**", ""])
            lines.extend(f"- {markdown_text(detail)}" for detail in details)
    lines += [
        "",
        "**Every run**",
        "",
        "| Case / method / run | ΔE (Eh) | max ΔG (Eh/bohr) | max Δν (cm⁻¹) | Status |",
        "|---|---:|---:|---:|---|",
        *per_run,
    ]
    return "\n".join(lines) + "\n", failed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", metavar="BASE.json", type=Path, help="Linux x64 baseline dump")
    parser.add_argument("other", metavar="OTHER.json", type=Path, help="dump to compare against BASE")
    parser.add_argument("--md", type=Path, metavar="FILE", help="also write the Markdown report to FILE")
    args = parser.parse_args(argv)
    try:
        base = json.loads(args.base.read_text(encoding="utf-8"))
        other = json.loads(args.other.read_text(encoding="utf-8"))
        report, failed = compare_dumps(base, other, args.base.name, args.other.name)
        if args.md:
            args.md.parent.mkdir(parents=True, exist_ok=True)
            args.md.write_text(report, encoding="utf-8")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"platform comparison failed: {exc}", file=sys.stderr)
        return 1
    print(report, end="")
    return int(failed > 0)


if __name__ == "__main__":
    raise SystemExit(main())
