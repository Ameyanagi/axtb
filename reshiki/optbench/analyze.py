# SPDX-License-Identifier: LGPL-3.0-or-later
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Analyze optimizer benchmark JSONL and write CSV, Markdown, and PNG reports.

Example::

    uv run analyze.py results.jsonl --refs references.jsonl --xtb xtb.jsonl --out report

Geometry paths are relative to their JSONL file. Summary distributions use
converged runs. Failures appear in the counts and per-run CSV. Ratios divide the
candidate by a converged run with the same case, set, start, and potential.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from itertools import pairwise
from pathlib import Path
from typing import Any

import numpy as np

IDENTITY = ("case", "set", "start", "potential")
MEASUREMENTS = (
    "steps",
    "evaluations",
    "hessian_evaluations",
    "potential_seconds",
    "hessian_seconds",
    "optimizer_seconds",
    "wall_seconds",
    "scf_iterations_total",
    "energy",
    "max_gradient",
)
COUNTS = {"steps", "evaluations", "hessian_evaluations", "scf_iterations_total", "gradient_evaluations"}
OPTIONAL_MEASUREMENTS = ("gradient_seconds", "gradient_evaluations", "hessian_gradient_equivalents")
XTB_BASELINE = "xtb-opt-normal"
Identity = tuple[str, str, str, str]
Run = dict[str, Any]


def identity(row: Run) -> Identity:
    return tuple(str(row[name]) for name in IDENTITY)


def finite_number(value: Any, name: str, *, nullable: bool = True, nonnegative: bool = False) -> None:
    if value is None and nullable:
        return
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number{' or null' if nullable else ''}")
    if nonnegative and value < 0:
        raise ValueError(f"{name} must not be negative")
    if name in COUNTS and int(value) != value:
        raise ValueError(f"{name} must be an integer or null")


def validate_identity(row: Run) -> None:
    for name in IDENTITY:
        value = row.get(name)
        if name == "start" and isinstance(value, int) and not isinstance(value, bool) and value >= 0:
            continue
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be a nonempty string (start may also be a nonnegative integer)")


def validate_run(row: Run) -> None:
    validate_identity(row)
    if not isinstance(row.get("config"), str) or not row["config"].strip():
        raise ValueError("config must be a nonempty string")
    atoms = row.get("atoms")
    if atoms is None and row.get("status") == "input_error":
        pass
    elif isinstance(atoms, bool) or not isinstance(atoms, int) or atoms <= 0:
        raise ValueError("atoms must be a positive integer (null for input_error)")
    if not isinstance(row.get("status"), str) or not re.fullmatch(r"[a-z][a-z0-9]*(?:_[a-z0-9]+)*", row["status"]):
        raise ValueError("status must use lower snake_case")
    for name in MEASUREMENTS:
        if name not in row:
            raise ValueError(f"missing measurement {name} (use null when unavailable)")
        finite_number(row[name], name, nonnegative=name != "energy")
    for name in OPTIONAL_MEASUREMENTS:
        if name in row:
            finite_number(row[name], name, nonnegative=True)
    if "final_xyz" not in row or (row["final_xyz"] is not None and not isinstance(row["final_xyz"], str)):
        raise ValueError("final_xyz must be a path string or null")
    if row["final_xyz"] == "":
        raise ValueError("final_xyz must not be an empty path")
    trace = row.get("trace")
    if not isinstance(trace, list):
        raise TypeError("trace must be a list")
    previous = (-1.0, -1.0)
    for point in trace:
        if not isinstance(point, list) or len(point) not in (5, 6):
            raise ValueError(
                "each trace sample must be [evaluations, seconds, energy, max_gradient, accuracy, gradient_norm?]"
            )
        for index, value in enumerate(point[:4]):
            finite_number(value, f"trace[{index}]", nullable=index == 3 and len(point) == 6, nonnegative=index != 2)
        accuracy = point[4]
        if not isinstance(accuracy, str) or accuracy not in ("coarse", "normal", "tight"):
            finite_number(accuracy, "trace accuracy", nonnegative=True)
        if len(point) == 6:
            finite_number(point[5], "trace gradient_norm", nonnegative=True)
        if int(point[0]) != point[0]:
            raise ValueError("trace evaluation counts must be integers")
        if point[0] < previous[0] or point[1] < previous[1]:
            raise ValueError("trace evaluation counts and times must be nondecreasing")
        previous = (point[0], point[1])


def read_jsonl(paths: list[Path], *, reference: bool = False) -> list[Run]:
    rows, seen = [], set()
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                location = f"{path}:{line_number}"
                try:
                    row = json.loads(line, parse_constant=_invalid_constant, parse_float=_finite_float)
                    if not isinstance(row, dict):
                        raise TypeError("each JSONL record must be an object")
                    if reference:
                        validate_identity(row)
                        finite_number(row.get("energy"), "energy", nullable=False)
                        if not isinstance(row.get("xyz"), str) or not row["xyz"]:
                            raise ValueError("reference xyz must be a nonempty path string")
                    else:
                        validate_run(row)
                    key = identity(row) if reference else (*identity(row), row["config"])
                    if key in seen:
                        raise ValueError(f"duplicate {'reference' if reference else 'run'} identity {key!r}")
                    seen.add(key)
                except (ValueError, TypeError) as exc:
                    raise ValueError(f"{location}: {exc}") from exc
                row["_source"] = path.resolve()
                rows.append(row)
    return rows


def _invalid_constant(value: str) -> None:
    raise ValueError(f"nonfinite JSON number {value}")


def _finite_float(value: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"nonfinite JSON number {value}")
    return number


def geometry_path(row: Run, name: str) -> Path | None:
    value = row.get(name)
    if value is None:
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (row["_source"].parent / path).resolve()


def read_xyz(path: Path) -> tuple[list[str], np.ndarray]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        atoms = int(lines[0].strip())
        if atoms <= 0 or len(lines) < atoms + 2 or any(line.strip() for line in lines[atoms + 2 :]):
            raise ValueError("expected one complete XYZ frame")
        symbols, positions = [], []
        for line in lines[2 : atoms + 2]:
            fields = line.split()
            if len(fields) < 4 or not re.fullmatch(r"[A-Za-z]{1,3}", fields[0]):
                raise ValueError("invalid atom line")
            xyz = [float(value) for value in fields[1:4]]
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError("nonfinite coordinate")
            symbols.append(fields[0].capitalize())
            positions.append(xyz)
        return symbols, np.asarray(positions)
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: invalid XYZ: {exc}") from exc


def heavy_atom_rmsd(left: tuple[list[str], np.ndarray], right: tuple[list[str], np.ndarray]) -> float:
    """Kabsch RMSD in Angstrom, without reflections or atom permutation.

    For hydrogen-only molecules use every atom so the comparison remains defined.
    """
    if left[0] != right[0]:
        raise ValueError("run and reference geometries must have the same ordered element symbols")
    indices = [i for i, symbol in enumerate(left[0]) if symbol not in {"H", "D", "T"}]
    if not indices:
        indices = list(range(len(left[0])))
    a, b = left[1][indices], right[1][indices]
    a, b = a - a.mean(axis=0), b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
    residual = a @ (u @ correction @ vt) - b
    return float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))


def positive_ratio(numerator: Any, denominator: Any) -> float | None:
    if numerator is None or denominator is None or numerator <= 0 or denominator <= 0:
        return None
    return float(numerator / denominator)


def hessian_gradient_equivalents(row: Run) -> float | None:
    direct = row.get("hessian_gradient_equivalents")
    if direct is not None:
        return float(direct)
    hessian_seconds = row["hessian_seconds"]
    if hessian_seconds is None:
        return None
    if hessian_seconds == 0:
        return 0.0
    gradient_seconds, gradient_evaluations = row.get("gradient_seconds"), row.get("gradient_evaluations")
    if gradient_seconds is None or gradient_evaluations is None:
        gradient_seconds, gradient_evaluations = row["potential_seconds"], row["evaluations"]
    gradient_cost = positive_ratio(gradient_seconds, gradient_evaluations)
    return hessian_seconds / gradient_cost if gradient_cost is not None else None


def xtb_criterion_evaluations(trace: list[list[Any]]) -> int | None:
    """Find the first trace pair that meets the xTB normal stopping rule."""
    for previous, current in pairwise(trace):
        if len(current) < 6 or current[5] is None or current[4] != previous[4]:
            continue
        if abs(current[2] - previous[2]) < 5e-6 and current[5] < 1e-3:
            return int(current[0])
    return None


def enrich_runs(rows: list[Run], references: list[Run], baseline: str) -> None:
    reference_index = {identity(row): row for row in references}
    run_index = {(*identity(row), row["config"]): row for row in rows}
    geometry_cache = {}

    def geometry(row: Run, name: str) -> tuple[list[str], np.ndarray] | None:
        path = geometry_path(row, name)
        if path is None:
            return None
        if path not in geometry_cache:
            geometry_cache[path] = read_xyz(path)
        return geometry_cache[path]

    for reference in references:
        geometry(reference, "xyz")
    for row in rows:
        final_geometry = geometry(row, "final_xyz")
        if final_geometry is not None and len(final_geometry[0]) != row["atoms"]:
            raise ValueError(f"{row['_source']}: atoms disagrees with final_xyz for {identity(row)!r}")
        reference = reference_index.get(identity(row))
        delta, rmsd, same = None, None, None
        if reference is not None:
            if row["energy"] is not None:
                delta = row["energy"] - reference["energy"]
            if final_geometry is not None:
                rmsd = heavy_atom_rmsd(final_geometry, geometry(reference, "xyz"))
            if delta is not None or rmsd is not None:
                same = (delta is not None and abs(delta) < 1e-6) or (rmsd is not None and rmsd < 0.1)
        row["_reference_energy"] = reference["energy"] if reference is not None else None
        row["reference_status"] = reference.get("status") if reference is not None else None
        row["delta_energy_hartree"] = delta
        row["heavy_atom_rmsd_angstrom"] = rmsd
        row["same_minimum"] = same
        overhead = positive_ratio(row["optimizer_seconds"], row["evaluations"])
        if row["optimizer_seconds"] == 0 and row["evaluations"] is not None and row["evaluations"] > 0:
            overhead = 0.0
        row["optimizer_ms_per_evaluation"] = 1000 * overhead if overhead is not None else None
        row["hessian_gradient_equivalents"] = hessian_gradient_equivalents(row)
        row["xtb_criterion_evaluations"] = xtb_criterion_evaluations(row["trace"])
        for label, config in (("xtb", XTB_BASELINE), ("baseline", baseline)):
            comparator = run_index.get((*identity(row), config))
            paired = comparator is not None and row["status"] == comparator["status"] == "converged"
            for metric in ("evaluations", "wall_seconds"):
                row[f"{label}_{metric}_ratio"] = positive_ratio(row[metric], comparator[metric]) if paired else None
            if label == "xtb":
                row["xtb_criterion_evaluations_ratio"] = (
                    positive_ratio(row["xtb_criterion_evaluations"], comparator["evaluations"]) if paired else None
                )

    for row in rows:
        comparator = run_index.get((*identity(row), baseline))
        if comparator is None:
            continue
        converged = row["status"] == "converged"
        baseline_converged = comparator["status"] == "converged"
        row["baseline_converged_difference"] = int(converged) - int(baseline_converged)
        if converged and baseline_converged:
            row["baseline_overhead_ratio"] = positive_ratio(
                row["optimizer_ms_per_evaluation"], comparator["optimizer_ms_per_evaluation"]
            )
            for metric in (
                "hessian_evaluations",
                "hessian_gradient_equivalents",
                "scf_iterations_total",
                "delta_energy_hartree",
                "same_minimum",
            ):
                if row[metric] is not None and comparator[metric] is not None:
                    row[f"baseline_{metric}_difference"] = row[metric] - comparator[metric]


def values(rows: list[Run], name: str) -> list[float]:
    return [float(row[name]) for row in rows if row.get(name) is not None]


def median(rows: list[Run], name: str) -> float | None:
    selected = values(rows, name)
    return statistics.median(selected) if selected else None


def mean(rows: list[Run], name: str) -> float | None:
    selected = values(rows, name)
    return statistics.mean(selected) if selected else None


def ratio_statistics(rows: list[Run], name: str) -> tuple[float | None, float | None, int]:
    selected = values(rows, name)
    if not selected:
        return None, None, 0
    return statistics.median(selected), math.exp(statistics.mean(math.log(value) for value in selected)), len(selected)


def fmt(value: float | None) -> str:
    return "—" if value is None else f"{value:.6g}"


def escape(value: Any) -> str:
    return str(value).replace("|", "\\|").replace("\n", " ")


def paired_cell(rows: list[Run], name: str) -> str:
    med, geometric, count = ratio_statistics(rows, name)
    return f"{fmt(med)} / {fmt(geometric)} ({count})"


def grouped_runs(rows: list[Run]) -> dict[tuple[str, str, str], list[Run]]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["potential"], row["set"], row["config"])].append(row)
    return dict(sorted(grouped.items()))


def write_summary(path: Path, rows: list[Run], baseline: str, molecules: list[Identity]) -> None:
    lines = [
        "# Optimizer benchmark analysis",
        "",
        (
            "Counts include every run. All other aggregate statistics use converged runs with the relevant measurement. "
            "A dash indicates an unavailable measurement. Evals count energy and gradient calls. Times are seconds."
        ),
        "",
        (
            "ΔE is final energy minus reference energy in Eh. Same minimum means "
            "Kabsch heavy-atom RMSD < 0.1 Å **or** |ΔE| < 1e−6 Eh. The fraction denominator includes converged "
            "runs with a comparable reference. Kabsch preserves atom order and excludes reflections. Hydrogen-only "
            "molecules use all atoms."
        ),
        "",
        (
            "Hessian cost is measured gradient equivalents when supplied, otherwise Hessian seconds divided by "
            "measured gradient seconds per call. Without separate gradient timing, potential seconds/evaluations "
            "is the estimate. Optimizer overhead is optimizer seconds/evaluations × 1000."
        ),
        "",
        (
            "Ratios are candidate / comparator. Each cell gives median / geometric mean (paired count). "
            "Only matching case × set × start × potential pairs with both runs converged and positive measurements "
            f"contribute. The official comparator is `{escape(XTB_BASELINE)}`."
        ),
        "",
        (
            "xTB criterion evals use the first pair of consecutive trace samples at the same accuracy with "
            "absolute energy change < 5e−6 Eh and gradient norm < 1e−3 Eh/bohr at the second sample. "
            "Both evaluation ratios to xTB divide by the official run's recorded evaluations (cycles in new files)."
        ),
        "",
        (
            "| Potential | Set | Config | Runs | Converged | Failures | Evals median / mean | "
            "xTB criterion evals median / mean | Wall s median / mean | "
            "Overhead ms/eval median | Hessian evals median | Hessian grad-equiv median | SCF iters median | ΔE median Eh | "
            "ΔE max Eh | Same minimum | Evals ratio vs xTB median / GM (n) | "
            "Criterion evals ratio vs xTB median / GM (n) | Wall ratio vs xTB median / GM (n) |"
        ),
        "|---|---|---|" + "---:|" * 16,
    ]
    groups = grouped_runs(rows)
    for key, group in groups.items():
        converged = [row for row in group if row["status"] == "converged"]
        same = [row["same_minimum"] for row in converged if row["same_minimum"] is not None]
        deltas = values(converged, "delta_energy_hartree")
        cells = [
            *map(escape, key),
            str(len(group)),
            str(len(converged)),
            str(len(group) - len(converged)),
            f"{fmt(median(converged, 'evaluations'))} / {fmt(mean(converged, 'evaluations'))}",
            f"{fmt(median(converged, 'xtb_criterion_evaluations'))} / {fmt(mean(converged, 'xtb_criterion_evaluations'))}",
            f"{fmt(median(converged, 'wall_seconds'))} / {fmt(mean(converged, 'wall_seconds'))}",
            fmt(median(converged, "optimizer_ms_per_evaluation")),
            fmt(median(converged, "hessian_evaluations")),
            fmt(median(converged, "hessian_gradient_equivalents")),
            fmt(median(converged, "scf_iterations_total")),
            fmt(statistics.median(deltas) if deltas else None),
            fmt(max(deltas) if deltas else None),
            f"{sum(same)}/{len(same)} ({sum(same) / len(same):.1%})" if same else "—",
            paired_cell(converged, "xtb_evaluations_ratio"),
            paired_cell(converged, "xtb_criterion_evaluations_ratio"),
            paired_cell(converged, "xtb_wall_seconds_ratio"),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += [
        "",
        f"## Paired metrics relative to `{escape(baseline)}`",
        "",
        (
            "Differences are candidate minus baseline. The converged-count difference uses every matched pair. "
            "Other comparisons use pairs where both converged. Same-minimum differences are in percentage points."
        ),
        "",
        (
            "| Potential | Set | Config | Matched | Converged diff | Evals ratio median / GM (n) | "
            "Wall ratio median / GM (n) | Overhead ratio median / GM (n) | Hessian evals diff median | "
            "Hessian grad-equiv diff median | SCF iters diff median | ΔE diff median / max Eh | Same-minimum diff pp |"
        ),
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key, group in groups.items():
        statuses = values(group, "baseline_converged_difference")
        energy_differences = values(group, "baseline_delta_energy_hartree_difference")
        same_differences = mean(group, "baseline_same_minimum_difference")
        cells = [
            *map(escape, key),
            str(len(statuses)),
            fmt(sum(statuses) if statuses else None),
            paired_cell(group, "baseline_evaluations_ratio"),
            paired_cell(group, "baseline_wall_seconds_ratio"),
            paired_cell(group, "baseline_overhead_ratio"),
            fmt(median(group, "baseline_hessian_evaluations_difference")),
            fmt(median(group, "baseline_hessian_gradient_equivalents_difference")),
            fmt(median(group, "baseline_scf_iterations_total_difference")),
            f"{fmt(statistics.median(energy_differences))} / {fmt(max(energy_differences))}"
            if energy_differences
            else "—",
            fmt(same_differences * 100 if same_differences is not None else None),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines += [
        "",
        "## Plots",
        "",
        "- `evaluations.png`: converged evaluations per molecule/start, ordered by atom count; logarithmic y axis.",
        (
            "- `convergence.png`: energy above the reference (Eh) versus cumulative evaluations and wall time. "
            "The logarithmic axis displays gaps at or below 1e−12 Eh at that floor. "
            "Preoptimization offsets are preserved from the recorded trace. Traces with timestamps unavailable "
            "appear only in the evaluation panel."
        ),
        "- `evaluation_cdf.png`: empirical CDF of converged evaluations. Missing values and failures are omitted.",
        "",
        "Trace selections: " + (", ".join("`" + escape("/".join(key)) + "`" for key in molecules) or "none"),
        "",
        (
            "`per_molecule.csv` contains every run, its original fields, source JSONL, and derived comparisons. "
            "Nested fields and traces are JSON-encoded. Unavailable values are empty CSV cells."
        ),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def write_csv(path: Path, rows: list[Run]) -> None:
    first = [*IDENTITY, "config", "atoms", "status", *MEASUREMENTS]
    fields = first + sorted(
        {name for row in rows for name in row if not name.startswith("_")} - set(first) - {"source_jsonl"}
    )
    fields.append("source_jsonl")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            output = {name: value for name, value in row.items() if not name.startswith("_")}
            output["source_jsonl"] = str(row["_source"])
            for name, value in output.items():
                if isinstance(value, (dict, list)):
                    output[name] = json.dumps(value, separators=(",", ":"), allow_nan=False)
            writer.writerow(output)


def choose_configs(rows: list[Run], requested: str | None) -> list[str]:
    available = {row["config"] for row in rows}
    selected = list(dict.fromkeys(requested.split(","))) if requested else sorted(available)
    unknown = set(selected) - available
    if unknown:
        raise ValueError(f"unknown configurations: {', '.join(sorted(unknown))}")
    if XTB_BASELINE in available and XTB_BASELINE not in selected:
        selected.append(XTB_BASELINE)
    return selected


def choose_molecules(rows: list[Run], requested: str | None) -> list[Identity]:
    available = {identity(row): row["atoms"] for row in rows if row["trace"]}
    ordered = sorted(available, key=lambda key: (available[key], key))
    if requested:
        selected = []
        for selector in requested.split(","):
            matches = [key for key in ordered if selector == key[0] or selector == "/".join(key)]
            if not matches:
                raise ValueError(f"no trace data for molecule selector {selector!r}")
            if matches[0] not in selected:
                selected.append(matches[0])
        if len(selected) > 4:
            raise ValueError("choose at most four molecules")
        return selected
    by_case = {}
    for key in ordered:
        by_case.setdefault(key[0], key)
    candidates = list(by_case.values())
    indices = np.linspace(0, len(candidates) - 1, min(4, len(candidates)), dtype=int)
    return [candidates[index] for index in indices]


def energy_gap_trace(row: Run) -> list[tuple[int, float, float]]:
    reference = row["_reference_energy"]
    if reference is None:
        return []
    return [(point[0], point[1], max(1e-12, point[2] - reference)) for point in row["trace"]]


def write_plots(out: Path, rows: list[Run], configs: list[str], molecules: list[Identity]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9})
    palette = plt.get_cmap("tab20").colors
    palette = palette[::2] + palette[1::2]
    colors = {config: palette[index % len(palette)] for index, config in enumerate(configs)}
    facets = sorted({(row["potential"], row["set"]) for row in rows})
    ncols = min(2, len(facets))
    nrows = math.ceil(len(facets) / ncols)
    legend = [Line2D([0], [0], color=colors[config], marker="o", label=config) for config in configs]

    for kind in ("evaluations", "evaluation_cdf"):
        figure, axes = plt.subplots(nrows, ncols, figsize=(max(10, ncols * 8), nrows * 4.7), squeeze=False)
        for axis, (potential, dataset) in zip(axes.flat, facets):
            group = [row for row in rows if (row["potential"], row["set"]) == (potential, dataset)]
            sizes = {identity(row): row["atoms"] for row in group if row["atoms"] is not None}
            keys = sorted(sizes, key=lambda key: (sizes[key], key))
            indices = {key: index for index, key in enumerate(keys)}
            for config in configs:
                selected = [
                    row
                    for row in group
                    if row["config"] == config
                    and row["status"] == "converged"
                    and row["evaluations"] is not None
                    and row["evaluations"] > 0
                ]
                selected.sort(key=lambda row: indices[identity(row)])
                if kind == "evaluations":
                    axis.plot(
                        [indices[identity(row)] for row in selected],
                        [row["evaluations"] for row in selected],
                        marker="o",
                        markersize=4,
                        linewidth=1.2,
                        color=colors[config],
                        label=config,
                    )
                elif selected:
                    data, counts = np.unique([row["evaluations"] for row in selected], return_counts=True)
                    fraction = np.cumsum(counts) / len(selected)
                    axis.step(np.r_[data[0], data], np.r_[0, fraction], where="post", color=colors[config])
                    axis.scatter(data, fraction, s=10, color=colors[config])
            axis.set_title(f"{potential} · {dataset}")
            axis.grid(True, alpha=0.25)
            if kind == "evaluations":
                axis.set_yscale("log")
                stride = max(1, math.ceil(len(keys) / 35))
                shown = list(range(0, len(keys), stride))
                axis.set_xticks(
                    shown,
                    [f"{keys[index][0]}/{keys[index][2]}" for index in shown],
                    rotation=55,
                    ha="right",
                    fontsize=8,
                )
                axis.set_xlabel("Molecule / start (increasing atom count)")
                axis.set_ylabel("Potential evaluations")
            else:
                axis.set_ylim(0, 1.03)
                axis.set_xlabel("Potential evaluations")
                axis.set_ylabel("Fraction of converged runs")
        for axis in list(axes.flat)[len(facets) :]:
            axis.set_visible(False)
        figure.legend(handles=legend, loc="upper center", ncol=min(5, len(configs)), frameon=False)
        figure.tight_layout(rect=(0, 0, 1, 1 - 0.035 * math.ceil(len(configs) / 5)))
        figure.savefig(out / f"{kind}.png", dpi=150)
        plt.close(figure)

    figure, axes = plt.subplots(max(1, len(molecules)), 2, figsize=(14, 3.6 * max(1, len(molecules))), squeeze=False)
    if not molecules:
        for axis in axes.flat:
            axis.text(0.5, 0.5, "No energy traces with matching references", transform=axis.transAxes, ha="center")
            axis.set_axis_off()
    for index, key in enumerate(molecules):
        selected = [row for row in rows if identity(row) == key and row["config"] in configs]
        for coordinate, axis in enumerate(axes[index]):
            for row in selected:
                trace = energy_gap_trace(row)
                if trace and (coordinate == 0 or any(point[1] > 0 for point in trace)):
                    axis.plot(
                        [point[coordinate] for point in trace],
                        [point[2] for point in trace],
                        color=colors[row["config"]],
                        linewidth=1.4,
                        label=row["config"],
                    )
            axis.set_yscale("log")
            axis.set_xlabel("Cumulative potential evaluations" if coordinate == 0 else "Cumulative wall time (s)")
            axis.set_ylabel("E − reference energy (Eh)")
            axis.set_title(" · ".join(key))
            axis.grid(True, alpha=0.25)
    figure.legend(handles=legend, loc="upper center", ncol=min(5, len(configs)), frameon=False)
    figure.tight_layout(rect=(0, 0, 1, 1 - 0.035 * math.ceil(len(configs) / 5)))
    figure.savefig(out / "convergence.png", dpi=150)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", type=Path, nargs="+", help="optimizer run JSONL files")
    parser.add_argument("--refs", type=Path, nargs="+", default=[], help="reference JSONL files")
    parser.add_argument("--xtb", type=Path, nargs="+", default=[], help="official xTB run JSONL files")
    parser.add_argument("--out", type=Path, required=True, help="report output directory")
    parser.add_argument("--baseline", default="anc", help="paired baseline config (default: anc)")
    parser.add_argument("--configs", help="comma-separated plot configs (default: all; xTB normal always included)")
    parser.add_argument(
        "--molecules", help="up to four comma-separated cases or case/set/start/potential trace selectors"
    )
    args = parser.parse_args(argv)
    try:
        rows = read_jsonl([*args.results, *args.xtb])
        if not rows:
            raise ValueError("no run records found")
        references = read_jsonl(args.refs, reference=True)
        enrich_runs(rows, references, args.baseline)
        rows.sort(
            key=lambda row: (
                row["potential"],
                row["set"],
                row["atoms"] or math.inf,
                row["case"],
                str(row["start"]),
                row["config"],
            )
        )
        configs = choose_configs(rows, args.configs)
        molecules = choose_molecules(
            [row for row in rows if row["config"] in configs and row["_reference_energy"] is not None], args.molecules
        )
        args.out.mkdir(parents=True, exist_ok=True)
        write_csv(args.out / "per_molecule.csv", rows)
        write_summary(args.out / "summary.md", rows, args.baseline, molecules)
        write_plots(args.out, rows, configs, molecules)
    except (ValueError, TypeError, OSError) as exc:
        print(f"analyze.py: error: {exc}", file=sys.stderr)
        return 2
    print(f"Analyzed {len(rows)} runs; report written to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
