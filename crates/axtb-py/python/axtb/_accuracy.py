# SPDX-License-Identifier: LGPL-3.0-or-later
"""Accuracy targets onto axtb's SCF accuracy factor (design 7; decision 22).

axtb's SCF thresholds all scale with one factor, ``acc`` (``accuracy_scale``): GFN2-xTB's
energy change 1e-6 acc Eh and RMS charge residual 2e-5 acc (xtb 6.7.1: per atom 1e-4 acc),
g-xTB 2.0.1's energy change 1e-6 acc Eh and Fock residual 2e-5 acc, v1.1's energy change 1e-8 acc
(open shells 4e-7 acc) and density change 2e-5 acc (5e-5 acc). A target maps onto ``acc`` through
an error model fitted to measurements (tests/accuracy_calibration.py, bench/results/accuracy.json):
per method, conventions kind (axtb or the reference program's) and shell (closed or open), the
upper envelope err_q <= k_q acc^p_q over molecules, fresh and warm starts, for the energy (Eh) and
the largest gradient component (Eh/bohr). A target's factor is the largest measured factor
(``CALIBRATION["accs"]``, ``ACC_MAX`` down to ``ACC_FLOOR``) whose bound times the safety factor
2 meets every target in every model that applies (the composite: its own fit and its methods';
an open shell: its fit and the closed shells'), so an estimate rests on a measurement at its
factor or above (between measured factors the errors plateau); an explicit ``accuracy_scale``
is a floor: ``acc = min(acc, accuracy_scale)``. g-xTB v1.1's open shells do not fall
monotonically with ``acc``: they use a table of the worst error at each measured factor. A
Hessian target sets ``hessian.accuracy_scale`` the same way (its own model per method and
Hessian method), never looser than the method's default or an explicit value.

``estimate`` is the model's error at the factor used: a bound fitted to measurements, not a
measurement of the run. A target below what the floor reaches runs at the floor, records an
estimate above the target and warns ``accuracy_limited``. GFN-FF has no convergence settings: it
meets every target (no thresholds, estimates 0). GFN2-xTB's integral screening stays at
``accuracy_scale``: a target changes the SCF thresholds only.
"""

from __future__ import annotations

from typing import Any

from qcdata import Accuracy

__all__ = ["ACC_FLOOR", "ACC_MAX", "SAFETY", "hessian_factor", "met", "scf_factor"]

SAFETY = 2.0
ACC_MAX = 10.0
ACC_FLOOR = 1e-3

# BEGIN CALIBRATION (written by tests/accuracy_calibration.py --write; do not edit by hand)
CALIBRATION: dict[str, Any] = {
    "source": "bench/results/accuracy.json",
    "axtb": "0.2.0",
    "commit": None,
    "acc_floor": 0.001,
    "accs": [10.0, 3.0, 1.0, 0.3, 0.1, 0.03, 0.01, 0.001],
    "hessian_accs": [0.1, 0.01, 0.001, 0.0001, 1e-05, 1e-06],
    "scf": {
        "composite axtb closed": {"energy": {"k": 2.6e-07, "p": 1.25}, "gradient": {"k": 4.32e-05, "p": 0.75}},
        "composite axtb open": {"energy": {"k": 6.51e-08, "p": 1.25}, "gradient": {"k": 1.97e-05, "p": 0.75}},
        "gfn2 axtb closed": {"energy": {"k": 3.04e-08, "p": 1.75}, "gradient": {"k": 1.36e-05, "p": 1.0}},
        "gfn2 axtb open": {"energy": {"k": 1.27e-09, "p": 2.0}, "gradient": {"k": 6.87e-06, "p": 0.75}},
        "gfn2 program closed": {"energy": {"k": 3.04e-08, "p": 1.75}, "gradient": {"k": 9.45e-06, "p": 1.0}},
        "gfn2 program open": {"energy": {"k": 1.28e-09, "p": 2.0}, "gradient": {"k": 4.76e-06, "p": 1.0}},
        "gxtb-v1.1 axtb closed": {"energy": {"k": 1.59e-07, "p": 0.75}, "gradient": {"k": 6.46e-05, "p": 0.5}},
        "gxtb-v1.1 program closed": {"energy": {"k": 8.27e-07, "p": 0.75}, "gradient": {"k": 0.000109, "p": 0.5}},
        "gxtb-v2.0.1 axtb closed": {"energy": {"k": 1.21e-07, "p": 1.25}, "gradient": {"k": 4.38e-05, "p": 0.75}},
        "gxtb-v2.0.1 axtb open": {"energy": {"k": 6.34e-08, "p": 1.25}, "gradient": {"k": 2.15e-05, "p": 0.75}},
        "gxtb-v2.0.1 program closed": {"energy": {"k": 1.35e-07, "p": 1.25}, "gradient": {"k": 7.8e-05, "p": 0.75}},
        "gxtb-v2.0.1 program open": {"energy": {"k": 2.19e-07, "p": 1.25}, "gradient": {"k": 3.47e-05, "p": 1.0}},
    },
    "tables": {
        "gxtb-v1.1 axtb open": {
            "energy": [
                [10.0, 8.28e-07],
                [3.0, 4.09e-08],
                [1.0, 3.04e-08],
                [0.3, 2.43e-08],
                [0.1, 8.46e-09],
                [0.03, 8.46e-09],
                [0.01, 8.43e-09],
                [0.001, 5.62e-10],
            ],
            "gradient": [
                [10.0, 0.00043],
                [3.0, 2.01e-05],
                [1.0, 1.21e-05],
                [0.3, 9.83e-06],
                [0.1, 1.8e-06],
                [0.03, 1.8e-06],
                [0.01, 1.71e-06],
                [0.001, 4.39e-07],
            ],
        },
        "gxtb-v1.1 program open": {
            "energy": [
                [10.0, 2.9e-07],
                [3.0, 5.64e-08],
                [1.0, 2.72e-08],
                [0.3, 2.17e-08],
                [0.1, 1.12e-08],
                [0.03, 1.12e-08],
                [0.01, 7.19e-09],
                [0.001, 5.64e-09],
            ],
            "gradient": [
                [10.0, 7.34e-05],
                [3.0, 1.69e-05],
                [1.0, 4.09e-06],
                [0.3, 3.66e-06],
                [0.1, 1.83e-06],
                [0.03, 1.83e-06],
                [0.01, 1.59e-06],
                [0.001, 7.99e-07],
            ],
        },
    },
    "hessian": {
        "gfn2 analytic": {"k": 0.001, "p": 1.0, "floor": 0.01},
        "gfn2 response": {"k": 0.001, "p": 1.0, "floor": 0.01},
        "gxtb-v1.1 numerical": {"k": 0.00345, "p": 0.5, "floor": 1e-05},
        "gxtb-v2.0.1 analytic": {"k": 2.24e-05, "p": 0.5, "floor": 0.001},
        "gxtb-v2.0.1 response": {"k": 2.24e-05, "p": 0.5, "floor": 0.001},
    },
}
# END CALIBRATION


def _kind(method: str, conventions: str) -> str:
    return "axtb" if conventions == "axtb" else "program"


def _thresholds(method: str, conventions: str, open_shell: bool, acc: float) -> dict[str, float]:
    """The SCF's convergence thresholds at a factor (design 7.1)."""
    if method == "gxtb-v1.1":
        return {
            "energy_tol": (4e-7 if open_shell else 1e-8) * acc,
            "residual_tol": (5e-5 if open_shell else 2e-5) * acc,
        }
    if method == "gfn2" and conventions == "xtb-6.7.1":
        return {"energy_tol": 1e-6 * acc, "residual_tol_per_atom": 1e-4 * acc}
    return {"energy_tol": 1e-6 * acc, "residual_tol": 2e-5 * acc}


def _models(method: str, conventions: str, composite: bool, gfn2_conventions: str | None, open_shell: bool):
    """The fitted models that apply: (kind, model) per quantity sources."""
    shell = "open" if open_shell else "closed"
    names = []
    if composite:
        if conventions == "axtb":
            names.append(f"composite axtb {shell}")
        names.append(f"{method} {_kind(method, conventions)} {shell}")
        names.append(f"gfn2 {_kind('gfn2', gfn2_conventions or 'axtb')} {shell}")
    else:
        names.append(f"{method} {_kind(method, conventions)} {shell}")
    # open shells: their own fit (few molecules) and the closed shells' alike, the larger bound wins
    if open_shell:
        names += [n.replace(" open", " closed") for n in names]
    out = []
    for name in names:
        if name in CALIBRATION.get("tables", {}):
            out.append(("table", CALIBRATION["tables"][name]))
        elif name in CALIBRATION.get("scf", {}):
            out.append(("power", CALIBRATION["scf"][name]))
    return out


def _grid() -> list[float]:
    """The measured factors, largest first."""
    return sorted(CALIBRATION.get("accs", [ACC_MAX, ACC_FLOOR]), reverse=True)


def _at_or_above(acc: float) -> float:
    """The smallest measured factor at or above `acc` (errors fall with the factor, so a model's
    bound there bounds `acc` too)."""
    above = [g for g in _grid() if g >= acc * (1 - 1e-12)]
    return above[-1] if above else _grid()[0]


def _acc_model(kind: str, model: Any, target: float) -> float:
    """The largest measured factor whose bound (times the safety factor) meets the target; the
    smallest measured one when none does."""
    for g in _grid():
        e = _estimate(kind, model, g)
        if e is not None and SAFETY * e <= target:
            return g
    return _grid()[-1]


def _estimate(kind: str, model: Any, acc: float) -> float | None:
    """A model's bound at the factor `acc`: at the smallest measured factor at or above it."""
    g = _at_or_above(acc)
    if kind == "power":
        return model["k"] * g ** model["p"]
    for a, err in model:  # a table, largest factor first
        if abs(a - g) <= 1e-12 * g:
            return err
    return model[0][1]


def scf_factor(
    method: str,
    conventions: str,
    composite: bool,
    gfn2_conventions: str | None,
    open_shell: bool,
    target: Accuracy,
    scale: float,
    scale_given: bool,
) -> tuple[float, dict[str, Any], bool]:
    """The SCF accuracy factor of a target (energy and gradient: the Hessian's goes through
    `hessian_factor`), the record's ``thresholds`` and ``estimate`` fields, and whether the floor
    stopped it short of the target.

    Args:
        method, conventions: canonical names; ``composite``: the g-xTB solvation composite with its
            GFN2-xTB parts in ``gfn2_conventions``.
        open_shell: Multiplicity above 1.
        target: The target.
        scale: ``accuracy_scale`` (the factor without a target); ``scale_given``: it was given (a
            floor under the target).
    """
    models = _models(method, conventions, composite, gfn2_conventions, open_shell)
    acc = _grid()[0]
    for q in ("energy", "gradient"):
        t = getattr(target, q)
        if t is None:
            continue
        for kind, model in models:
            if q in model:
                acc = min(acc, _acc_model(kind, model[q], t))
    if scale_given:
        acc = min(acc, scale)
    estimate: dict[str, float] = {}
    limited = False
    for q in ("energy", "gradient"):
        if getattr(target, q) is None:
            continue
        values = [e for kind, model in models if q in model and (e := _estimate(kind, model[q], acc)) is not None]
        if values:
            estimate[q] = max(values)
            limited = limited or estimate[q] > getattr(target, q)
    thresholds: dict[str, Any] = {"accuracy_scale": acc, **_thresholds(method, conventions, open_shell, acc)}
    return acc, {"thresholds": thresholds, "estimate": estimate}, limited


def hessian_factor(
    method: str, hessian_method: str, target: float, default: float, explicit: float | None
) -> tuple[float, float | None, bool]:
    """``hessian.accuracy_scale`` of a Hessian target: from the Hessian's model, never looser than
    the method's default or an explicit value; the estimate and whether the floor stopped it."""
    fit = CALIBRATION.get("hessian", {}).get(f"{method} {hessian_method}")
    if fit is None:
        return min(default, explicit if explicit is not None else default), None, False
    grid = sorted(CALIBRATION.get("hessian_accs", [default]), reverse=True)
    h = next((g for g in grid if SAFETY * fit["k"] * g ** fit["p"] <= target), grid[-1])
    h = min(h, default)
    if explicit is not None:
        h = min(h, explicit)
    above = [g for g in grid if g >= h * (1 - 1e-12)]
    estimate = fit["k"] * (above[-1] if above else grid[0]) ** fit["p"]
    return h, estimate, estimate > target


def met(target: Accuracy) -> dict[str, Any]:
    """The record of a target a method without convergence settings meets (GFN-FF; decision 22)."""
    q = target.model_dump(mode="json", exclude_none=True, exclude={"qcdata_version"})
    return {"target": q, "thresholds": {}, "estimate": {k: 0.0 for k in q}}
