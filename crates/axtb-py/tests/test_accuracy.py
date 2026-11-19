# SPDX-License-Identifier: LGPL-3.0-or-later
"""Accuracy targets (axtb/_accuracy.py, design 7): the factor a target sets, the record, explicit
settings as floors, compute and sessions alike, GFN-FF meeting every target, the floor's
warning, and the model's estimates at or above the errors measured against tight SCFs (three
molecules per method, fresh and warm)."""

import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import Accuracy, ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402
from axtb import _accuracy, _axtb, _keywords  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")
CU = Structure.open(ROOT / "tests/structures/cu-nh3-4.xyz").model_copy(update={"charge": 2, "multiplicity": 2})
REF = {"gfn2": 1e-6, "gxtb-v2.0.1": 1e-4, "gxtb-v1.1": 1e-5}


def inp(method="gxtb-v2.0.1", calctype="gradient", keywords=None, structure=CAFFEINE):
    return ProgramInput(
        program="axtb",
        calctype=calctype,
        structure=structure,
        model={"method": method},
        keywords={"threads": 1, **(keywords or {})},
    )


def test_calibration_is_loaded():
    cal = _accuracy.CALIBRATION
    assert cal["source"] == "bench/results/accuracy.json"
    assert {"gfn2 axtb closed", "gxtb-v2.0.1 axtb closed", "gxtb-v1.1 axtb closed"} <= set(cal["scf"])
    assert (ROOT / cal["source"]).exists()


def test_a_tight_gradient_target_tightens_the_scf():
    out = axtb.compute(inp(keywords={"accuracy": {"gradient": 1e-5}}))
    acc = out.record["accuracy"]
    assert acc["target"] == {"gradient": 1e-5}
    assert acc["thresholds"]["accuracy_scale"] < 1.0
    assert acc["estimate"]["gradient"] <= 1e-5
    assert acc["thresholds"]["energy_tol"] == pytest.approx(1e-6 * acc["thresholds"]["accuracy_scale"])
    assert not out.warnings
    # a loose target loosens the SCF, at most to 10
    loose = axtb.compute(inp(keywords={"accuracy": {"energy": 1e-3}})).record["accuracy"]
    assert 1.0 < loose["thresholds"]["accuracy_scale"] <= 10.0


def test_explicit_accuracy_scale_is_a_floor():
    out = axtb.compute(inp(keywords={"accuracy": {"energy": 1e-3}, "accuracy_scale": 0.5}))
    assert out.record["accuracy"]["thresholds"]["accuracy_scale"] == 0.5
    out = axtb.compute(inp(keywords={"accuracy": {"gradient": 1e-6}, "accuracy_scale": 0.5}))
    assert out.record["accuracy"]["thresholds"]["accuracy_scale"] < 0.5


def test_the_floor_warns():
    out = axtb.compute(inp(keywords={"accuracy": {"gradient": 1e-14}}))
    acc = out.record["accuracy"]
    assert acc["thresholds"]["accuracy_scale"] == _accuracy.ACC_FLOOR
    assert acc["estimate"]["gradient"] > 1e-14
    assert [w.code for w in out.warnings] == ["accuracy_limited"]


def test_compute_and_session_set_the_same_thresholds():
    target = {"gradient": 2e-5, "energy": 1e-7}
    out = axtb.compute(inp(keywords={"accuracy": target}))
    with axtb.Session(inp()) as s:
        f = s.compute(accuracy=Accuracy(**target))
        g = s.compute(lineage="other")  # no target: the session's own factor again
    assert f.record["accuracy"] == out.record["accuracy"]
    assert f.results.energy == out.results.energy
    assert "accuracy" not in g.record
    # a bare number targets the calctype's highest derivative
    with axtb.Session(inp()) as s:
        assert s.compute(accuracy=1e-5).record["accuracy"]["target"] == {"gradient": 1e-5}


def test_hessian_target():
    out = axtb.compute(inp("gfn2", "hessian", {"accuracy": {"hessian": 1e-7}}, structure=WATER))
    acc = out.record["accuracy"]
    assert acc["thresholds"]["hessian.accuracy_scale"] <= 1e-3
    assert "hessian" in acc["estimate"]


def test_gfnff_meets_every_target():
    out = axtb.compute(inp("gfnff", keywords={"accuracy": {"gradient": 1e-9, "energy": 1e-12}}))
    assert out.record["accuracy"] == {
        "target": {"energy": 1e-12, "gradient": 1e-9},
        "thresholds": {},
        "estimate": {"energy": 0.0, "gradient": 0.0},
    }


def test_composite_takes_the_tightest_model():
    kw = {"accuracy": {"gradient": 1e-5}, "solvent": {"name": "water", "correction": "gfn2"}}
    comp = axtb.compute(inp(keywords=kw)).record["accuracy"]["thresholds"]["accuracy_scale"]
    plain = axtb.compute(inp(keywords={"accuracy": {"gradient": 1e-5}})).record["accuracy"]["thresholds"]
    assert comp <= plain["accuracy_scale"]


def errors(r, ref):
    return abs(r["energy"] - ref["energy"]), float(np.abs(r["gradient"] - ref["gradient"]).max())


@pytest.mark.parametrize(
    "method,conventions",
    [
        ("gfn2", "axtb"),
        ("gfn2", "xtb-6.8.0"),
        ("gxtb-v2.0.1", "axtb"),
        ("gxtb-v2.0.1", "gxtb-2.0.1"),
        ("gxtb-v1.1", "axtb"),
        ("gxtb-v1.1", "gxtb-1.1.0"),
    ],
)
@pytest.mark.parametrize("structure", [WATER, CAFFEINE, CU], ids=["water", "caffeine", "cu-doublet"])
@pytest.mark.parametrize("target", [{"gradient": 1e-5}, {"energy": 1e-8, "gradient": 1e-4}])
def test_estimates_bound_the_errors(method, conventions, structure, target):
    """The estimate of a target is at or above the error measured against an SCF at the method's
    floor, fresh and warm (from the structure to a 0.02-bohr displacement)."""
    i = inp(method, keywords={"conventions": conventions, "accuracy": target}, structure=structure)
    out = axtb.compute(i)
    acc = out.record["accuracy"]
    p = _keywords.plan_input(i)
    req = json.dumps({**p.request, "threads": 1})
    task = json.dumps({"numbers": [], "positions": [], "gradient": True})
    ref = _axtb.compute_once(req, p.numbers, p.positions, task, accuracy=REF[method])
    got = {"energy": out.results.energy, "gradient": out.results.gradient}
    de, dg = errors(got, ref)
    est = acc["estimate"]
    if "energy" in est:
        assert de <= est["energy"], (de, est)
    assert dg <= est["gradient"], (dg, est)
    # warm: a session at the structure, then at a displacement, against a tight fresh SCF there
    x = np.asarray(structure.geometry, dtype=float)
    n = len(x)
    x1 = x + 0.01 * np.array([[((7 * a + 3 * c) % 5) - 2 for c in range(3)] for a in range(n)], dtype=float)
    with axtb.Session(i) as s:
        s.compute(x)
        f = s.compute(x1)
    ref1 = _axtb.compute_once(req, p.numbers, x1, task, accuracy=REF[method])
    de, dg = errors({"energy": f.results.energy, "gradient": f.results.gradient}, ref1)
    if "energy" in est:
        assert de <= est["energy"], ("warm", de, est)
    assert dg <= est["gradient"], ("warm", dg, est)
