# SPDX-License-Identifier: LGPL-3.0-or-later
"""``scf.retries``, the SCF retry ladder: a rung's result is the calculation at the rung's
settings bit for bit, a ladder that does not fire changes no number, a Hessian runs whole at
one rung, a session frame from a rung commits and the next frame tries the calculation's own
settings first, the composite climbs the given rungs per part, every rung failing is
``scf_not_converged`` with every attempt recorded, and the refusals (the floor of
``scf.etemp``, the rungs' shapes, GFN-FF).

Every SCF that must fail fails on its iteration limit, far below what it needs, and every one
that must converge has iterations to spare: water and caffeine at their equilibrium geometries
(one SCF solution) take 8 to 11 iterations from the method's own start in every method, and 2 or
3 from a converged state of the same geometry, so no outcome depends on round-off (the CPU)."""

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramError, ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
AXTB_BIN = Path(os.environ.get("AXTB_BIN", ROOT / "target/release/axtb"))
WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")
# 5 iterations: too few for an SCF from the method's own start, enough for one from a converged
# state of the same geometry
SHORT = 5
LADDER = [{"etemp": 1000, "max_iter": 250}, {"etemp": 3000, "max_iter": 250}]
COMPOSITE = {"solvent": {"name": "water", "correction": "gfn2"}}
METHODS = ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"]


def inp(method, calctype="gradient", keywords=None, structure=WATER):
    return ProgramInput(program="axtb", calctype=calctype, structure=structure, model={"method": method},
                        keywords={"threads": 1, **(keywords or {})})


def short(retries=LADDER, **extra):
    """Keywords with the short iteration limit and `retries` (None: without), and `extra`."""
    kw = {"scf": {"max_iter": SHORT}}
    if retries is not None:
        kw["scf"]["retries"] = retries
    for k, v in extra.items():
        if isinstance(v, dict):
            kw.setdefault(k, {}).update(v)
        else:
            kw[k] = v
    return kw


def at_rung(**extra):
    """The keywords of a calculation at the first rung's settings, without a ladder."""
    return short(None, scf={"etemp": 1000, "max_iter": 250}, **extra)


def bits(r):
    out = [r.energy, r.atomic_charges.tobytes(), r.scf_iterations]
    for x in (r.gradient, r.hessian):
        if x is not None:
            out.append(np.asarray(x).tobytes())
    return out


@pytest.mark.parametrize("method", METHODS)
def test_a_rung_is_the_calculation_at_its_settings(method):
    out = axtb.compute(inp(method, keywords=short()))
    rec = out.record
    assert rec["contingent"]["scf.retries"] == {"fired": True, "rung": 1}
    assert out.resolved.settings["scf"]["retries"] == [{"etemp": 1000.0, "max_iter": 250},
                                                       {"etemp": 3000.0, "max_iter": 250}]
    attempts = rec["scf"]["attempts"]
    base_etemp = 300.0 if method == "gfn2" else 0.0
    first = attempts[0]
    assert (first["retry"], first["etemp"], first["max_iter"], first["converged"], first["error"]) == (
        0, base_etemp, SHORT, False, "scf_not_converged")
    last = attempts[-1]
    assert (last["retry"], last["etemp"], last["max_iter"], last["converged"]) == (1, 1000.0, 250, True)
    assert last["start"] == rec["scf"]["start"] and last["iterations"] == out.results.scf_iterations
    assert [w.code for w in out.warnings] == ["scf_retry"]
    assert out.warnings[0].message == "converged on retry 1 of 2: etemp 1000 K, max_iter 250"
    assert out.results.wavefunction.orbitals.electronic_temperature == 1000.0
    # bit for bit the calculation at the rung's settings
    direct = axtb.compute(inp(method, keywords=at_rung()))
    assert bits(out.results) == bits(direct.results)


def test_a_failed_rung_climbs_to_the_next():
    # the first rung keeps the short iteration limit: it fails too
    out = axtb.compute(inp("gfn2", keywords=short([{"etemp": 500}, {"etemp": 1000, "max_iter": 250}])))
    assert out.record["contingent"]["scf.retries"] == {"fired": True, "rung": 2}
    a = out.record["scf"]["attempts"]
    assert [(x["retry"], x["etemp"], x["max_iter"]) for x in a if not x["converged"]] == [
        (0, 300.0, SHORT), (1, 500.0, SHORT)]
    assert (a[-1]["retry"], a[-1]["converged"]) == (2, True)
    assert out.warnings[0].message == "converged on retry 2 of 2: etemp 1000 K, max_iter 250"


@pytest.mark.parametrize("method", ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1", "composite"])
@pytest.mark.parametrize("calctype", ["gradient", "hessian"])
def test_a_ladder_that_does_not_fire_changes_no_number(method, calctype):
    if method == "gxtb-v1.1" and calctype == "hessian":
        pytest.skip("v1.1's numerical Hessian: covered by the gradient and the other methods")
    m = "gxtb-v2.0.1" if method == "composite" else method
    extra = COMPOSITE if method == "composite" else {}
    plain = axtb.compute(inp(m, calctype, extra))
    laddered = axtb.compute(inp(m, calctype, {**extra, "scf": {"retries": LADDER}}))
    assert bits(plain.results) == bits(laddered.results)
    assert laddered.record["contingent"]["scf.retries"]["fired"] is False
    assert not [w for w in laddered.warnings if w.code == "scf_retry"]
    assert "retries" not in plain.resolved.settings["scf"]


@pytest.mark.parametrize("method", ["gfn2", "gxtb-v2.0.1"])
def test_a_hessian_runs_whole_at_one_rung(method):
    out = axtb.compute(inp(method, "hessian", short()))
    assert out.record["contingent"]["scf.retries"] == {"fired": True, "rung": 1}
    # the Hessian and its single point at the rung, bit for bit the calculation at its settings
    direct = axtb.compute(inp(method, "hessian", at_rung()))
    assert bits(out.results) == bits(direct.results)
    assert np.array_equal(out.results.vibrational_analysis.frequencies,
                          direct.results.vibrational_analysis.frequencies)


@pytest.mark.parametrize("method", METHODS)
def test_session_frames_commit_and_try_the_own_settings_first(method):
    x = np.asarray(WATER.geometry, dtype=float)
    with axtb.Session(inp(method, keywords=short(scf={"restart": "previous"}))) as s:
        f1 = s.compute(x)
        f2 = s.compute(x)
    # the first frame came from the rung, as compute's does
    want = axtb.compute(inp(method, keywords=short()))
    assert bits(f1.results) == bits(want.results)
    assert f1.record["contingent"]["scf.retries"]["fired"] is True
    assert [w.code for w in f1.warnings] == ["scf_retry"]
    # the second starts from the first's committed state at the calculation's own settings,
    # where it converges
    a = f2.record["scf"]["attempts"]
    assert len(a) == 1 or method == "gfn2"  # GFN2-xTB: one entry per mixer of the try
    assert (a[-1]["retry"], a[-1]["etemp"], a[-1]["max_iter"], a[-1]["converged"]) == (
        0, 300.0 if method == "gfn2" else 0.0, SHORT, True)
    assert f2.record["scf"]["start"] == "previous" and f2.record["warm_start"] is True
    assert f2.record["contingent"]["scf.retries"] == {"fired": False, "rung": 0}
    assert not f2.warnings


def test_composite_climbs_the_given_rungs_per_part():
    out = axtb.compute(inp("gxtb-v2.0.1", "hessian", short(**COMPOSITE)))
    rec = out.record
    assert "composite.retries" not in rec["contingent"]
    c = rec["contingent"]["scf.retries"]
    assert c["fired"] is True and sorted(c["parts"]) == ["gfn2_gas", "gfn2_solvated", "gxtb"]
    definition = rec["composite"]["definition"]
    assert definition["retries"] is False
    assert definition["scf_retries"] == "etemp 1000 K, max_iter 250; etemp 3000 K, max_iter 250"
    parts, hparts = rec["composite"]["parts"], rec["composite"]["hessian_parts"]
    for key in ("gxtb", "gfn2_gas", "gfn2_solvated"):
        # the single point after the Hessian starts at the Hessian's rung: one rung per part
        assert parts[key]["rung"] == hparts[key]["rung"] == 1
        assert parts[key]["retry"] == "retry 1 of 2: etemp 1000 K, max_iter 250"
    attempts = rec["scf"]["attempts"]
    assert {a["part"] for a in attempts} == {"gxtb", "gfn2_gas", "gfn2_solvated"}
    assert any(a.get("hessian") for a in attempts)
    assert sorted(w.message for w in out.warnings if w.code == "scf_retry") == sorted(
        f"the {k} part converged on retry 1 of 2: etemp 1000 K, max_iter 250" for k in parts
    )
    # the composite at the rung's settings (every part at 1000 K): its Hessian bit for bit; the
    # single point after it starts the parts fresh at the rung, not from their Hessians'
    # references, so it agrees within the SCF thresholds
    direct = axtb.compute(inp("gxtb-v2.0.1", "hessian", at_rung(**COMPOSITE)))
    assert np.array_equal(out.results.hessian, direct.results.hessian)
    assert abs(out.results.energy - direct.results.energy) < 1e-10
    assert np.abs(out.results.gradient - direct.results.gradient).max() < 1e-6
    # a gradient alone: bit for bit
    gradient = axtb.compute(inp("gxtb-v2.0.1", "gradient", short(**COMPOSITE)))
    assert bits(gradient.results) == bits(axtb.compute(inp("gxtb-v2.0.1", "gradient", at_rung(**COMPOSITE))).results)
    # an explicit temperature: the ladder starts there
    hot = axtb.compute(inp("gxtb-v2.0.1", "gradient", short([{"max_iter": 250}], scf={"etemp": 1000}, **COMPOSITE)))
    assert {p["etemp"] for p in hot.record["composite"]["parts"].values()} == {1000.0}
    assert {p["rung"] for p in hot.record["composite"]["parts"].values()} == {1}


def test_every_rung_failing_is_scf_not_converged():
    # caffeine needs 10 or 11 iterations in every method: 3, 4 and 5 all fail
    kw = {"scf": {"max_iter": 3, "retries": [{"max_iter": 4}, {"etemp": 1000, "max_iter": 5}]}}
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gfn2", keywords=kw, structure=CAFFEINE))
    err = e.value.error
    assert (err.kind.value, err.code) == ("convergence", "scf_not_converged")
    out = e.value.output
    attempts = out.record["scf"]["attempts"]
    assert [a["retry"] for a in attempts] == [0, 1, 2] and not any(a["converged"] for a in attempts)
    assert [(a["etemp"], a["max_iter"]) for a in attempts] == [(300.0, 3), (300.0, 4), (1000.0, 5)]
    assert {a["error"] for a in attempts} == {"scf_not_converged"}
    # the details are the last rung's
    assert err.details["iterations"] == attempts[-1]["iterations"]
    assert out.record["error"] == err.details
    # a session's failed frame too
    with axtb.Session(inp("gfn2", keywords=kw, structure=CAFFEINE)) as s:
        with pytest.raises(ProgramError) as e:
            s.compute()
    assert e.value.frame.record["scf"]["attempts"] == attempts


@pytest.mark.parametrize("keywords,code,path", [
    ({"scf": {"etemp": 1000, "retries": [{"etemp": 3000}, {"etemp": 500}]}}, "no_effect", "scf.retries[1].etemp"),
    ({"scf": {"retries": []}}, "invalid_value", "scf.retries"),
    ({"scf": {"retries": [{}]}}, "invalid_value", "scf.retries[0]"),
    ({"scf": {"retries": [1000]}}, "invalid_value", "scf.retries[0]"),
    ({"scf": {"retries": {"etemp": 1000}}}, "invalid_value", "scf.retries"),
    ({"scf": {"retries": [{"etemp": -1}]}}, "invalid_value", "scf.retries[0].etemp"),
    ({"scf": {"retries": [{"max_iter": 0}]}}, "invalid_value", "scf.retries[0].max_iter"),
    ({"scf": {"retries": [{"max_iter": 10.5}]}}, "invalid_value", "scf.retries[0].max_iter"),
    ({"scf": {"retries": [{"temp": 1000}]}}, "unknown_key", "scf.retries[0].temp"),
    ({"retries": [{"etemp": 1000}]}, "unknown_key", "retries"),
])
def test_refusals(keywords, code, path):
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gfn2", keywords=keywords))
    assert (e.value.error.kind.value, e.value.error.code) == ("input", code)
    assert e.value.error.message.startswith(path)


def test_the_floor_admits_rungs_at_and_above_it():
    out = axtb.compute(inp("gfn2", keywords={"scf": {"etemp": 1000, "retries": [{"etemp": 1000}, {"max_iter": 500}]}}))
    assert out.resolved.settings["scf"]["retries"] == [{"etemp": 1000.0}, {"max_iter": 500}]


def test_gfnff_has_no_scf():
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gfnff", keywords={"scf": {"retries": LADDER}}))
    assert e.value.error.code == "no_effect" and e.value.error.message.startswith("scf.retries")


def test_schema_and_plan():
    props = axtb.program.settings_schema("gradient")["properties"]["scf"]["properties"]
    assert props["retries"]["items"]["properties"]["etemp"]["x-unit"] == "K"
    plan = axtb.program.plan(inp("gfn2", keywords={"scf": {"retries": LADDER}}))
    assert "scf.retries" in plan["record"]["contingent"]


@pytest.mark.skipif(not AXTB_BIN.exists(), reason="needs the CLI (target/release/axtb)")
def test_the_cli_flag():
    xyz = ROOT / "tests/structures/water.xyz"
    p = subprocess.run([str(AXTB_BIN), str(xyz), "--grad", "--max-iter", str(SHORT),
                        "--scf-retry", "etemp=1000,max_iter=250", "--scf-retry", "etemp=3000,max_iter=250"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    r = json.loads(p.stdout)
    assert r["scf_retry"] == {"rung": 1, "of": 2, "label": "retry 1 of 2: etemp 1000 K, max_iter 250"}
    assert "warning: converged on retry 1 of 2: etemp 1000 K, max_iter 250" in p.stderr
    bad = subprocess.run([str(AXTB_BIN), str(xyz), "--scf-retry", "temp=3"], capture_output=True, text=True)
    assert bad.returncode == 2 and "unknown key" in bad.stderr
