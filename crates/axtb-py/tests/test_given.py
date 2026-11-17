# SPDX-License-Identifier: LGPL-3.0-or-later
"""Given results as the SCF's start (decision 19, design 4.6): a converged frame given back
reproduces its result within the accuracy model's estimate in at most three iterations (g-xTB
from the orbitals' density, GFN2-xTB from the shell populations and atomic multipoles), sessions
start every request without history from it, and misfits are refused."""

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramError, ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402
from axtb import _accuracy  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")
CU = Structure.open(ROOT / "tests/structures/cu-nh3-4.xyz").model_copy(update={"charge": 2, "multiplicity": 2})


def inp(method, structure=CAFFEINE, conventions="axtb", given=None, output=True, calctype="gradient", **kw):
    keywords = {"threads": 1, "conventions": conventions, **kw}
    if output and method != "gfnff":
        keywords["output"] = {"orbitals": True}
    if method == "composite":
        method = "gxtb-v2.0.1"
        keywords["solvent"] = {"name": "water", "correction": "gfn2"}
    return ProgramInput(
        program="axtb", calctype=calctype, structure=structure, model={"method": method}, keywords=keywords, given=given
    )


def bound(method, conventions, structure, q):
    """The accuracy model's estimate at the default factor (1)."""
    if method == "composite":
        method = "gxtb-v2.0.1"
    kind = "axtb" if conventions == "axtb" else "program"
    shell = "open" if structure.multiplicity > 1 else "closed"
    fits = [_accuracy.CALIBRATION["scf"].get(f"{method} {kind} {s}") for s in (shell, "closed")]
    tables = [_accuracy.CALIBRATION["tables"].get(f"{method} {kind} {shell}")]
    values = [f[q]["k"] for f in fits if f] + [dict(map(tuple, t[q]))[1.0] for t in tables if t]
    return max(values)


@pytest.mark.parametrize(
    "method,conventions,structure",
    [
        ("gfn2", "axtb", CAFFEINE),
        ("gfn2", "xtb-6.8.0", CAFFEINE),
        ("gfn2", "axtb", CU),
        ("gxtb-v2.0.1", "axtb", CAFFEINE),
        ("gxtb-v2.0.1", "gxtb-2.0.1", CAFFEINE),
        ("gxtb-v2.0.1", "axtb", CU),
        ("gxtb-v1.1", "axtb", CAFFEINE),
        ("gxtb-v1.1", "gxtb-1.1.0", CAFFEINE),
        ("gxtb-v1.1", "axtb", CU),
        ("composite", "axtb", CAFFEINE),
    ],
    ids=lambda x: x if isinstance(x, str) else ("doublet" if x.multiplicity > 1 else "caffeine"),
)
def test_a_converged_frame_given_back(method, conventions, structure):
    first = axtb.compute(inp(method, structure, conventions))
    again = axtb.compute(inp(method, structure, conventions, given=first.results))
    assert again.record["scf"]["start"] == "given"
    if method != "composite":
        assert again.resolved.settings["scf"]["guess"] == "given"
    assert again.results.scf_iterations <= 3, again.results.scf_iterations
    de = abs(again.results.energy - first.results.energy)
    dg = np.abs(again.results.gradient - first.results.gradient).max()
    assert de <= bound(method, conventions, structure, "energy"), de
    assert dg <= bound(method, conventions, structure, "gradient"), dg


def test_sessions_start_requests_without_history_from_it():
    first = axtb.compute(inp("gxtb-v2.0.1"))
    with axtb.Session(inp("gxtb-v2.0.1", given=first.results)) as s:
        f = s.compute(lineage="a")
        assert f.record["scf"]["start"] == "given" and f.record["warm_start"] is False
        g = s.compute(np.asarray(CAFFEINE.geometry) + 0.01, lineage="a")
        assert g.record["scf"]["start"] == "previous" and g.record["warm_start"] is True
        h = s.compute(lineage="b")  # another lineage without history: the given start again
        assert h.record["scf"]["start"] == "given" and h.results.energy == f.results.energy


@pytest.mark.parametrize(
    "method,given_by,code",
    [
        ("gfnff", "gxtb-v2.0.1", "not_applicable"),
        ("gxtb-v2.0.1", "gfn2", "not_applicable"),  # GFN2-xTB's basis set
        ("gfn2", "gxtb-v2.0.1", "not_applicable"),  # g-xTB's shells (Cu: s, p, d; GFN2-xTB d, s, p)
    ],
)
def test_misfits_are_refused(method, given_by, code):
    given = axtb.compute(inp(given_by, CU)).results
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp(method, CU, given=given))
    assert (e.value.error.kind.value, e.value.error.code) == ("input", code)


def test_another_methods_moments_on_the_same_shells_are_a_guess():
    """A starting guess needs no match of method (qcdata's conventions): g-xTB's populations and
    multipoles of caffeine (H s; C, N, O s and p, as GFN2-xTB's shells) start GFN2-xTB."""
    given = axtb.compute(inp("gxtb-v2.0.1")).results
    out = axtb.compute(inp("gfn2", given=given))
    assert out.record["scf"]["start"] == "given"
    plain = axtb.compute(inp("gfn2"))
    assert abs(out.results.energy - plain.results.energy) < 1e-6


def test_nothing_usable_is_refused():
    no_orbitals = axtb.compute(inp("gxtb-v2.0.1", output=False)).results
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gxtb-v2.0.1", given=no_orbitals))
    assert e.value.error.code == "not_applicable" and "orbitals" in e.value.error.message
    # other spin channels
    doublet = axtb.compute(inp("gxtb-v2.0.1", CU)).results
    cation = CU.model_copy(update={"charge": 1, "multiplicity": 1})
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gxtb-v2.0.1", cation, given=doublet))
    assert e.value.error.code == "not_applicable"
    # a given from another charge: the density holds other electrons (the driver's check)
    neutral = Structure.open(ROOT / "tests/structures/caffeine.xyz")
    dication = neutral.model_copy(update={"charge": 2})
    given = axtb.compute(inp("gxtb-v2.0.1")).results
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gxtb-v2.0.1", dication, given=given))
    assert e.value.error.kind.value == "input" and "electrons" in e.value.error.message


def test_scf_guess_given():
    given = axtb.compute(inp("gfn2")).results
    out = axtb.compute(inp("gfn2", given=given, scf={"guess": "given"}))
    assert out.record["scf"]["start"] == "given"
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gfn2", given=given, scf={"guess": "sad"}))
    assert e.value.error.code == "no_effect"
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp("gfn2", scf={"guess": "given"}))
    assert e.value.error.code == "missing_setting"
