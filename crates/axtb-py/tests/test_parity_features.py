# SPDX-License-Identifier: LGPL-3.0-or-later
"""SCF guess reuse, Mayer-Wiberg bond orders, orbital levels, HOMO-LUMO gap and the molecular
quadrupole through the Python package (axtb 0.1's calls through tests/flat.py, and sessions'
lineages) and the CLI.

    python -m pytest crates/axtb-py/tests/test_parity_features.py

Needs the extension built from this tree and the CLI at ``target/release/axtb`` (or
``$AXTB_BIN``). The comparison with the programs is ``tests/gxtb_parity_features.py``.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import flat  # tests/flat.py: axtb 0.1's calls through axtb.compute and Session

import axtb

ROOT = Path(__file__).resolve().parents[3]
AXTB_BIN = Path(os.environ.get("AXTB_BIN", ROOT / "target/release/axtb"))
A = 1.0 / 0.529177210903
# water (tests/structures/water.xyz), bohr
NUMBERS = [1, 1, 8]
POS = A * np.array([[-0.0121648342, 0.0, -0.9323517940], [0.8750977221, 0.0, 0.3219980096],
                    [-0.0543635745, 0.0, 0.0384514415]])
METHODS = ["gfn2", "gxtb-v1.1", "gxtb-v2.0.1"]
KEYS = {"guess", "orbital_energies_alpha", "orbital_energies_beta", "occupations_alpha", "occupations_beta",
        "homo_lumo_gap", "quadrupole"}


@pytest.mark.parametrize("method", METHODS)
def test_singlepoint_properties(method):
    r = flat.singlepoint(NUMBERS, POS, method=method)
    assert KEYS <= set(r) and "bond_orders" not in r
    assert r["guess"] == ("eeq" if method == "gfn2" else "eeq-bc")
    # closed shell: the same orbitals for both spins, half the occupations each
    assert np.array_equal(r["orbital_energies_alpha"], r["orbital_energies_beta"])
    assert abs(r["occupations_alpha"].sum() - 4.0) < 1e-6 and abs(r["occupations_beta"].sum() - 4.0) < 1e-6
    e = r["orbital_energies_alpha"]
    assert abs(r["homo_lumo_gap"] - (e[4] - e[3])) < 1e-14
    q = r["quadrupole"]
    assert q.shape == (6,) and abs(q[0] + q[2] + q[5]) < 1e-12
    r = flat.singlepoint(NUMBERS, POS, method=method, bond_orders=True)
    bo = r["bond_orders"]
    assert bo.shape == (3, 3) and np.allclose(bo, bo.T) and np.all(np.diag(bo) == 0.0)
    assert 0.85 < bo[0, 2] < 0.95 and bo[0, 1] < 0.01


def test_open_shell_levels_and_bond_orders():
    o2 = [8, 8], A * np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.21]])
    for method in METHODS:
        r = flat.singlepoint(*o2, uhf=2, method=method, bond_orders=True)
        assert abs(r["occupations_alpha"].sum() - 7.0) < 1e-6 and abs(r["occupations_beta"].sum() - 5.0) < 1e-6
        if method == "gfn2":  # one set of orbitals
            assert np.array_equal(r["orbital_energies_alpha"], r["orbital_energies_beta"])
        else:
            assert not np.array_equal(r["orbital_energies_alpha"], r["orbital_energies_beta"])
        # Mayer's open-shell bond order 2 sum_s (P^s S)(P^s S): close to 2 for triplet O2
        assert abs(r["bond_orders"][0, 1] - 2.0) < 0.05


@pytest.mark.parametrize("method", ["gxtb-v1.1", "gxtb-v2.0.1"])
def test_guess_reuse(method):
    moved = POS + 0.01
    calc = flat.Calculator(method, NUMBERS, POS)
    a = calc.singlepoint(POS)
    b = calc.singlepoint(moved, gradient=True)
    assert (a["guess"], b["guess"]) == ("eeq-bc", "previous")
    fresh = flat.singlepoint(NUMBERS, moved, method=method, gradient=True)
    assert b["iterations"] < fresh["iterations"]
    assert abs(b["energy"] - fresh["energy"]) < 1e-6
    assert np.max(np.abs(b["gradient"] - fresh["gradient"])) < 1e-4
    # program conventions: the program's start every time, unless asked
    calc = flat.Calculator(method, NUMBERS, POS, conventions="gxtb")
    assert [calc.singlepoint(p)["guess"] for p in (POS, moved)] == ["eeq-bc", "eeq-bc"]
    calc = flat.Calculator(method, NUMBERS, POS, conventions="gxtb", guess="previous")
    assert [calc.singlepoint(p)["guess"] for p in (POS, POS)] == ["eeq-bc", "previous"]
    calc = flat.Calculator(method, NUMBERS, POS, guess="eeq-bc")
    assert [calc.singlepoint(p)["guess"] for p in (POS, moved)] == ["eeq-bc", "eeq-bc"]
    # a converged wavefunction seeds another lineage (0.1: guess_state / set_guess_state)
    inp = flat.program_input({"numbers": NUMBERS, "positions": POS.tolist(), "method": method})
    with axtb.Session(inp) as s:
        s.compute(POS, lineage="src")
        f = s.compute(moved, lineage="dst", seed=["src"])
        assert (f.record["scf"]["start"], f.record["from_lineage"]) == ("previous", "src")
        assert s.compute(moved, lineage="other").record["scf"]["start"] == "eeq-bc"


def test_gfn2_guess_option():
    calc = flat.Calculator("gfn2", NUMBERS, POS)
    assert [calc.singlepoint(p)["guess"] for p in (POS, POS + 0.01)] == ["eeq", "previous"]
    calc = flat.Calculator("gfn2", NUMBERS, POS, guess="eeq")
    assert [calc.singlepoint(p)["guess"] for p in (POS, POS + 0.01)] == ["eeq", "eeq"]
    with pytest.raises(ValueError, match="guess"):
        flat.Calculator("gfn2", NUMBERS, POS, guess="huckel")
    with pytest.raises(ValueError, match="scf.restart: GFN-FF has no SCF"):
        flat.Calculator("gfnff", NUMBERS, POS, guess="previous")


def test_composite_reuses_the_gxtb_part():
    kw = {"alpb": "water", "solvation_correction": "gfn2"}
    calc = flat.Calculator("gxtb-v2.0.1", NUMBERS, POS, **kw)
    a, b = calc.singlepoint(POS), calc.singlepoint(POS + 0.01)
    assert (a["guess"], b["guess"]) == ("eeq-bc", "previous")
    assert b["composite"]["gxtb_reuse_guess"] and b["composite"]["reuse_guess"]
    calc = flat.Calculator("gxtb-v2.0.1", NUMBERS, POS, conventions="gxtb", **kw)
    a, b = calc.singlepoint(POS), calc.singlepoint(POS + 0.01)
    assert (a["guess"], b["guess"]) == ("eeq-bc", "eeq-bc")
    assert not b["composite"]["gxtb_reuse_guess"]
    # the state of all three parts seeds another lineage
    inp = flat.program_input({"numbers": NUMBERS, "positions": POS.tolist(), "method": "gxtb-v2.0.1", **kw})
    with axtb.Session(inp) as s:
        s.compute(POS, lineage="src")
        f = s.compute(POS, lineage="dst", seed=["src"])
        assert f.record["scf"]["start"] == "previous" and f.record["from_lineage"] == "src"
