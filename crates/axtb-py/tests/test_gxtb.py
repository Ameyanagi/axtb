# SPDX-License-Identifier: LGPL-3.0-or-later
"""g-xTB v1.1 through the Python package (axtb 0.1's calls through tests/flat.py).

Run with the extension built from this tree (``maturin develop -m crates/axtb-py/Cargo.toml``)
and the CLI at ``target/release/axtb`` (or ``$AXTB_BIN``)::

    python -m pytest crates/axtb-py/tests/test_gxtb.py

The energies and gradients must equal the CLI's bit for bit (the CLI's equal the Rust API's,
see crates/axtb-cli/tests/gxtb.rs), and with conventions="gxtb" the energies the g-xTB
1.1.0 program's to 1e-11 Eh.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import flat  # tests/flat.py: axtb 0.1's calls through axtb.compute and Session


ROOT = Path(__file__).resolve().parents[3]
REF = ROOT / "crates/gxtb/tests/reference"
AXTB_BIN = Path(os.environ.get("AXTB_BIN", ROOT / "target/release/axtb"))
BOHR_ANGSTROM = float(np.float32(0.52917726))  # the program's factor (axtb_gxtb::BOHR_ANGSTROM)
CASES = [("water", 0, 0), ("dioxygen_triplet", 0, 2), ("tetrachlorocuprate", -2, 1)]


def read(name):
    lines = (REF / name / "mol.xyz").read_text().splitlines()
    nat = int(lines[0])
    numbers, pos = [], []
    for line in lines[2:2 + nat]:
        f = line.split()
        numbers.append({"H": 1, "C": 6, "N": 7, "O": 8, "Cl": 17, "Cu": 29}[f[0]])
        pos.append([float(x) / BOHR_ANGSTROM for x in f[1:4]])
    program = json.loads((REF / name / "summary.json").read_text())["total_energy"]
    return numbers, np.array(pos), program


def cli(numbers, pos, charge, uhf, gradient=False, conventions="gxtb"):
    req = {"numbers": numbers, "positions": pos.tolist(), "charge": charge, "uhf": uhf, "method": "gxtb-v1.1",
           "gradient": gradient, "conventions": conventions}
    out = subprocess.run([str(AXTB_BIN), "--json"], input=json.dumps(req), capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


@pytest.mark.parametrize("name,charge,uhf", CASES)
def test_python_equals_cli_and_program(name, charge, uhf):
    numbers, pos, program = read(name)
    r = flat.singlepoint(numbers, pos, charge=charge, uhf=uhf, method="gxtb-v1.1", conventions="gxtb")
    assert abs(r["energy"] - program) < 1e-11
    calc = flat.Calculator("gxtb-v1.1", numbers, pos, charge=charge, uhf=uhf, threads=2, conventions="gxtb")
    rc = calc.singlepoint(pos, gradient=False)
    assert rc["energy"] == r["energy"]
    assert "gradient" not in r and "gradient" not in rc
    assert np.array_equal(rc["charges"], r["charges"])
    # the program's reference shell occupations are single-precision values, so its Mulliken
    # charges miss the total charge by up to ~1e-7 e per atom (reproduced)
    assert abs(sum(r["charges"]) - charge) < 1e-6
    rg = calc.singlepoint(pos, gradient=True)
    g = flat.singlepoint(numbers, pos, charge=charge, uhf=uhf, method="gxtb-v1.1", gradient=True, conventions="gxtb")
    assert rg["energy"] == g["energy"] == r["energy"]
    assert np.array_equal(rg["gradient"], g["gradient"])
    assert np.asarray(g["gradient"]).shape == (len(numbers), 3)
    assert np.abs(np.asarray(g["gradient"]).sum(axis=0)).max() < 1e-9
    if AXTB_BIN.exists():
        c = cli(numbers, pos, charge, uhf)
        assert c["energy"] == r["energy"]
        assert c["terms"] == r["terms"]
        assert c["charges"] == list(r["charges"])
        c = cli(numbers, pos, charge, uhf, gradient=True)
        assert np.array_equal(np.asarray(c["gradient"]), np.asarray(g["gradient"]))
    # default conventions (exact Fock matrix): slightly lower energy than the program's
    d = flat.singlepoint(numbers, pos, charge=charge, uhf=uhf, method="gxtb-v1.1", gradient=True)
    assert program - 1e-4 < d["energy"] < program
    if AXTB_BIN.exists():
        c = cli(numbers, pos, charge, uhf, gradient=True, conventions="axtb")
        assert c["energy"] == d["energy"]
        assert np.array_equal(np.asarray(c["gradient"]), np.asarray(d["gradient"]))


def test_hessian():
    numbers, pos, _ = read("water")
    calc = flat.Calculator("gxtb-v1.1", numbers, pos, threads=2)
    h = calc.hessian(pos)
    hess = np.asarray(h["hessian"])
    assert hess.shape == (9, 9)
    assert np.array_equal(hess, hess.T)
    assert h["hessian_max_asymmetry"] < 2e-5
    # translational invariance of every row
    assert np.abs(hess.reshape(9, 3, 3).sum(axis=1)).max() < 1e-5
    serial = flat.Calculator("gxtb-v1.1", numbers, pos).hessian(pos)
    assert np.allclose(serial["hessian"], hess, atol=1e-12)


def test_errors():
    numbers, pos, _ = read("water")
    with pytest.raises(ValueError, match="no implicit-solvent parameters.*solvent.correction"):
        flat.singlepoint(numbers, pos, method="gxtb-v1.1", alpb="water")
    # Fermi smearing (the 1.1.0 program's -tel): program energy at 1000 K
    rows = json.loads((ROOT / "crates/gxtb/tests/reference/etemp.json").read_text())["rows"]
    ref = next(r for r in rows if r["name"] == "water" and r["etemp"] == 1000.0)
    hot = flat.singlepoint(numbers, pos, method="gxtb-v1.1", etemp=1000.0, conventions="gxtb")
    assert abs(hot["energy"] - ref["total_energy"]) < 1e-11
    assert hot["iterations"] == ref["scf_iterations"]
    with pytest.raises(ValueError, match="conventions: unknown value"):
        flat.singlepoint(numbers, pos, method="gxtb-v1.1", conventions="xtb")
    # GFN2 keeps its defaults (etemp 300 K, max_iter 250, conventions axtb)
    default = flat.singlepoint(numbers, pos, method="gfn2", gradient=True)
    explicit = flat.singlepoint(numbers, pos, method="gfn2", gradient=True, etemp=300.0, max_iter=250,
                                conventions="axtb", accuracy=1.0)
    assert default["energy"] == explicit["energy"]
    assert np.array_equal(default["gradient"], explicit["gradient"])
    # a solvent's reference state without the solvent is refused, not ignored
    with pytest.raises(flat.Failed, match="solvent.state"):
        flat.singlepoint(numbers, pos, method="gfn2", solvation_state="gsolv")
