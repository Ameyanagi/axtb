# SPDX-License-Identifier: LGPL-3.0-or-later
"""Vibrational analysis through the Python package (axtb 0.1's calls through tests/flat.py) and the CLI.

Run with the extension built from this tree (a wheel built from ``crates/axtb-py``) and the CLI
at ``target/release/axtb`` (or ``$AXTB_BIN``)::

    python -m pytest crates/axtb-py/tests/test_vibrations.py

The comparison with the xtb program (frequencies, thermochemistry, symmetry numbers) is
``tests/vibrational_check.py``; here: the interfaces agree with each other and with the
textbook formulas.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

import flat  # tests/flat.py: axtb 0.1's calls through axtb.compute and Session


ROOT = Path(__file__).resolve().parents[3]
AXTB_BIN = Path(os.environ.get("AXTB_BIN", ROOT / "target/release/axtb"))
KB = 3.1668115634564e-6  # CODATA 2022 kelvin-hartree relationship

WATER = ([8, 1, 1], np.array([[0.0, 0.0, 0.2214], [0.0, 1.4304, -0.8857], [0.0, -1.4304, -0.8857]]))
CO2 = ([6, 8, 8], np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.1953], [0.0, 0.0, -2.1953]]))


def unit(modes):
    return np.allclose(np.linalg.norm(modes.reshape(len(modes), -1), axis=1), 1.0)


def test_calculator_hessian_has_vibrations():
    numbers, pos = WATER
    h = flat.Calculator("gfn2", numbers, pos).hessian(pos)
    assert h["frequencies"].shape == (3,)
    assert h["normal_modes"].shape == (3, 3, 3) and unit(h["normal_modes"])
    assert h["dipole_derivatives"].shape == (9, 3)
    assert h["ir_intensities"].shape == (3,) and (h["ir_intensities"] > 0).all()
    assert h["vibrational_conventions"] == "axtb" and not h["linear"]
    t = h["thermo"]
    assert (t["point_group"], t["symmetry_number"]) == ("C2v", 2)
    assert t["free_energy"] == pytest.approx(h["energy"] + t["g_rrho"], abs=1e-12)
    assert t["g_rrho"] == pytest.approx(t["enthalpy"] - t["ts"], abs=1e-15)
    assert t["ts"] == pytest.approx(t["temperature"] * t["entropy"], rel=1e-14)
    # the standalone analysis of the same Hessian is the same analysis
    v = flat.analysis(numbers, pos, h["hessian"], dipole_derivatives=h["dipole_derivatives"],
                                  energy=h["energy"])
    assert np.array_equal(v["frequencies"], h["frequencies"])
    assert np.array_equal(v["ir_intensities"], h["ir_intensities"])
    assert v["thermo"]["free_energy"] == t["free_energy"]
    # xtb conventions of the calculation select the xtb analysis
    hx = flat.Calculator("gfn2", numbers, pos, conventions="xtb").hessian(pos)
    assert hx["vibrational_conventions"] == "xtb-6.8.0"
    assert "frequencies" not in flat.Calculator("gfn2", numbers, pos).hessian(pos, vibrations=False)


def test_cli_matches_python():
    if not AXTB_BIN.exists():
        pytest.skip("no CLI")
    numbers, pos = WATER
    req = {"numbers": numbers, "positions": pos.tolist(), "hessian": True, "temperature": [250.0, 300.0]}
    d = json.loads(subprocess.run([str(AXTB_BIN), "--json"], input=json.dumps(req), capture_output=True, text=True,
                                  check=True).stdout)
    h = flat.Calculator("gfn2", numbers, pos).hessian(pos, temperature=[250.0, 300.0])
    assert np.allclose(d["frequencies"], h["frequencies"], rtol=0, atol=1e-9)
    assert d["thermo"]["temperature"] == 300.0 and len(d["thermo"]["points"]) == 2
    assert d["thermo"]["g_rrho"] == pytest.approx(h["thermo"]["g_rrho"], abs=1e-12)
    # analysis only, like `xtb thermo`
    a = json.loads(subprocess.run([str(AXTB_BIN), "--json"], capture_output=True, text=True, check=True,
                                  input=json.dumps({"numbers": numbers, "positions": pos.tolist(),
                                                    "analyze_hessian": d["hessian"], "analysis_conventions": "xtb"})).stdout)
    assert a["vibrational_conventions"] == "xtb" and len(a["frequencies"]) == 3


def test_thermo_options():
    numbers, pos = WATER
    h = flat.Calculator("gfn2", numbers, pos).hessian(pos, vibrations=False)["hessian"]
    base = flat.analysis(numbers, pos, h)["thermo"]
    # symmetry number: G changes by kT ln(sigma) (rotational entropy -k ln sigma)
    sym1 = flat.analysis(numbers, pos, h, symmetry_number=1)["thermo"]
    assert sym1["g_rrho"] - base["g_rrho"] == pytest.approx(-KB * 298.15 * math.log(2), rel=1e-10)
    # pressure: translational entropy -k ln p
    p2 = flat.analysis(numbers, pos, h, pressure=2 * 101325.0)["thermo"]
    assert p2["g_rrho"] - base["g_rrho"] == pytest.approx(KB * 298.15 * math.log(2), rel=1e-10)
    # several temperatures: the last is reported, all are kept
    many = flat.analysis(numbers, pos, h, temperature=[298.15, 500.0])["thermo"]
    assert many["temperature"] == 500.0 and many["points"][0]["g_rrho"] == pytest.approx(base["g_rrho"], abs=1e-15)
    # rotor cutoff 0: harmonic oscillators only; frequencies here are far above 50 cm-1 anyway
    ho = flat.analysis(numbers, pos, h, rotor_cutoff=0.0)["thermo"]
    assert ho["g_rrho"] == pytest.approx(base["g_rrho"], abs=1e-9)
    # isotopes: D2O frequencies drop by about 1/sqrt(2)
    d2o = flat.analysis(numbers, pos, h, masses=[15.999, 2.014, 2.014])
    assert d2o["frequencies"][2] / flat.analysis(numbers, pos, h)["frequencies"][2] == pytest.approx(0.73, abs=0.02)
    assert "thermo" not in flat.analysis(numbers, pos, h, thermo=False)
    with pytest.raises(ValueError):
        flat.analysis(numbers, pos, h, temperature=0.0)
    with pytest.raises(ValueError):
        flat.analysis(numbers, pos, h[:6, :6])
    with pytest.raises(ValueError):
        flat.analysis(numbers, pos, h, conventions="tblite")


def test_linear_and_gxtb():
    numbers, pos = CO2
    h = flat.Calculator("gfn2", numbers, pos).hessian(pos)
    assert h["linear"] and h["frequencies"].shape == (4,)
    assert (h["thermo"]["point_group"], h["thermo"]["symmetry_number"]) == ("Dinfh", 2)
    numbers, pos = WATER
    g = flat.Calculator("gxtb-v2.0.1", numbers, pos, conventions="gxtb").hessian(pos)
    assert g["vibrational_conventions"] == "xtb-6.7.1" and g["frequencies"].shape == (3,)
    assert g["ir_intensities"].shape == (3,) and g["thermo"]["symmetry_number"] == 2
    g1 = flat.Calculator("gxtb-v1.1", numbers, pos).hessian(pos)
    assert g1["frequencies"].shape == (3,) and g1["ir_intensities"].shape == (3,)
