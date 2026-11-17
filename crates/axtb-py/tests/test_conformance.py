# SPDX-License-Identifier: LGPL-3.0-or-later
"""qccompute's conformance checks (``qccompute.testing.check_program``) of ``axtb.program``,
per method."""

import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
testing = pytest.importorskip("qccompute.testing")
from qcdata import ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402

WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])


@pytest.mark.parametrize("method,keywords", [
    ("gfn2", {}),
    ("gfn2", {"conventions": "xtb", "solvent": {"name": "water"}}),
    ("gxtb-v2.0.1", {}),
    ("gxtb-v2.0.1", {"conventions": "gxtb", "solvent": {"name": "water", "correction": "gfn2"}}),
    ("gxtb-v1.1", {}),
    ("gfnff", {}),
    ("gfnff", {"accuracy": {"gradient": 1e-6}}),
])
def test_check_program(method, keywords):
    inputs = [ProgramInput(program="axtb", calctype=c, structure=WATER, model={"method": method}, keywords=keywords)
              for c in ("gradient", "energy", "hessian")]
    ran = testing.check_program(axtb.program, inputs)
    assert "energy components add up" in " ".join(ran)
