# SPDX-License-Identifier: LGPL-3.0-or-later
"""The release workflow's smoke test of an installed axtb (a wheel or a build from the sdist):
one gradient per method through ``axtb.compute``, the qccompute entry point, and no qccompute
import at run time. Run it outside the repository, so that ``import axtb`` finds the install."""

import sys
from importlib.metadata import entry_points

import numpy as np
from qcdata import ProgramInput, Structure

import axtb

water = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
for method in ("gfn2", "gxtb-v1.1", "gxtb-v2.0.1", "gfnff"):
    out = axtb.compute(ProgramInput(program="axtb", calctype="gradient", structure=water, model={"method": method}))
    assert out.success and np.asarray(out.results.gradient).shape == (3, 3), out
    print(axtb.__version__, method, out.results.energy, out.record["numerics"])
assert [e.value for e in entry_points(group="qccompute.programs", name="axtb")] == ["axtb:program"]
assert "qccompute" not in sys.modules
