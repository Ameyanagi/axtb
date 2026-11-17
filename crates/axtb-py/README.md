# axtb

Extended tight binding (GFN2-xTB; g-xTB v1.1 and 2.0.1) and the GFN-FF force field, in process,
as a program of the qc* stack: a `qcdata.ProgramInput` goes in, a `qcdata.ProgramOutput` comes
out. Energies, analytic gradients, Hessians with vibrational analysis and thermochemistry,
implicit solvation (ALPB/GBSA for GFN2-xTB, and for g-xTB through a composite), orbitals,
populations and multipoles, bond orders. Molecules only. The engine is
[axtb](https://github.com/atomsforhumanity/axtb), written in Rust.

In a reference program's conventions axtb reproduces that program's numbers to the last
digits (xtb 6.8.0 and 6.7.1, tblite 0.8.0, the g-xTB 1.1.0 and 2.0.1 programs); its default
conventions use CODATA 2022 constants and correct the program bugs found on the way. A
calculation runs on exactly the threads it is given and never spawns more, so many can run
side by side.

```python
import axtb
from qcdata import ProgramInput, Structure

water = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])  # bohr
inp = ProgramInput(program="axtb", calctype="gradient", structure=water,
                   model={"method": "gxtb-v2.0.1"}, keywords={"threads": 4})
out = axtb.compute(inp)                      # raises qcdata.ProgramError on failure
out.results.energy, out.results.gradient

with axtb.Session(inp) as s:                 # one molecule at many geometries, warm-started
    frame = s.compute(water.geometry + 0.01)
    frame.results.energy, frame.record["warm_start"]

hess = axtb.compute(ProgramInput(program="axtb", calctype="hessian", structure=water, model={"method": "gfn2"}))
va = hess.results.vibrational_analysis      # frequencies (cm⁻¹), normal modes, IR intensities
hess.results.thermochemistry[-1].total.gibbs
```

qccompute finds axtb through its entry point (`axtb.program`), with no import of qccompute
in axtb: `qccompute.compute(inp)`, `qccompute.session(inp)` and `python -m qccompute.serve
axtb` run it like any other program.

**Install.** This is a pre-release: `pip install --pre axtb` (or `pip install axtb==0.2.0a1`).
It needs qcdata 0.24 (a pre-release, which pip installs with it) and numpy. The wheels are for
x86-64 Linux with glibc 2.28 or newer and Python 3.11 or newer (one abi3 wheel), and need an
x86-64-v3 CPU (AVX2, FMA: Intel since Haswell, AMD since Excavator); `import axtb` checks the
CPU and raises `ImportError` on one without these instructions. On other platforms pip builds
the package from the source distribution, which needs Rust 1.89 or newer.

**Documentation.** [docs/python.md](https://github.com/atomsforhumanity/axtb/blob/main/docs/python.md)
is the reference of the Python interface: keywords, outputs, the record, errors, sessions,
threads, accuracy targets. The [README](https://github.com/atomsforhumanity/axtb#readme)
describes the methods, conventions and settings;
[docs/validation.md](https://github.com/atomsforhumanity/axtb/blob/main/docs/validation.md) has
the agreement with each program; the
[changelog](https://github.com/atomsforhumanity/axtb/blob/main/CHANGELOG.md) lists what each
release contains, with the migration from axtb 0.1's `singlepoint` and `Calculator`.

**License.** LGPL-3.0-or-later (`COPYING.LESSER`, with `COPYING`), g-xTB included: g-xTB itself
is GPL-3.0 (grimme-lab/g-xtb), and its authors gave explicit permission to distribute axtb's
implementation, with the g-xTB parameters, under LGPL-3.0-or-later. `LICENSE` lists the parts;
`NOTICE` quotes the permission and credits the upstream work (xtb, tblite, dftd4, mctc-lib,
multicharge, g-xTB).
