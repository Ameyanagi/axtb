# SPDX-License-Identifier: LGPL-3.0-or-later
"""axtb: extended tight binding (GFN2-xTB; g-xTB v1.1 and 2.0.1) and GFN-FF in process, as a
program of the qc* stack: `qcdata.ProgramInput` in, `qcdata.ProgramOutput` out.

>>> import axtb
>>> from qcdata import ProgramInput, Structure
>>> water = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
>>> out = axtb.compute(ProgramInput(program="axtb", calctype="gradient", structure=water,
...                                 model={"method": "gxtb-v2.0.1"}))
>>> out.results.energy, out.results.gradient
>>> with axtb.Session(ProgramInput(program="axtb", calctype="gradient", structure=water,
...                                model={"method": "gfn2"})) as s:
...     frame = s.compute(water.geometry + 0.01)            # warm-started per lineage

- `compute(inp, ctx=None)`: one ``energy``, ``gradient`` or ``hessian`` calculation (a Hessian
  with its vibrational analysis and thermochemistry); failures raise `qcdata.ProgramError` with
  the failed output attached.
- `Session(inp, *, calctypes=None, ctx=None)`: one molecule at many geometries on one engine,
  with lineages, seeds and uncommitted requests (qccompute's session protocol).
- `program`: the object registered under the ``qccompute.programs`` entry points (no qccompute
  import at run time).
- `__version__`: axtb's version as pip reports it.

Keywords, outputs, the record, errors, threads and accuracy targets: ``docs/python.md``. Units
are atomic (bohr, hartree), as qcdata's conventions say.

The x86-64 Linux builds need an x86-64-v3 CPU (AVX2, FMA, BMI2, ...). ``import axtb`` checks the
CPU flags in ``/proc/cpuinfo`` before it loads the extension module and raises ``ImportError``
on a CPU without them. Where ``/proc/cpuinfo`` cannot be read or lists no flags (other
operating systems, some emulators) the check is skipped, and on such a CPU the import would end
the process with "Illegal instruction". User-mode QEMU shows the host's ``/proc/cpuinfo``, not
the emulated CPU's. A build with ``RUSTFLAGS="-C target-cpu=x86-64"`` runs on any x86-64 CPU and
skips the check.
"""

import platform as _platform
import sys as _sys

# The CLI's message (axtb_driver's CPU_MISMATCH in crates/driver/src/lib.rs; a test compares them).
_CPU_MISMATCH = (
    "this axtb was built for x86-64-v3 CPUs (AVX2, FMA, BMI2: Intel since Haswell (2013), AMD "
    "since Excavator (2015); not many low-end Intel Atom, Celeron and Pentium N-series CPUs), and "
    "this CPU (or its operating system) does not support them. Build axtb from source for any "
    'x86-64 CPU with RUSTFLAGS="-C target-cpu=x86-64" (CLI: cargo install axtb-cli; '
    "Python: pip install --no-binary axtb axtb)"
)

# x86-64-v3 as Linux lists it in /proc/cpuinfo; LZCNT is "abm". The kernel lists AVX (and
# everything built on it) only when it saves the AVX registers (OSXSAVE, which it does not list).
_V3_FLAGS = ("avx", "avx2", "fma", "bmi1", "bmi2", "f16c", "movbe", "xsave")


def _missing_cpu_flags(lines):
    """The x86-64-v3 flags missing from the first ``flags`` line of ``/proc/cpuinfo``'s
    ``lines``, or ``None`` when there is no such line."""
    for line in lines:
        key, sep, value = line.partition(":")
        if sep and key.strip() == "flags":
            flags = set(value.split())
            missing = [f for f in _V3_FLAGS if f not in flags]
            if not flags & {"abm", "lzcnt"}:
                missing.append("lzcnt")
            return missing
    return None


def _check_cpu():
    """``ImportError`` on a Linux x86-64 CPU without what this build assumes, before the
    extension module (and its first AVX instruction) loads."""
    try:
        from ._build import X86_64_V3  # written by build.rs
    except ImportError:
        X86_64_V3 = True  # the default build
    if not X86_64_V3 or not _sys.platform.startswith("linux") or _platform.machine() != "x86_64":
        return
    try:
        with open("/proc/cpuinfo", encoding="ascii", errors="replace") as f:
            missing = _missing_cpu_flags(f)
    except OSError:
        return
    if missing:
        raise ImportError(_CPU_MISMATCH)


_check_cpu()

from ._program import program  # noqa: E402  (after the CPU check)
from ._qc import VERSION as __version__  # noqa: E402
from ._qc import compute  # noqa: E402
from ._session import Session  # noqa: E402

__all__ = ["Session", "__version__", "compute", "program"]
