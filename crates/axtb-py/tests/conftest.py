# SPDX-License-Identifier: LGPL-3.0-or-later
"""Shared setup of the Python tests: the repository's tests/ on the path, for tests/flat.py (a
flat 0.1 request as a qcdata.ProgramInput, and axtb 0.1's calls through axtb.compute and
axtb.Session for the tests written against them)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "tests"))
