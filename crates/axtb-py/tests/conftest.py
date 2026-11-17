# SPDX-License-Identifier: LGPL-3.0-or-later
"""Shared helpers of the Python tests."""

import importlib.metadata

import pytest


def old_qccompute(module_level=False):
    """qccompute for the tests of the qccodec codec (``axtb.qccodec``), which needs qccompute
    before 0.23; skips with the 0.23 alphas, whose programs register under
    ``qccompute.programs`` (the codec and these tests go with the old Python API)."""
    qccompute = pytest.importorskip("qccompute")
    version = importlib.metadata.version("qccompute")
    if tuple(int(x) for x in version.split(".")[:2]) >= (0, 23):
        pytest.skip(f"the qccodec codec needs qccompute < 0.23, not {version}", allow_module_level=module_level)
    return qccompute
