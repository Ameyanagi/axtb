# SPDX-License-Identifier: LGPL-3.0-or-later
"""axtb in a forked child: a child that never computed opens a session and computes, whether the
parent computed before the fork (threaded single points, a Hessian, a session) or not. Workflows
whose workers are forked (mepd) rely on it; an earlier vibrational analysis on a global thread
pool hung there."""

import os
import signal
import sys
import time
from pathlib import Path

import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402

pytestmark = pytest.mark.skipif(not hasattr(os, "fork") or sys.platform == "darwin", reason="needs os.fork")

ROOT = Path(__file__).resolve().parents[3]
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")


def work():
    for calctype, threads in (("gradient", 1), ("gradient", 4), ("hessian", 4)):
        for method in ("gxtb-v2.0.1", "gfn2"):
            axtb.compute(
                ProgramInput(
                    program="axtb",
                    calctype=calctype,
                    structure=CAFFEINE,
                    model={"method": method},
                    keywords={"threads": threads},
                )
            )
    inp = ProgramInput(
        program="axtb",
        calctype="gradient",
        structure=CAFFEINE,
        model={"method": "gxtb-v2.0.1"},
        keywords={"threads": 4},
    )
    with axtb.Session(inp) as s:
        for _ in range(3):
            s.compute()


def in_child(timeout=120.0):
    pid = os.fork()
    if pid == 0:  # the child: compute, then leave without running the parent's cleanup
        signal.alarm(int(timeout))
        try:
            work()
        except BaseException:
            os._exit(3)
        os._exit(0)
    start = time.monotonic()
    while time.monotonic() - start < timeout + 30:
        done, status = os.waitpid(pid, os.WNOHANG)
        if done:
            return os.waitstatus_to_exitcode(status)
        time.sleep(0.1)
    os.kill(pid, signal.SIGKILL)
    os.waitpid(pid, 0)
    return "hung"


def test_child_of_a_parent_that_never_computed():
    assert in_child() == 0


def test_child_of_a_parent_that_computed():
    work()
    assert in_child() == 0
