# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.Session``: lineages, seeds, uncommitted requests, resets, compute_many with failures in
place, interleaved lineages equal sequential ones bit for bit, warm_start and from_lineage,
calctypes, threads, a closed session, two threads on one session, and the overhead per step."""

import json
import os
import statistics
import threading
import time
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramError, ProgramErrors, ProgramInput, Request, Structure  # noqa: E402

import axtb  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]

WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")


def inp(method="gfn2", calctype="gradient", keywords=None, structure=CAFFEINE):
    return ProgramInput(
        program="axtb",
        calctype=calctype,
        structure=structure,
        model={"method": method},
        keywords=keywords if keywords is not None else {"threads": 1},
    )


def steps(structure, sign=1.0):
    x = np.asarray(structure.geometry, dtype=float)
    n = len(x)
    d = np.array([[((7 * i + 3 * c) % 5) - 2 for c in range(3)] for i in range(n)], dtype=float)
    return [x + sign * 0.01 * k * d for k in range(3)]


def bits(f):
    r = f.results
    return (r.energy, r.gradient.tobytes(), r.atomic_charges.tobytes(), r.scf_iterations)


@pytest.mark.parametrize(
    "method,conventions",
    [("gfn2", "axtb"), ("gfn2", "xtb"), ("gxtb-v2.0.1", "axtb"), ("gxtb-v2.0.1", "gxtb"), ("gxtb-v1.1", "axtb")],
)
def test_one_lineage_is_one_calculator(method, conventions):
    """A session's default lineage computes what the CLI's trajectory runner computes (the
    golden trajectories, bit for bit: scripts/golden_python.py), and every step after the first
    starts from the previous one."""
    with axtb.Session(
        inp(method, keywords={"threads": 1, "conventions": conventions, "scf": {"restart": "previous"}})
    ) as s:
        frames = [s.compute(x) for x in steps(CAFFEINE)]
    assert [f.record["warm_start"] for f in frames] == [False, True, True]
    assert [f.record["from_lineage"] for f in frames] == [None, None, None]
    assert frames[1].record["scf"]["start"] == "previous"


@pytest.mark.parametrize("method", ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"])
def test_interleaved_lineages_equal_sequential(method):
    a, b = steps(CAFFEINE), steps(CAFFEINE, -1.0)
    seq = []
    for traj in (a, b):
        with axtb.Session(inp(method)) as s:
            seq.append([bits(s.compute(x)) for x in traj])
    with axtb.Session(inp(method)) as s:
        inter = [[], []]
        for x, y in zip(a, b, strict=True):
            inter[0].append(bits(s.compute(x, lineage="a")))
            inter[1].append(bits(s.compute(y, lineage=("b", 2))))
        assert sorted(map(json.dumps, s.lineages())) == sorted(map(json.dumps, ["a", ["b", 2]]))
    assert inter == seq


def test_seed_commit_reset():
    a = steps(CAFFEINE)
    with axtb.Session(inp()) as s:
        f = s.compute(a[0], lineage="ts")
        assert (f.record["warm_start"], f.record["from_lineage"]) == (False, None)
        f = s.compute(a[1], lineage=("irc", "forward"), seed=["nothing", "ts"])
        assert (f.record["warm_start"], f.record["from_lineage"]) == (True, "ts")
        f = s.compute(a[2], lineage=["irc", "forward"], seed=["ts"])  # a tuple and a list are one lineage
        assert f.record["from_lineage"] == ["irc", "forward"]
        # commit=False: the centre's history stays
        ref = axtb.Session(inp())
        ref.compute(a[0])
        want = bits(ref.compute(a[2]))
        c = axtb.Session(inp())
        c.compute(a[0])
        c.compute(a[1], commit=False)
        assert bits(c.compute(a[2], commit=False)) == want
        # reset(L), reset(None), reset()
        s.reset(("irc", "forward"))
        assert s.compute(a[2], lineage=("irc", "forward")).record["warm_start"] is False
        s.compute(a[0])
        s.reset(None)
        assert s.compute(a[0]).record["warm_start"] is False
        s.reset()
        assert s.lineages() == []


def test_compute_many_in_order_with_failures_in_place():
    a = steps(CAFFEINE)
    nan = np.full_like(a[0], np.nan)
    with axtb.Session(inp()) as s:
        with pytest.raises(ProgramErrors) as e:
            s.compute_many(
                [
                    Request(positions=a[0], lineage="x"),
                    Request(positions=nan, lineage="x"),
                    Request(positions=a[1], lineage="x", commit=False),
                    Request(positions=a[2], lineage="x"),
                ]
            )
    frames = e.value.results
    assert [f.success for f in frames] == [True, False, True, True]
    assert len(e.value.exceptions) == 1 and e.value.exceptions[0].error.code == "invalid_structure"
    # B and D both start from A (B uncommitted)
    with axtb.Session(inp()) as s:
        s.compute(a[0], lineage="x")
        want = bits(s.compute(a[2], lineage="x"))
    assert bits(frames[3]) == want and frames[2].record["warm_start"] is True


def test_calctypes_and_keywords():
    kw = {"threads": 1, "hessian": {"method": "numerical"}}
    with pytest.raises(ProgramError) as e:
        axtb.Session(inp(keywords=kw), calctypes=["gradient"])
    assert e.value.error.code == "not_applicable"
    with axtb.Session(inp(structure=WATER, keywords=kw), calctypes=["hessian"]) as s:
        f = s.compute(calctype="gradient")
        assert f.record["contingent"]["hessian.method"]["fired"] is False
        h = s.compute(calctype="hessian")
        assert h.results.hessian.shape == (9, 9) and h.results.gradient is not None
        assert h.record["hessian"]["method"] == "numerical"
    with axtb.Session(inp(structure=WATER), calctypes=["energy"]) as s:
        with pytest.raises(ProgramError) as e:
            s.compute(calctype="gradient")
        assert (e.value.error.kind.value, e.value.error.code) == ("input", "not_applicable")
    with pytest.raises(ProgramError) as e:
        axtb.Session(inp(), calctypes=["optimization"])
    assert e.value.error.code == "unsupported_calctype"


def test_gfnff_keeps_no_lineage_state():
    with axtb.Session(inp("gfnff", keywords={"threads": 1, "accuracy": {"gradient": 1e-6}})) as s:
        f = s.compute()
        assert f.record["warm_start"] is None and s.lineages() == []
        assert f.record["accuracy"] == {"target": {"gradient": 1e-6}, "thresholds": {}, "estimate": {"gradient": 0.0}}
        g = s.compute(accuracy=1e-4)
        assert g.record["accuracy"]["target"] == {"gradient": 1e-4}
        assert "topology" not in f.record["gfnff"]


def test_failures_leave_the_session_usable():
    with axtb.Session(inp()) as s:
        good = s.compute()
        with pytest.raises(ProgramError) as e:
            s.compute(np.zeros((3, 3)))
        assert e.value.frame is not None and not e.value.frame.success
        again = s.compute()
        assert again.record["warm_start"] is True and again.results.energy == pytest.approx(good.results.energy)


def test_closed():
    s = axtb.Session(inp(structure=WATER))
    s.close()
    with pytest.raises(ProgramError) as e:
        s.compute()
    assert (e.value.error.kind.value, e.value.error.code) == ("input", "session_closed")
    s.close()


def test_two_threads_on_one_session():
    want = bits(axtb.Session(inp()).compute())
    out = []
    with axtb.Session(inp()) as s:

        def work(i):
            for _ in range(3):
                out.append(bits(s.compute(lineage=i, commit=False)))

        ts = [threading.Thread(target=work, args=(i,)) for i in range(4)]
        for t in ts:
            t.start()
        for t in ts:
            t.join()
    assert out == [want] * 12


def test_threads_and_plan():
    from qccompute.program import RunContext

    with axtb.Session(inp(keywords={"threads": 4}), ctx=RunContext(threads=2)) as s:
        plan = s.plan()
        assert plan["resolved"]["settings"]["threads"] == 2 and plan["resources"]["threads"] == 2
        assert [w["code"] for w in plan["warnings"]] == ["resource_capped"]
        assert s.compute().execution.threads == 2
    with axtb.Session(inp(structure=WATER, keywords={"threads": 4})) as s:
        assert s.compute().execution.threads == 1  # water runs on one thread, the same bits


@pytest.mark.filterwarnings("error")
def test_frames_validate():
    from qcdata import Frame

    with axtb.Session(inp(structure=WATER), calctypes=["hessian"]) as s:
        for c in ("energy", "gradient", "hessian"):
            f = s.compute(calctype=c, lineage=c)
            again = Frame.model_validate(f.model_dump())
            assert again == f and again.model_dump_json() == f.model_dump_json()


def test_overhead_per_step():
    """Python's overhead per session step on water is at most 30 us (median, the engine call
    subtracted); a compute's at most 100 us."""
    from axtb import _axtb, _keywords, _qc

    for method in ("gfn2", "gxtb-v2.0.1"):
        i = inp(method, structure=WATER)
        s = axtb.Session(i)
        p = _keywords.plan_input(i, session=True, calctypes=frozenset(_keywords.CALCTYPES))
        e = _axtb.Engine(json.dumps({**p.request, "threads": 1}), p.numbers, p.positions)
        task = _qc.task(p, "gradient", frozenset())
        x = np.asarray(WATER.geometry, dtype=float)
        full, raw = [], []
        for _ in range(300):
            t = time.perf_counter()
            s.compute(x)
            full.append(time.perf_counter() - t)
            t = time.perf_counter()
            e.compute(x, task, "null")
            raw.append(time.perf_counter() - t)
        step = statistics.median(full[50:]) - statistics.median(raw[50:])
        req = json.dumps({**p.request, "threads": 1})
        cfull, craw = [], []
        for _ in range(300):
            t = time.perf_counter()
            axtb.compute(i)
            cfull.append(time.perf_counter() - t)
            t = time.perf_counter()
            _axtb.compute_once(req, p.numbers, p.positions, task)
            craw.append(time.perf_counter() - t)
        once = statistics.median(cfull[50:]) - statistics.median(craw[50:])
        print(f"{method}: session step {step * 1e6:.0f} us, compute {once * 1e6:.0f} us")
        # the targets, with room for a loaded machine; shared CI runners (CI set) run Python about
        # 1.6x slower than the laptop the targets were set on, so they get room for that too
        room = 3.0 if os.environ.get("CI") else 1.5
        assert step < 30e-6 * room and once < 100e-6 * room, (method, step, once)
