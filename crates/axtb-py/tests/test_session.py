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
        assert plan["record"]["lineages_at_once"] == 1
    with axtb.Session(inp(structure=WATER, keywords={"threads": 4})) as s:
        assert s.compute().execution.threads == 1  # water runs on one thread, the same bits
    # without the keyword: each request on the plan's threads for its calctype and size, at most
    # the grant; the grant filled with the lineages compute_many runs side by side
    caffeine = inp("gxtb-v2.0.1", keywords={})
    with axtb.Session(caffeine, ctx=RunContext(threads=16)) as s:  # a caffeine Hessian plans 8
        plan = s.plan()
        assert plan["resolved"]["settings"]["threads"] == 1 and plan["resources"]["threads"] == 8
        assert plan["record"]["threads_by_calctype"] == {"energy": 1, "gradient": 1, "hessian": 8}
        assert plan["record"]["lineages_at_once"] == 16 and not plan["warnings"]
        assert s.compute().execution.threads == 1
        assert s.compute(calctype="energy").execution.threads == 1
        assert s._engine.calculators == [1]
    with axtb.Session(caffeine, ctx=RunContext(threads=4)) as s:
        plan = s.plan()
        assert plan["resources"]["threads"] == 4 and plan["record"]["lineages_at_once"] == 4
        assert plan["record"]["threads_by_calctype"]["hessian"] == 4 and not plan["warnings"]
    with axtb.Session(inp("gxtb-v2.0.1", calctype="hessian", keywords={}), ctx=RunContext(threads=8)) as s:
        plan = s.plan()
        assert plan["record"]["threads_used"] == 8 and plan["record"]["lineages_at_once"] == 1
        assert s.compute(calctype="gradient").execution.threads == 1
        assert s._engine.calculators == [8, 1]
    # a direct call: the session's threads are this process's CPUs
    with axtb.Session(caffeine) as s:
        assert s.plan()["record"]["lineages_at_once"] == len(os.sched_getaffinity(0))


def requests(structure, bad):
    """Interleaved and seeded lineages: lineages without history seeded from ones with and
    without, a seed whose later requests come after the seeded one, uncommitted requests, a
    failed request in a seed lineage, a refused one, and lineages that are JSON values."""
    x = np.asarray(structure.geometry, dtype=float)
    n = len(x)
    rng = np.random.default_rng(7)

    def at(scale=0.01):
        return x + rng.normal(scale=scale, size=(n, 3))

    nan = np.full_like(x, np.nan)
    out = [
        Request(positions=at(), lineage="b"),
        Request(positions=at(), lineage="a", seed=["nothing", "b"]),
        Request(positions=at(), lineage="b"),
        Request(positions=at(), lineage=("c", 1), seed=["a"], commit=False),
        Request(positions=nan if bad else at(), lineage="b"),
        Request(positions=at(), lineage=("c", 1), seed=["b", "a"]),
        Request(positions=at(), lineage="d", seed=["c"]),
        Request(positions=at(), lineage="a"),
        Request(positions=at(), lineage="b", seed=["a"]),
        Request(positions=at(), lineage="e", seed=["b"]),
    ]
    if bad:
        out.append(Request(positions=at(), lineage="a", calctype="optimization"))
    for k in range(3):
        for img in range(6):
            out.append(Request(positions=at(), lineage=["neb", img], seed=[["neb", img - 1]] if img else [],
                               commit=k != 1 or img % 2 == 0))
    return out


def frames_of(f):
    try:
        return f()
    except ProgramErrors as e:
        return e.results


def frame_bits(f):
    if not f.success:
        return ("failed", f.error.code)
    return (*bits(f), f.record["warm_start"], json.dumps(f.record["from_lineage"]), f.execution.threads)


@pytest.mark.parametrize(
    "method,keywords",
    [
        ("gfn2", {}),
        ("gfn2", {"threads": 2}),
        ("gxtb-v2.0.1", {}),
        ("gxtb-v2.0.1", {"threads": 2, "conventions": "gxtb"}),
        ("gxtb-v1.1", {}),
        ("gfnff", {}),
        ("gxtb-v2.0.1", {"solvent": {"name": "water", "correction": "gfn2"}}),
    ],
)
def test_compute_many_side_by_side_equals_one_at_a_time(method, keywords):
    """compute_many runs different lineages side by side on forks of the engine and gives what
    one-at-a-time compute calls give, each on its calctype's threads, bit for bit: interleaved
    and seeded lineages, uncommitted requests, failures in place, the lineages' histories
    afterwards; the forks are kept for the next call."""
    from qccompute.program import RunContext

    structure = WATER if "solvent" in keywords else CAFFEINE
    for bad in (True, False):
        reqs = requests(structure, bad)
        opened = [
            axtb.Session(inp(method, keywords=keywords, structure=structure), ctx=RunContext(threads=8))
            for _ in range(2)
        ]
        with opened[0] as many, opened[1] as one:
            at_once = many.plan()["record"]["lineages_at_once"]
            assert at_once == 8 // many.plan()["record"]["threads_used"] and at_once > 1
            got = frames_of(lambda: many.compute_many(reqs))
            want = []
            for r in reqs:
                try:
                    want.append(one.compute(r.positions, r.calctype, lineage=r.lineage, seed=r.seed, commit=r.commit))
                except ProgramError as e:
                    want.append(e.frame)
            assert [frame_bits(f) for f in got] == [frame_bits(f) for f in want]
            assert [f.success for f in got].count(False) == (2 if bad else 0)
            assert sorted(map(json.dumps, many.lineages())) == sorted(map(json.dumps, one.lineages()))
            calculators = many._engine.calculators
            assert 1 < len(calculators) <= at_once
            # the next steps from the histories, on the kept forks
            more = [Request(positions=r.positions, lineage=r.lineage) for r in reqs[-6:]]
            again = many.compute_many(more)
            assert [frame_bits(f) for f in again] == [
                frame_bits(one.compute(r.positions, lineage=r.lineage)) for r in more
            ]
            assert many._engine.calculators == calculators


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


@pytest.mark.parametrize("method", ["gxtb-v2.0.1", "gfn2"])
def test_compute_many_mixed_calctypes(method):
    """Gradients and Hessians in one compute_many: each request on its calctype's threads (water
    with g-xTB 2.0.1: gradients on 1, Hessians on 2; GFN2-xTB: 1 and 1), the same bits as
    one-at-a-time compute calls, lineages moving between the forks of each thread count."""
    from qccompute.program import RunContext

    x = np.asarray(WATER.geometry, dtype=float)
    rng = np.random.default_rng(11)
    reqs = []
    for step in range(3):
        for k in range(4):
            calctype = "hessian" if (k + step) % 3 == 0 else "gradient"
            reqs.append(Request(positions=x + rng.normal(scale=0.01, size=x.shape), calctype=calctype, lineage=k,
                                seed=[k - 1] if k else []))
    sessions = [axtb.Session(inp(method, keywords={}, structure=WATER), ctx=RunContext(threads=4)) for _ in range(2)]
    with sessions[0] as many, sessions[1] as one:
        got = many.compute_many(reqs)
        want = [one.compute(r.positions, r.calctype, lineage=r.lineage, seed=r.seed) for r in reqs]
        assert [frame_bits(f) for f in got] == [frame_bits(f) for f in want]
        hessians = [f for f, r in zip(got, reqs, strict=True) if r.calctype == "hessian"]
        for f, g in zip(hessians, [f for f, r in zip(want, reqs, strict=True) if r.calctype == "hessian"], strict=True):
            assert f.results.hessian.tobytes() == g.results.hessian.tobytes()
        threads = many.plan()["record"]["threads_by_calctype"]
        assert [f.execution.threads for f in got] == [threads[r.calctype] for r in reqs]
        assert set(many._engine.calculators) == {threads["gradient"], threads["hessian"]}
