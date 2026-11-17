# SPDX-License-Identifier: LGPL-3.0-or-later
"""The private engine of ``axtb.compute`` and ``axtb.Session`` (``axtb._axtb``): one-off
calculations as the CLI runs them, lineages with their committed histories, seeds, uncommitted
requests, resets, errors with kinds and codes, a closed engine, and calls from several threads.
"""

import json
import threading
from pathlib import Path

import numpy as np
import pytest

import flat  # tests/flat.py (the CLI's response)
from axtb import _axtb

ROOT = Path(__file__).resolve().parents[3]
BOHR = 0.529177210544
SYMBOLS = {"H": 1, "C": 6, "N": 7, "O": 8}


def xyz(name):
    lines = (ROOT / "tests/structures" / name).read_text().split("\n")
    n = int(lines[0])
    atoms = [line.split() for line in lines[2 : 2 + n]]
    return [SYMBOLS[a[0]] for a in atoms], np.array([[float(x) / BOHR for x in a[1:4]] for a in atoms])


def request(**settings):
    return json.dumps({"numbers": [], "positions": [], **settings})


def task(**fields):
    return json.dumps({"numbers": [], "positions": [], **fields})


GRADIENT = task(gradient=True)


def steps(positions, sign=1.0):
    """The structure and two displacements of it (golden.py's trajectories)."""
    n = len(positions)
    d = np.array([[((7 * i + 3 * c) % 5) - 2 for c in range(3)] for i in range(n)], dtype=float)
    return [positions + sign * 0.01 * k * d for k in range(3)]


def bits(r):
    return (r["energy"], r["gradient"].tobytes(), r["charges"].tobytes(), r["scf"]["guess"], r["scf"]["iterations"])


@pytest.fixture(scope="module")
def caffeine():
    return xyz("caffeine.xyz")


@pytest.mark.skipif(not Path(flat.AXTB_BIN).exists(), reason="needs the release CLI")
def test_compute_once_is_the_cli_calculation(caffeine):
    numbers, pos = caffeine
    for method in ("gfn2", "gxtb-v2.0.1"):
        r = _axtb.compute_once(request(method=method, threads=1), numbers, pos, GRADIENT)
        cli = flat.cli({"numbers": numbers, "positions": pos, "method": method, "gradient": True, "threads": 1})
        assert r["energy"] == cli["energy"]
        assert np.array_equal(r["gradient"], cli["gradient"])
        assert r["threads"] == 1
        assert r["numerics"] == _axtb.numerics(method, False)


def test_lineages_interleaved_equal_sequential(caffeine):
    numbers, pos = caffeine
    for method, conventions in (("gfn2", "axtb"), ("gxtb-v2.0.1", "axtb"), ("gxtb-v1.1", "gxtb")):
        req = request(method=method, conventions=conventions, guess="previous", threads=1)
        a, b = steps(pos), steps(pos, -1.0)
        seq = []
        for traj in (a, b):
            e = _axtb.Engine(req, numbers, pos)
            seq.append([bits(e.compute(x, GRADIENT, "null")) for x in traj])
        e = _axtb.Engine(req, numbers, pos)
        inter = [[], []]
        for x, y in zip(a, b):
            inter[0].append(bits(e.compute(x, GRADIENT, '"a"')))
            inter[1].append(bits(e.compute(y, GRADIENT, '"b"')))
        assert inter == seq, (method, conventions)
        assert sorted(e.lineages()) == ['"a"', '"b"']


def test_seed_commit_and_reset(caffeine):
    numbers, pos = caffeine
    e = _axtb.Engine(request(method="gfn2", threads=1), numbers, pos)
    a = steps(pos)
    r0 = e.compute(a[0], GRADIENT, '"ts"')
    assert r0["history"] is None and r0["scf"]["guess"] == "eeq"
    # a new lineage seeded from "ts" starts from its history; the first seed with history wins
    r1 = e.compute(a[1], GRADIENT, '"fw"', seeds=['"none"', '"ts"'])
    assert r1["history"] == '"ts"' and r1["scf"]["guess"] == "previous"
    # its own history wins over its seeds
    r2 = e.compute(a[2], GRADIENT, '"fw"', seeds=['"ts"'])
    assert r2["history"] == '"fw"'
    # commit=False stores nothing: two displacements from the centre start from the centre alike
    ref = _axtb.Engine(request(method="gfn2", threads=1), numbers, pos)
    ref.compute(a[0], GRADIENT, "null")
    want = bits(ref.compute(a[2], GRADIENT, "null"))
    c = _axtb.Engine(request(method="gfn2", threads=1), numbers, pos)
    c.compute(a[0], GRADIENT, "null")
    c.compute(a[1], GRADIENT, "null", commit=False)
    assert bits(c.compute(a[2], GRADIENT, "null", commit=False)) == want
    assert c.lineages() == ["null"]
    # reset(key) and reset(): the next request starts fresh
    e.reset('"fw"')
    assert sorted(e.lineages()) == ['"ts"']
    r = e.compute(a[2], GRADIENT, '"fw"')
    assert r["history"] is None and r["scf"]["guess"] == "eeq"
    e.reset()
    assert e.lineages() == []
    assert e.compute(a[2], GRADIENT, '"ts"')["history"] is None


def test_a_fresh_request_after_an_uncommitted_one_starts_fresh(caffeine):
    numbers, pos = caffeine
    e = _axtb.Engine(request(method="gxtb-v2.0.1", threads=1), numbers, pos)
    fresh = bits(_axtb.Engine(request(method="gxtb-v2.0.1", threads=1), numbers, pos).compute(pos, GRADIENT, "null"))
    e.compute(pos + 0.01, GRADIENT, '"x"', commit=False)
    assert bits(e.compute(pos, GRADIENT, '"y"')) == fresh


def test_scf_accuracy_per_call(caffeine):
    numbers, pos = caffeine
    tight = _axtb.compute_once(request(method="gxtb-v2.0.1", accuracy=0.01, threads=1), numbers, pos, GRADIENT)
    e = _axtb.Engine(request(method="gxtb-v2.0.1", threads=1), numbers, pos)
    loose = e.compute(pos, GRADIENT, '"a"')
    r = e.compute(pos, GRADIENT, '"b"', accuracy=0.01)
    assert bits(r) == bits(tight)
    # back to the engine's own accuracy
    assert bits(e.compute(pos, GRADIENT, '"c"', accuracy=1.0)) == bits(loose)
    # one-off: the accuracy set after the method is built (GFN2-xTB's screening stays)
    once = _axtb.compute_once(request(method="gxtb-v2.0.1", threads=1), numbers, pos, GRADIENT, accuracy=0.01)
    assert bits(once) == bits(tight)


def test_errors_carry_kind_code_and_fields(caffeine):
    numbers, pos = caffeine
    with pytest.raises(_axtb.Error) as e:
        _axtb.compute_once(request(method="gfn2", max_iter=2, threads=1), numbers, pos, GRADIENT)
    assert (e.value.kind, e.value.code) == ("scf_not_converged", "scf_not_converged")
    assert set(e.value.fields) == {"iterations", "energy_change", "residual"}
    with pytest.raises(_axtb.Error) as e:
        _axtb.Engine(request(method="gxtb-v1.1"), [8, 1, 1, 60], np.zeros((4, 3)) + np.arange(4)[:, None])
    assert (e.value.kind, e.value.code, e.value.fields) == ("unsupported_element", "unsupported_element", {"z": 60})
    with pytest.raises(_axtb.Error) as e:
        _axtb.compute_once(request(method="gfnff"), numbers, pos, task(bond_orders=True))
    assert (e.value.kind, e.value.code) == ("not_applicable", "not_applicable")
    with pytest.raises(_axtb.Error) as e:
        _axtb.compute_once(request(method="gfn3"), numbers, pos, GRADIENT)
    assert (e.value.kind, e.value.code) == ("input", "invalid_input")
    assert e.value.message.startswith("invalid input: unknown method")
    # a failed request leaves the lineage's history untouched
    eng = _axtb.Engine(request(method="gfn2", threads=1), numbers, pos)
    good = eng.compute(pos, GRADIENT, "null")
    with pytest.raises(_axtb.Error):
        eng.compute(np.full_like(pos, np.nan), GRADIENT, "null")
    assert eng.lineages() == ["null"]
    again = eng.compute(pos, GRADIENT, "null")
    assert again["history"] == "null" and again["energy"] == pytest.approx(good["energy"], abs=1e-9)


def test_closed(caffeine):
    numbers, pos = caffeine
    e = _axtb.Engine(request(method="gfnff"), numbers, pos)
    assert e.compute(pos, GRADIENT, "null")["history"] is None
    assert e.lineages() == []  # GFN-FF keeps no wavefunction
    e.close()
    assert e.closed
    with pytest.raises(_axtb.Error) as x:
        e.compute(pos, GRADIENT, "null")
    assert (x.value.kind, x.value.code) == ("closed", "session_closed")


def test_threads_share_one_engine(caffeine):
    numbers, pos = caffeine
    e = _axtb.Engine(request(method="gfn2", threads=1), numbers, pos)
    want = bits(_axtb.compute_once(request(method="gfn2", threads=1), numbers, pos, GRADIENT))
    out = []

    def work(i):
        for _ in range(3):
            out.append(bits(e.compute(pos, GRADIENT, json.dumps(i), commit=False)))

    ts = [threading.Thread(target=work, args=(i,)) for i in range(4)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    assert out == [want] * 12


def test_name_helpers():
    assert _axtb.canonical("method", "GFN2-xTB") == "gfn2"
    assert _axtb.canonical("conventions", "xtb-6.8.0", "gfn2") == "xtb"
    assert _axtb.canonical("solvent", "H2O") == "water"
    assert _axtb.canonical_solvent("H2O", "gbsa") == "water"
    assert _axtb.canonical_solvent("ammonia", "alpb") is None or isinstance(_axtb.canonical_solvent("ammonia", "alpb"), str)
    assert "water" in _axtb.solvents("alpb")
    assert _axtb.supported_elements("gfn2") == list(range(1, 87))
    assert _axtb.numerics("gxtb-v2.0.1", True) == "gxtb-v2.0.1:1+gfn2:1"
    assert _axtb.single_thread_atoms("gfnff", False) == 100
    assert _axtb.ANALYSIS_NUMERICS == "analysis:1"
