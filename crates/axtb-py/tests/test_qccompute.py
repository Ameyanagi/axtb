# SPDX-License-Identifier: LGPL-3.0-or-later
"""axtb under qccompute: the registered program found by name, the runner's compute and sessions,
the grant, identities, lists with failures, and a worker process (``python -m qccompute.serve
axtb``) giving the in-process numbers."""

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
qccompute = pytest.importorskip("qccompute", minversion="0.23.0a4")
from qcdata import ProgramError, ProgramErrors, ProgramInput, Request, Structure  # noqa: E402

import axtb  # noqa: E402

WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])


def inp(method="gfn2", calctype="gradient", keywords=None):
    return ProgramInput(
        program="axtb", calctype=calctype, structure=WATER, model={"method": method}, keywords=keywords or {}
    )


def test_entry_point():
    from qccompute.program import available_programs, capabilities_of, find_program

    assert "axtb" in available_programs()
    program = find_program("axtb")
    assert program is axtb.program
    caps = capabilities_of(program)
    assert caps.sessions and caps.in_process and caps.thread_safe and not caps.batches
    assert set(caps.methods) == {"gfn2", "gfnff", "gxtb-v1.1", "gxtb-v2.0.1"}


def test_numerics():
    p = axtb.program
    assert p.numerics(inp()) == "gfn2:3"
    assert p.numerics(inp("GFN2-xTB", keywords={"solvent": {"name": "water"}})) == "gfn2:3"
    assert (
        p.numerics(inp("gxtb-v2.0.1", keywords={"solvent": {"name": "water", "correction": "gfn2"}}))
        == "gxtb-v2.0.1:3+gfn2:3"
    )
    assert p.numerics(inp("gxtb-v1.1", keywords={"conventions": "gxtb"})) == "gxtb-v1.1:5"
    assert p.numerics(inp("gfn3")) == f"axtb {axtb.__version__}"
    from qcdata import ProgramSpec

    assert p.numerics(ProgramSpec(program="axtb", calctype="gradient")) == f"axtb {axtb.__version__}"


def test_runner_compute_and_identities():
    out = qccompute.compute(inp())
    assert out.success and out.execution.identities == {"axtb": "gfn2:3"}
    assert out.execution.threads == 1
    direct = axtb.compute(inp(keywords={"threads": 1}))
    assert out.results.energy == direct.results.energy
    assert np.array_equal(out.results.gradient, direct.results.gradient)


def test_runner_list_with_a_failure():
    bad = inp(keywords={"scf": {"etmp": 1}})
    with pytest.raises(ProgramErrors) as e:
        qccompute.compute([inp(), bad, inp("gfnff")])
    outs = e.value.results
    assert [o.success for o in outs] == [True, False, True]
    assert outs[1].error.code == "unknown_key"
    with pytest.raises(ProgramError) as e1:
        qccompute.compute(bad)
    assert e1.value.output is not None and e1.value.output.error.code == "unknown_key"


def test_runner_grant():
    runner = qccompute.Runner(threads=2)
    # qccompute 0.23.0 lowers a threads keyword above the runner's count and says so
    with pytest.warns(RuntimeWarning, match="threads=8"):
        out = runner.compute(inp(keywords={"threads": 8}))
    assert [w.code for w in out.warnings] == ["resource_capped"]
    assert out.resolved.settings["threads"] == 2


def test_runner_session():
    x = np.asarray(WATER.geometry)
    with qccompute.session(inp()) as s:
        f0 = s.compute(x, lineage="a")
        f1 = s.compute(x + 0.01, lineage="a")
        frames = s.compute_many([Request(positions=x, lineage="b"), Request(positions=x + 0.01, lineage="b")])
        assert f1.record["warm_start"] is True
        assert frames[1].results.energy == f1.results.energy
        assert s.lineages() is not None
    assert f0.record["numerics"] == "gfn2:3"


def test_worker_round_trip():
    runner = qccompute.Runner(threads=2, isolation="process")
    try:
        out = runner.compute(inp(keywords={"threads": 1}))
        direct = axtb.compute(inp(keywords={"threads": 1}))
        assert out.results.energy == direct.results.energy
        assert np.array_equal(out.results.gradient, direct.results.gradient)
        with runner.session(inp(keywords={"threads": 1})) as s:
            f = s.compute(np.asarray(WATER.geometry) + 0.01)
            g = s.compute(np.asarray(WATER.geometry) + 0.02)
            assert g.record["warm_start"] is True and f.results.energy is not None
        with pytest.raises(ProgramError) as e:
            runner.compute(inp(keywords={"scf": {"max_iter": 0}}))
        assert e.value.error.code == "invalid_value"
    finally:
        runner.close()
