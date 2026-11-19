# SPDX-License-Identifier: LGPL-3.0-or-later
"""``axtb.compute``: every method and calctype, the qcdata fields filled and left out, the
partition of the energy, the vibrational analysis and thermochemistry, the record, threads and
grants, failures with their outputs, and the CLI's numbers bit for bit."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import ProgramError, ProgramInput, Structure  # noqa: E402

import axtb  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests"))
import flat  # noqa: E402

WATER = Structure(symbols=["O", "H", "H"], geometry=[[0, 0, 0.22], [0, 1.43, -0.89], [0, -1.43, -0.89]])
CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")
AXTB_BIN = ROOT / "target/release/axtb"


def inp(method="gfn2", calctype="energy", keywords=None, structure=WATER):
    return ProgramInput(program="axtb", calctype=calctype, structure=structure, model={"method": method},
                        keywords=keywords if keywords is not None else {"threads": 1})


def cli(req):
    p = subprocess.run([str(AXTB_BIN), "--json"], input=json.dumps(req), capture_output=True, text=True)
    return json.loads(p.stdout)


METHODS = ["gfn2", "gxtb-v1.1", "gxtb-v2.0.1", "gfnff"]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("calctype", ["energy", "gradient", "hessian"])
@pytest.mark.parametrize("structure", [WATER, CAFFEINE], ids=["water", "caffeine"])
def test_every_method_and_calctype(method, calctype, structure):
    if calctype == "hessian" and structure is CAFFEINE and method == "gxtb-v1.1":
        pytest.skip("numerical v1.1 Hessian of caffeine: slow, covered on water")
    out = axtb.compute(inp(method, calctype, structure=structure))
    r = out.results
    n = len(structure.symbols)
    assert out.success and r.provenance.program == "axtb" and r.provenance.program_version == axtb.__version__
    assert r.calcinfo_natoms == n
    assert abs(sum(r.energy_components.values()) - r.energy) <= 1e-12 * max(1.0, abs(r.energy))
    assert r.atomic_charges.shape == (n,) and len(r.dipole_moment) == 3
    assert r.atomic_charges_kind == ("eeq" if method == "gfnff" else "mulliken")
    assert (r.gradient is None) == (calctype == "energy")
    if calctype != "energy":
        assert r.gradient.shape == (n, 3)
    if calctype == "hessian":
        assert r.hessian.shape == (3 * n, 3 * n) and np.array_equal(r.hessian, r.hessian.T)
        va = r.vibrational_analysis
        assert va.masses.shape == (n,) and len(va.frequencies) == 3 * n - 6
        assert va.rotor_type == "nonlinear" and va.projection == "translations_rotations"
        assert len(r.thermochemistry) == 1
        t = r.thermochemistry[0]
        assert t.temperature == 298.15 and t.pressure == 101325.0 and t.electronic_energy == r.energy
        assert t.total.enthalpy - r.energy == pytest.approx(
            t.translational.enthalpy + t.rotational.enthalpy + t.vibrational.enthalpy, abs=1e-12)
        assert t.total.gibbs == pytest.approx(t.total.enthalpy - t.temperature * t.total.entropy, abs=1e-10)
        assert out.record["hessian"]["method"] in ("analytic", "numerical", "response")
        assert any(k.startswith("hessian:") for k in out.execution.timings)
    else:
        assert r.hessian is None and r.vibrational_analysis is None and not r.thermochemistry
    if method == "gfnff":
        assert r.scf_iterations is None and r.quadrupole_moment is None
        assert "topology" in out.record["gfnff"]
    else:
        assert r.scf_iterations > 0 and r.quadrupole_moment.shape == (3, 3)
        assert abs(np.trace(r.quadrupole_moment)) < 1e-10
        assert out.record["scf"]["iterations"] == r.scf_iterations
        assert "electronic" in out.record["energy_subtotals"]
        assert r.calcinfo_nbasis == r.calcinfo_nmo
    assert out.record["numerics"] == axtb._axtb.numerics(method, False)
    assert out.execution.threads == 1 and out.execution.timings["total"] > 0
    assert out.resolved.settings["threads"] == 1


@pytest.mark.parametrize("method,calctype,keywords", [
    ("gfn2", "hessian", {"freq": {"temperature": [200, 300]}}),
    ("gxtb-v2.0.1", "gradient", {"solvent": {"name": "water", "correction": "gfn2"}, "output": {"bond_orders": True}}),
    ("gxtb-v1.1", "energy", {"scf": {"etemp": 300}}),
    ("gfnff", "hessian", {}),
])
@pytest.mark.filterwarnings("error")
def test_outputs_validate(method, calctype, keywords):
    """Outputs are built without validation (axtb._models.build): they validate as they are."""
    from qcdata import ProgramOutput

    out = axtb.compute(inp(method, calctype, {"threads": 1, **keywords}))
    again = ProgramOutput.model_validate(out.model_dump())
    assert again == out
    assert type(out).__name__ == type(again).__name__
    assert out.model_dump_json() == again.model_dump_json()


def test_bitwise_the_cli():
    if not AXTB_BIN.exists():
        pytest.skip("the release CLI is not built")
    for method in METHODS:
        req = {"numbers": [8, 1, 1], "positions": WATER.geometry.tolist(), "method": method, "hessian": True}
        o = cli(req)
        r = axtb.compute(flat.program_input(req)).results
        assert r.energy == o["energy"] and np.array_equal(r.gradient, o["gradient"])
        assert np.array_equal(r.hessian, o["hessian"])
        assert np.array_equal(r.vibrational_analysis.frequencies, o["frequencies"])


def test_composite_and_its_record():
    kw = {"threads": 1, "solvent": {"name": "water", "correction": "gfn2"}}
    out = axtb.compute(inp("gxtb-v2.0.1", "gradient", kw))
    rec = out.record
    assert set(rec["composite"]) >= {"definition", "components", "parts"}
    assert rec["wavefunction"] == {"part": "gxtb"}
    assert rec["contingent"]["composite.retries"]["fired"] is False
    assert rec["numerics"] == "gxtb-v2.0.1:2+gfn2:2"
    assert abs(sum(out.results.energy_components.values()) - out.results.energy) < 1e-12


def test_several_temperatures_in_order():
    kw = {"threads": 1, "freq": {"temperature": [400, 200, 298.15], "symmetry_number": 2}}
    r = axtb.compute(inp("gfn2", "hessian", kw)).results
    assert [t.temperature for t in r.thermochemistry] == [400.0, 200.0, 298.15]
    assert r.vibrational_analysis.symmetry_number == 2
    assert r.vibrational_analysis.symmetry_number_source == "input"
    for t in r.thermochemistry:
        assert t.vibrational.enthalpy == t.zero_point_energy + (t.vibrational.enthalpy - t.zero_point_energy)


def test_analysis_off():
    out = axtb.compute(inp("gfn2", "hessian", {"threads": 1, "freq": {"analysis": "off"}}))
    assert out.results.vibrational_analysis is None and out.results.hessian is not None
    out = axtb.compute(inp("gfn2", "hessian", {"threads": 1, "freq": {"thermo": "off"}}))
    assert out.results.vibrational_analysis is not None and not out.results.thermochemistry


def test_bond_orders():
    r = axtb.compute(inp("gfn2", "energy", {"threads": 1, "output": {"bond_orders": True}})).results
    assert r.bond_orders.shape == (3, 3) and r.bond_orders_kind == "mayer"
    assert axtb.compute(inp()).results.bond_orders is None


def test_threads_and_grants():
    from qccompute.program import RunContext

    out = axtb.compute(inp(keywords={"threads": 4}), RunContext(threads=2))
    assert out.resolved.settings["threads"] == 2
    assert [w.code for w in out.warnings] == ["resource_capped"]
    # water: one thread whatever the grant (below the one-thread size), the same bits
    assert out.execution.threads == 1
    out = axtb.compute(inp(keywords={}), RunContext(threads=6))
    assert out.resolved.settings["threads"] == 6 and not out.warnings
    # the keyword is a maximum: a caffeine gradient uses 1 of 3 (axtb's cap), its Hessian all 3
    grad = axtb.compute(inp("gfn2", "gradient", {"threads": 3}, structure=CAFFEINE))
    assert grad.resolved.settings["threads"] == 3 and grad.execution.threads == 1
    hess = axtb.compute(inp("gfn2", "hessian", {"threads": 3}, structure=CAFFEINE))
    assert hess.execution.threads == 3
    # a direct call without the keyword: this process's CPUs, at most the cap
    import os

    direct = axtb.compute(inp("gfn2", "gradient", {}, structure=CAFFEINE))
    assert direct.resolved.settings["threads"] == len(os.sched_getaffinity(0)) and direct.execution.threads == 1


def test_refusals_carry_the_failed_output():
    with pytest.raises(ProgramError) as e:
        axtb.compute(inp(keywords={"etemp": 300}))
    out = e.value.output
    assert out is not None and not out.success and out.error.code == "unknown_key"
    assert out.results.provenance.program == "axtb" and out.resolved is None
    with pytest.raises(TypeError):
        axtb.compute({"program": "axtb"})


def test_driver_failures_carry_their_fields():
    caf = inp("gfn2", "gradient", {"threads": 1, "scf": {"max_iter": 2}}, structure=CAFFEINE)
    with pytest.raises(ProgramError) as e:
        axtb.compute(caf)
    err = e.value.error
    assert (err.kind.value, err.code) == ("convergence", "scf_not_converged")
    assert set(err.details) == {"iterations", "energy_change", "residual"}
    out = e.value.output
    assert out.resolved.settings["scf"]["max_iter"] == 2
    assert out.record["error"] == err.details and out.execution.threads == 1


def test_no_qccompute_import():
    code = ("import sys, axtb; from qcdata import ProgramInput, Structure; "
            "s = Structure(symbols=['H', 'H'], geometry=[[0, 0, 0], [0, 0, 1.4]]); "
            "axtb.compute(ProgramInput(program='axtb', calctype='energy', structure=s, model={'method': 'gfn2'})); "
            "print('qccompute' in sys.modules)")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert p.stdout.strip() == "False"
