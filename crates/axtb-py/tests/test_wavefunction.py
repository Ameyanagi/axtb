# SPDX-License-Identifier: LGPL-3.0-or-later
"""The wavefunction outputs (design 5.2-5.4) for every SCF method and conventions: the orbital
levels and core electrons, the orbitals with their basis set, the overlap and the Hamiltonian
(C^T S C = 1, H C = S C e), an overlap and atomic dipoles computed independently from the basis
set (Obara-Saika, tests/molden_check.py), the atomic multipoles adding up to the molecular
moments, shell populations, restricted and unrestricted orbitals, refusals, and that requesting
an output changes no other number."""

import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("qcdata", minversion="0.24.0a5")
from qcdata import MOLDEN_CONVENTION, ProgramError, ProgramInput, ProgramOutput, Structure  # noqa: E402

import axtb  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests"))
import molden_check as mc  # noqa: E402

CAFFEINE = Structure.open(ROOT / "tests/structures/caffeine.xyz")
CU = Structure.open(ROOT / "tests/structures/cu-nh3-4.xyz").model_copy(update={"charge": 2, "multiplicity": 2})
SETUPS = [("gfn2", c) for c in ("axtb", "xtb-6.8.0", "xtb-6.7.1", "tblite-0.8.0")] + [
    ("gxtb-v2.0.1", "axtb"),
    ("gxtb-v2.0.1", "gxtb-2.0.1"),
    ("gxtb-v1.1", "axtb"),
    ("gxtb-v1.1", "gxtb-1.1.0"),
    ("composite", "axtb"),
    ("composite", "gxtb-2.0.1"),
]
OUT = {"orbitals": True, "overlap": True, "hamiltonian": True}


def inp(method, conventions, structure, output=None, calctype="energy"):
    kw = {"threads": 1, "conventions": conventions}
    if method == "composite":
        method = "gxtb-v2.0.1"
        kw["solvent"] = {"name": "water", "correction": "gfn2"}
    if output:
        kw["output"] = output
    return ProgramInput(program="axtb", calctype=calctype, structure=structure, model={"method": method}, keywords=kw)


def program_hamiltonian(method, conventions):
    return not (method in ("gxtb-v2.0.1", "composite") and conventions != "axtb")


def independent(wf, structure):
    """The overlap and the dipole integrals of the basis set, computed apart from axtb (Molden
    order and normalization), and the wavefunction in that convention."""
    m = wf.with_convention(MOLDEN_CONVENTION)
    pos = np.asarray(structure.geometry, dtype=float)
    mol = {
        "atoms": [(s, pos[i]) for i, s in enumerate(structure.symbols)],
        "shells": [
            (sh.center, sh.angular_momentum, np.asarray(sh.exponents), np.asarray(sh.coefficients))
            for sh in m.basis.shells
        ],
        "pure": {2: True, 3: True, 4: True},
    }
    b = mc.Basis(mol, prims="raw")
    return m, b


def dipole_integrals(b, k):
    """D_nu,mu = <nu| (r - R_mu)_k |mu>, the integrals about the atom of the column function."""
    out = np.zeros((b.nbf, b.nbf))
    i0 = 0
    for sa, ta in zip(b.shells, b.trafo, strict=True):
        j0 = 0
        for sb, tb in zip(b.shells, b.trafo, strict=True):
            (ca, la, aa, wa), (cb, lb, ab, wb) = sa, sb
            tx = [mc.s1d(la, lb + 1, aa, ab, ca[d], cb[d]) for d in range(3)]
            cc = np.outer(wa, wb)
            blk = np.empty((len(mc.CART[la]), len(mc.CART[lb])))
            for i, ea in enumerate(mc.CART[la]):
                for j, eb in enumerate(mc.CART[lb]):
                    e = list(eb)
                    e[k] += 1
                    blk[i, j] = np.sum(cc * tx[0][ea[0]][e[0]] * tx[1][ea[1]][e[1]] * tx[2][ea[2]][e[2]])
            out[i0 : i0 + ta.shape[0], j0 : j0 + tb.shape[0]] = ta @ blk @ tb.T
            j0 += tb.shape[0]
        i0 += ta.shape[0]
    return out


@pytest.mark.parametrize("structure", [CAFFEINE, CU], ids=["caffeine", "cu-doublet"])
@pytest.mark.parametrize("method,conventions", SETUPS)
def test_wavefunction_outputs(method, conventions, structure):
    output = dict(OUT)
    if not program_hamiltonian(method, conventions):
        del output["hamiltonian"]
    out = axtb.compute(inp(method, conventions, structure, output))
    ProgramOutput.model_validate(out.model_dump())  # qcdata's checks: electrons, shells, shapes
    r = out.results
    wf = r.wavefunction
    orb = wf.orbitals
    n = len(structure.symbols)
    restricted = method == "gfn2" or structure.multiplicity == 1
    assert orb.spin.value == ("restricted" if restricted else "unrestricted")
    assert (orb.energies_b is None) == restricted and (orb.coefficients_b is None) == restricted
    assert (
        len(orb.core_electrons) == n
        and (orb.core_electrons[structure.symbols.index("C")] if "C" in structure.symbols else 2) == 2
    )
    held = orb.occupations_a.sum() + orb.occupations_b.sum()
    assert abs(held - (sum(structure.atomic_numbers) - sum(orb.core_electrons) - structure.charge)) < 1e-6
    s = wf.overlap
    channels = [(orb.coefficients_a, orb.energies_a, wf.hamiltonian_a)]
    if not restricted:
        channels.append((orb.coefficients_b, orb.energies_b, wf.hamiltonian_b))
    for c, e, h in channels:
        assert np.abs(c.T @ s @ c - np.eye(c.shape[1])).max() < 1e-12
        if h is not None:
            assert np.abs(h @ c - s @ c * e).max() < 1e-10
    assert [sh.center for sh in wf.basis.shells] == r.shell_atoms
    assert [sh.angular_momentum for sh in wf.basis.shells] == r.shell_angular_momenta
    # the overlap of the basis set, computed apart (the program's screening: ~3e-9)
    m, b = independent(wf, structure)
    s_ind = b.overlap(b)
    assert np.abs(s_ind - m.overlap).max() < 1e-9
    # atomic dipoles from the density and dipole integrals: mu_A = -sum_{mu on A} (P D^A)_mu,mu
    ca = m.orbitals.coefficients_a
    cb = m.orbitals.coefficients_b if m.orbitals.coefficients_b is not None else ca
    dens = (ca * m.orbitals.occupations_a) @ ca.T + (cb * m.orbitals.occupations_b) @ cb.T
    atom_of = np.repeat([sh.center for sh in m.basis.shells], [2 * sh.angular_momentum + 1 for sh in m.basis.shells])
    mu = np.zeros((n, 3))
    for k in range(3):
        diag = np.einsum("mn,nm->m", dens, dipole_integrals(b, k))
        np.add.at(mu[:, k], atom_of, -diag)
    # (measured: 3e-11 at most over these setups; the overlap 1.5e-10, the program's screening)
    assert np.abs(mu - r.atomic_dipoles).max() < 1e-9
    # populations: Mulliken shell populations of the orbitals' density
    pops = np.zeros(len(r.shell_populations))
    shell_of = np.repeat(np.arange(len(m.basis.shells)), [2 * sh.angular_momentum + 1 for sh in m.basis.shells])
    np.add.at(pops, shell_of, np.einsum("mn,nm->m", dens, m.overlap))
    assert np.abs(pops - r.shell_populations).max() < 1e-9
    # the atomic multipoles add up to the molecular moments (exactly: one set of numbers)
    pos = np.asarray(structure.geometry, dtype=float)
    q, d = r.atomic_charges, r.atomic_dipoles
    assert np.abs((q[:, None] * pos).sum(0) + d.sum(0) - r.dipole_moment).max() < 1e-12
    quad = sum(
        r.atomic_quadrupoles[a]
        + 1.5 * (q[a] * np.outer(pos[a], pos[a]) + np.outer(pos[a], d[a]) + np.outer(d[a], pos[a]))
        - 0.5 * (q[a] * pos[a] @ pos[a] + 2 * pos[a] @ d[a]) * np.eye(3)
        for a in range(n)
    )
    assert np.abs(quad - r.quadrupole_moment).max() < 1e-11
    assert np.abs(np.trace(r.atomic_quadrupoles, axis1=1, axis2=2)).max() < 1e-12


@pytest.mark.parametrize("method,conventions", SETUPS)
def test_outputs_change_no_other_number(method, conventions):
    plain = axtb.compute(inp(method, conventions, CAFFEINE, calctype="gradient")).results
    output = {**OUT, "bond_orders": True}
    if not program_hamiltonian(method, conventions):
        del output["hamiltonian"]
    full = axtb.compute(inp(method, conventions, CAFFEINE, output, calctype="gradient")).results
    for name in ("energy", "gradient", "atomic_charges", "dipole_moment", "atomic_dipoles", "shell_populations"):
        assert np.array_equal(getattr(plain, name), getattr(full, name)), name
    for name in ("energies_a", "occupations_a", "occupations_b"):
        assert np.array_equal(getattr(plain.wavefunction.orbitals, name), getattr(full.wavefunction.orbitals, name))


def test_levels_always_coefficients_on_request():
    r = axtb.compute(inp("gfn2", "axtb", CAFFEINE)).results
    assert r.wavefunction.orbitals.coefficients_a is None and r.wavefunction.basis is None
    assert r.wavefunction.orbitals.electronic_temperature == 300.0
    assert r.wavefunction.orbitals.fermi_level_a is not None
    r = axtb.compute(inp("gxtb-v2.0.1", "axtb", CAFFEINE)).results
    assert r.wavefunction.orbitals.electronic_temperature == 0.0 and r.wavefunction.orbitals.fermi_level_a is None
    # an overlap alone ships the basis set and no coefficients (decision 25a)
    r = axtb.compute(inp("gfn2", "axtb", CAFFEINE, {"overlap": True})).results
    assert r.wavefunction.basis is not None and r.wavefunction.overlap is not None
    assert r.wavefunction.orbitals.coefficients_a is None
    # the hamiltonian brings the overlap and the orbitals
    r = axtb.compute(inp("gfn2", "axtb", CAFFEINE, {"hamiltonian": True})).results
    assert r.wavefunction.overlap is not None and r.wavefunction.orbitals.coefficients_a is not None


def test_composite_wavefunction_is_the_gxtb_parts():
    comp = axtb.compute(inp("composite", "axtb", CAFFEINE, {"orbitals": True}))
    base = axtb.compute(inp("gxtb-v2.0.1", "axtb", CAFFEINE, {"orbitals": True}))
    assert comp.record["wavefunction"] == {"part": "gxtb"}
    assert np.array_equal(comp.results.wavefunction.orbitals.energies_a, base.results.wavefunction.orbitals.energies_a)


def test_session_outputs_per_request():
    with axtb.Session(inp("gxtb-v2.0.1", "axtb", CU)) as s:
        f = s.compute(outputs=["hamiltonian"])
        assert f.results.wavefunction.hamiltonian_b is not None
        g = s.compute()
        assert g.results.wavefunction.basis is None
        with pytest.raises(ProgramError) as e:
            s.compute(outputs=["core_hamiltonian"])
        assert (e.value.error.kind.value, e.value.error.code) == ("unsupported", "not_applicable")
    with axtb.Session(inp("gxtb-v2.0.1", "gxtb-2.0.1", CAFFEINE)) as s:
        with pytest.raises(ProgramError) as e:
            s.compute(outputs=["hamiltonian"])
        assert e.value.error.code == "not_applicable"
    gfnff = ProgramInput(program="axtb", calctype="energy", structure=CAFFEINE, model={"method": "gfnff"})
    with axtb.Session(gfnff) as s:
        assert s.compute().results.wavefunction is None
        with pytest.raises(ProgramError) as e:
            s.compute(outputs=["orbitals"])
        assert (e.value.error.kind.value, e.value.error.code) == ("unsupported", "not_applicable")


def test_spin_without_electrons_at_finite_temperature():
    h = Structure(symbols=["H"], geometry=[[0, 0, 0]], multiplicity=2)
    out = axtb.compute(ProgramInput(program="axtb", calctype="energy", structure=h, model={"method": "gfn2"}))
    assert out.results.wavefunction is None
    assert [w.code for w in out.warnings] == ["orbitals_omitted"]
