// SPDX-License-Identifier: LGPL-3.0-or-later
//! The results of every method: [`SinglePoint`] and [`Hessian`].

use super::composite::{CompositeReport, PartReport};
use super::HessianMethod;
use crate::orbitals::{OrbitalLevels, Orbitals};

/// Whether an energy term is a piece of the partition of the energy or a subtotal of pieces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TermKind {
    /// A piece: the pieces of a single point sum to its energy (to rounding), in every method.
    Piece,
    /// A sum of other contributions (e.g. the electronic energy).
    Subtotal,
}

/// One energy term (Eh) of a single point, named as the method's reference program reports it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Term {
    pub name: &'static str,
    pub value: f64,
    pub kind: TermKind,
}

impl Term {
    /// A piece of the partition.
    pub fn piece(name: &'static str, value: f64) -> Self {
        Self { name, value, kind: TermKind::Piece }
    }

    /// A subtotal.
    pub fn subtotal(name: &'static str, value: f64) -> Self {
        Self { name, value, kind: TermKind::Subtotal }
    }
}

/// The results of a single point.
#[derive(Clone, Debug)]
pub struct SinglePoint {
    /// Total energy (Eh).
    pub energy: f64,
    /// The method's energy terms in its reference program's order, pieces and subtotals
    /// together ([`TermKind`]).
    pub terms: Vec<Term>,
    /// Eh/bohr, N x 3 (when asked for).
    pub gradient: Option<Vec<[f64; 3]>>,
    /// Atomic charges (e): Mulliken (SCF methods) or the force field's.
    pub charges: Vec<f64>,
    /// Dipole moment about the origin (e bohr).
    pub dipole: [f64; 3],
    /// The electronic structure of an SCF method; `None` without electrons.
    pub electronic: Option<Electronic>,
    /// What only some methods report.
    pub extras: Extras,
    /// Wall time per phase (s).
    pub timings: Vec<(&'static str, f64)>,
}

impl SinglePoint {
    /// The value of the term `name`.
    pub fn term(&self, name: &str) -> Option<f64> {
        self.terms.iter().find(|t| t.name == name).map(|t| t.value)
    }

    /// The sum of the pieces, in order.
    pub fn pieces(&self) -> f64 {
        self.terms.iter().filter(|t| t.kind == TermKind::Piece).map(|t| t.value).sum()
    }

    /// How the SCF went (`None` without an SCF).
    pub fn scf(&self) -> Option<&ScfReport> {
        self.electronic.as_ref().map(|e| &e.scf)
    }
}

/// The electronic structure of an SCF single point.
#[derive(Clone, Debug)]
pub struct Electronic {
    /// Mulliken shell charges, atomic dipoles (N x 3) and quadrupoles (N x 6: xx, xy, yy, xz,
    /// yz, zz, traceless).
    pub shell_charges: Vec<f64>,
    pub atomic_dipoles: Vec<[f64; 3]>,
    pub atomic_quadrupoles: Vec<[f64; 6]>,
    /// The molecular quadrupole about the origin (e bohr²), in the form of
    /// [`MethodInfo::quadrupole`](super::MethodInfo::quadrupole).
    pub quadrupole: [f64; 6],
    /// Orbital energies and occupations per spin channel: one for a closed shell (occupations
    /// up to 2), alpha and beta for an open shell.
    pub levels: Vec<OrbitalLevels>,
    /// HOMO–LUMO gap (Eh, the method's and conventions' rule); `None` without a LUMO.
    pub homo_lumo_gap: Option<f64>,
    /// The Fermi levels (Eh) of the alpha and beta electrons; `None` where the occupations are
    /// the Aufbau filling (at 0 K, and for a spin without electrons). A closed shell has the same
    /// level in both.
    pub fermi_levels: [Option<f64>; 2],
    /// The electronic temperature (K) of the occupations; 0 where they are the Aufbau filling
    /// (0 K, or GFN2-xTB's program conventions up to 0.1 K).
    pub electronic_temperature: f64,
    /// The converged orbitals with their basis ([`Want::orbitals`](super::Want::orbitals)).
    pub orbitals: Option<Orbitals>,
    /// Mayer–Wiberg bond orders, N x N row-major ([`Want::bond_orders`](super::Want::bond_orders)).
    pub bond_orders: Option<Vec<f64>>,
    /// The AO overlap matrix S of the generalized eigenproblem the SCF solved, with the
    /// screening it used (g-xTB: of the SCF basis), n x n row-major in the AO order of
    /// [`Orbitals`] (atom by atom, shell by shell, m = −l..l; p as y, z, x)
    /// ([`Want::overlap`](super::Want::overlap), and with the Hamiltonian).
    pub overlap: Option<Vec<f64>>,
    /// The Hamiltonian H whose generalized eigenvectors with [`Electronic::overlap`] are the
    /// orbitals, H C = S C ε: the operator of one build at the last input, n x n row-major, one
    /// per spin channel of the orbitals (GFN2-xTB: one, open shells included; g-xTB: one for a
    /// closed shell, alpha and beta for an open shell)
    /// ([`Want::hamiltonian`](super::Want::hamiltonian)).
    pub hamiltonian: Option<Vec<Vec<f64>>>,
    /// Mulliken shell populations Σ_{μ∈shell} (P S)_μμ of the total density (electrons), with
    /// the atom and angular momentum of each shell and the reference shell occupations the shell
    /// charges are relative to (charge = reference − population), in the order of the basis'
    /// shells.
    pub shell_populations: Vec<f64>,
    pub shell_atoms: Vec<usize>,
    pub shell_angular_momenta: Vec<usize>,
    pub reference_shell_occupations: Vec<f64>,
    /// The electrons of each atom outside the valence-only orbitals: Z less the reference
    /// occupations of its shells ([`core_electrons`]; Zn: 28 in GFN2-xTB and g-xTB).
    pub core_electrons: Vec<usize>,
    /// The molecular quadrupole about the origin as the traceless (Buckingham) tensor
    /// Θ = Σ_A [θ_A + 3/2 (q r r + r d + d r) − 1/2 (q r² + 2 r·d) 1] (xx, xy, yy, xz, yz, zz; e
    /// bohr²) in every conventions; equal to `quadrupole` where that is the traceless form.
    pub quadrupole_traceless: [f64; 6],
    pub scf: ScfReport,
}

/// The core electrons of each atom of the atomic numbers `numbers`: Z less the reference
/// occupations `reference` of its shells (`shell_atoms`: the atom of each shell), which must be
/// an integer to 1e-6 (a numerical error otherwise).
pub fn core_electrons(numbers: &[usize], shell_atoms: &[usize], reference: &[f64]) -> crate::Result<Vec<usize>> {
    let mut valence = vec![0.0; numbers.len()];
    for (&a, &n) in shell_atoms.iter().zip(reference) {
        valence[a] += n;
    }
    numbers
        .iter()
        .zip(&valence)
        .map(|(&z, &v)| {
            let core = z as f64 - v;
            if (core - core.round()).abs() <= 1e-6 && core.round() >= 0.0 {
                Ok(core.round() as usize)
            } else {
                Err(crate::Error::Numerical(format!("Z = {z} less {v} valence electrons is not a core electron count")))
            }
        })
        .collect()
}

/// How the SCF of a single point went.
#[derive(Clone, Debug, PartialEq)]
pub struct ScfReport {
    /// Iterations of the attempt that converged.
    pub iterations: usize,
    /// Where it started: `previous` (a given or the previous geometry's wavefunction) or the
    /// method's own start (`eeq`, `sad`, `eeq-bc`).
    pub start: &'static str,
    /// The mixer that produced the result (GFN2-xTB: `broyden`, `diis`, `diis-fallback`).
    pub mixer_used: Option<&'static str>,
    /// Every attempt of the method's own plan (GFN2-xTB's mixers); empty without one.
    pub attempts: Vec<Attempt>,
    /// The eigensolver: `faer`, or the run-time LAPACK library.
    pub eigensolver: String,
}

/// One SCF attempt.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Attempt {
    pub mixer: &'static str,
    pub iterations: usize,
    pub converged: bool,
}

/// What only some methods report.
#[derive(Clone, Debug, Default)]
pub struct Extras {
    /// GFN-FF: its topology and the hydrogen and halogen bonds of this geometry.
    pub gfnff: Option<GfnffReport>,
    /// A composite: its definition, components and parts.
    pub composite: Option<Box<CompositeReport>>,
}

/// GFN-FF's report of a single point.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GfnffReport {
    /// Hydrogen bonds (case 1, then cases 2 and 3) and halogen bonds at this geometry.
    pub n_hbonds: [usize; 2],
    pub n_xbonds: usize,
    pub topology: GfnffTopology,
}

/// GFN-FF's topology (set up once; 0-based atom indices).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GfnffTopology {
    pub bonds: Vec<[usize; 2]>,
    pub pi_bond_orders: Vec<f64>,
    pub bond_types: Vec<i32>,
    pub hybridization: Vec<i32>,
    pub pi_atoms: Vec<bool>,
    pub fragments: Vec<usize>,
    pub fragment_charges: Vec<f64>,
    pub topology_charges: Vec<f64>,
    pub n_angles: usize,
    pub n_torsions: usize,
}

/// A nuclear Hessian.
#[derive(Clone, Debug)]
pub struct Hessian {
    /// Energy and gradient at the reference geometry.
    pub energy: f64,
    pub gradient: Vec<[f64; 3]>,
    /// Row-major `3N x 3N` matrix of second derivatives (Eh/bohr²), symmetrized.
    pub hessian: Vec<f64>,
    /// Largest asymmetry |H_ij - H_ji| before symmetrization.
    pub max_asymmetry: f64,
    /// Dipole derivatives ∂μ/∂R_x (e; one row per Cartesian coordinate x = 3·atom + k), for
    /// IR intensities; `None` when not available.
    pub dipole_derivatives: Option<Vec<[f64; 3]>>,
    /// How it was computed.
    pub method: HessianMethod,
    /// Wall time per phase in seconds.
    pub timings: Vec<(&'static str, f64)>,
    /// A composite's parts (with their Hessian methods); empty otherwise.
    pub parts: Vec<PartReport>,
}
