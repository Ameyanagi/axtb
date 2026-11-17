// SPDX-License-Identifier: LGPL-3.0-or-later
//! Requesting an output changes no other returned number: a single point with the gradient,
//! with and without each optional output (orbitals as a Molden file, bond orders), gives the
//! same bits in the energy, the terms, the gradient, the charges, the multipoles, the orbital
//! energies and occupations and the iterations, for every SCF method in each of its conventions
//! (and the g-xTB solvation composite), closed and open shells, small systems (below the 200
//! rows of the eigensolver's two-sided reduction) and large ones (above the ~313 rows of faer's
//! parallel tridiagonalization, on four threads).

use axtb_driver::{run, Request, Response};
use std::path::Path;

/// Atomic numbers and positions (bohr) of a structure of the repository's `tests/structures`.
fn structure(name: &str) -> Request {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/structures").join(name);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let (numbers, positions) = axtb_driver::read_xyz(&text, axtb_driver::Method::Gfn2, None).unwrap();
    Request::new(numbers, positions)
}

/// The bits of every number a response carries besides the requested outputs.
fn numbers(r: &Response) -> Vec<(&'static str, Vec<u64>)> {
    let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
    let flat3 = |v: &[[f64; 3]]| bits(&v.iter().flatten().copied().collect::<Vec<_>>());
    let flat6 = |v: &[[f64; 6]]| bits(&v.iter().flatten().copied().collect::<Vec<_>>());
    let s = r.scf.as_ref().expect("an SCF method");
    let terms: Vec<f64> = r.terms.as_ref().unwrap().0.iter().map(|t| t.1).collect();
    vec![
        ("energy", bits(&[r.energy.unwrap()])),
        ("terms", bits(&terms)),
        ("gradient", flat3(r.gradient.as_ref().unwrap())),
        ("charges", bits(r.charges.as_ref().unwrap())),
        ("dipole", bits(&r.dipole.unwrap())),
        ("shell_charges", bits(&s.shell_charges)),
        ("atomic_dipoles", flat3(&s.atomic_dipoles)),
        ("atomic_quadrupoles", flat6(&s.atomic_quadrupoles)),
        ("quadrupole", bits(&s.quadrupole)),
        ("orbital_energies_alpha", bits(&s.orbital_energies_alpha)),
        ("orbital_energies_beta", bits(&s.orbital_energies_beta)),
        ("occupations_alpha", bits(&s.occupations_alpha)),
        ("occupations_beta", bits(&s.occupations_beta)),
        ("homo_lumo_gap", bits(&[s.homo_lumo_gap.unwrap_or(f64::NAN)])),
        ("iterations", vec![s.iterations as u64]),
    ]
}

/// The optional outputs of a single point, each as a change of the request.
#[allow(clippy::type_complexity)]
fn outputs() -> Vec<(&'static str, fn(&mut Request))> {
    vec![("orbitals", |r| r.molden = true), ("bond_orders", |r| r.bond_orders = true)]
}

/// Every number of `req` with the gradient alone equals those with each optional output.
fn check(label: &str, req: &Request) {
    let plain = Request { gradient: true, ..req.clone() };
    let base = run(&plain).unwrap_or_else(|e| panic!("{label}: {e}"));
    let want = numbers(&base);
    for (name, set) in outputs() {
        let mut r = plain.clone();
        set(&mut r);
        let got = run(&r).unwrap_or_else(|e| panic!("{label} with {name}: {e}"));
        for ((q, a), (_, b)) in want.iter().zip(&numbers(&got)) {
            assert!(a == b, "{label}: {q} changes when {name} is requested");
        }
    }
}

fn with(req: &Request, method: &str, conventions: &str) -> Request {
    Request { method: Some(method.into()), conventions: Some(conventions.into()), ..req.clone() }
}

/// GFN2-xTB in each of its conventions: water and caffeine (closed shells, 6 and 66 basis
/// functions), O2 and the Cu(NH3)4 2+ doublet (open shells), (H2O)100 (600) on four threads.
#[test]
fn gfn2_outputs_change_no_number() {
    let water = structure("water.xyz");
    let caffeine = structure("caffeine.xyz");
    let o2 = Request {
        uhf: Some(2),
        ..structure("../../crates/gxtb/tests/reference_v2/reference/dioxygen_triplet/mol.xyz")
    };
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let large = Request { threads: Some(4), ..structure("large/water100.xyz") };
    for conv in ["axtb", "xtb", "xtb-6.7.1", "tblite"] {
        for (name, mol) in [("water", &water), ("caffeine", &caffeine), ("o2", &o2), ("cu", &cu), ("(H2O)100", &large)]
        {
            check(&format!("gfn2 {conv} {name}"), &with(mol, "gfn2", conv));
        }
    }
    // a finite temperature: fractional occupations above the occupied orbitals' cutoff
    check("gfn2 caffeine 3000 K", &Request { etemp: Some(3000.0), ..with(&caffeine, "gfn2", "axtb") });
}

/// g-xTB 2.0.1 and v1.1 in both conventions: water and caffeine (closed shells, 166 basis
/// functions), the Cu(NH3)4 2+ doublet (two spin channels), taxol (762) on four threads; the
/// solvation composite.
#[test]
fn gxtb_outputs_change_no_number() {
    let caffeine = structure("caffeine.xyz");
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let large = Request { threads: Some(4), ..structure("taxol.xyz") };
    for method in ["gxtb-v2.0.1", "gxtb-v1.1"] {
        for conv in ["axtb", "gxtb"] {
            for (name, mol) in [("water", &structure("water.xyz")), ("caffeine", &caffeine), ("cu", &cu)] {
                check(&format!("{method} {conv} {name}"), &with(mol, method, conv));
            }
        }
        check(&format!("{method} axtb taxol"), &with(&large, method, "axtb"));
        let etemp = Request { etemp: Some(3000.0), ..with(&caffeine, method, "axtb") };
        check(&format!("{method} caffeine 3000 K"), &etemp);
    }
    check("gxtb-v2.0.1 gxtb taxol", &with(&large, "gxtb-v2.0.1", "gxtb"));
    let composite = Request {
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        ..with(&caffeine, "gxtb-v2.0.1", "axtb")
    };
    check("gxtb-v2.0.1 composite caffeine", &composite);
}
