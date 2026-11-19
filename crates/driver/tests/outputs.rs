// SPDX-License-Identifier: LGPL-3.0-or-later
//! Requesting an output changes no other returned number: a single point with the gradient,
//! with and without each optional output (orbitals as a Molden file, bond orders, the AO overlap,
//! the Hamiltonian), gives the
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
        ("fermi_levels", bits(&s.fermi_levels.map(|x| x.unwrap_or(f64::NAN)))),
        ("electronic_temperature", bits(&[s.electronic_temperature])),
        ("iterations", vec![s.iterations as u64]),
    ]
}

/// The optional outputs of a single point, each as a change of the request.
#[allow(clippy::type_complexity)]
fn outputs() -> Vec<(&'static str, fn(&mut Request))> {
    vec![
        ("orbitals", |r| r.molden = true),
        ("bond_orders", |r| r.bond_orders = true),
        ("overlap", |r| r.overlap = true),
        ("hamiltonian", |r| r.hamiltonian = true),
    ]
}

/// Every number of `req` with the gradient alone equals those with each optional output.
fn check(label: &str, req: &Request) {
    let plain = Request { gradient: true, ..req.clone() };
    let base = run(&plain).unwrap_or_else(|e| panic!("{label}: {e}"));
    let want = numbers(&base);
    for (name, set) in outputs() {
        let mut r = plain.clone();
        set(&mut r);
        let got = match run(&r) {
            Ok(got) => got,
            // g-xTB 2.0.1 in the program's conventions has no Hamiltonian
            Err(e) if name == "hamiltonian" && label.starts_with("gxtb-v2.0.1 gxtb") => {
                assert_eq!(e.kind(), "not_applicable", "{label}: {e}");
                continue;
            }
            Err(e) => panic!("{label} with {name}: {e}"),
        };
        match name {
            "overlap" => assert!(got.scf.as_ref().is_some_and(|s| s.overlap.is_some()), "{label}: no overlap"),
            "hamiltonian" => {
                let s = got.scf.as_ref().unwrap();
                assert!(s.hamiltonian.is_some() && s.overlap.is_some(), "{label}: the Hamiltonian brings the overlap");
            }
            _ => {}
        }
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

/// The Fermi levels per spin and the electronic temperature used: `None` and 0 K for the Aufbau
/// filling (g-xTB's default 0 K, GFN2-xTB at 0 K, and GFN2-xTB's program conventions up to
/// 0.1 K), the same level twice for a closed shell, two levels for an open shell.
#[test]
fn fermi_levels_and_temperature() {
    let water = structure("water.xyz");
    let o2 = Request {
        uhf: Some(2),
        ..structure("../../crates/gxtb/tests/reference_v2/reference/dioxygen_triplet/mol.xyz")
    };
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let at = |req: &Request, etemp: Option<f64>| {
        let s = run(&Request { etemp, ..req.clone() }).unwrap().scf.unwrap();
        (s.fermi_levels, s.electronic_temperature)
    };
    let closed = |(f, t): ([Option<f64>; 2], f64), want: f64, label: &str| {
        assert!(f[0].is_some() && f[0] == f[1], "{label}: {f:?}");
        assert_eq!(t, want, "{label}");
    };
    let aufbau = |r: ([Option<f64>; 2], f64), label: &str| assert_eq!(r, ([None, None], 0.0), "{label}");
    let open = |(f, t): ([Option<f64>; 2], f64), want: f64, label: &str| {
        assert!(f[0].is_some() && f[1].is_some() && f[0] != f[1], "{label}: {f:?}");
        assert_eq!(t, want, "{label}");
    };
    closed(at(&with(&water, "gfn2", "axtb"), None), 300.0, "gfn2 water (300 K)");
    aufbau(at(&with(&water, "gfn2", "axtb"), Some(0.0)), "gfn2 water 0 K");
    aufbau(at(&with(&water, "gfn2", "xtb"), Some(0.05)), "gfn2 xtb water 0.05 K");
    closed(at(&with(&water, "gfn2", "xtb"), Some(0.5)), 0.5, "gfn2 xtb water 0.5 K");
    open(at(&with(&o2, "gfn2", "axtb"), None), 300.0, "gfn2 O2 triplet");
    for method in ["gxtb-v2.0.1", "gxtb-v1.1"] {
        for conv in ["axtb", "gxtb"] {
            aufbau(at(&with(&water, method, conv), None), &format!("{method} {conv} water (0 K)"));
            closed(at(&with(&water, method, conv), Some(3000.0)), 3000.0, &format!("{method} {conv} water 3000 K"));
            open(at(&with(&cu, method, conv), Some(3000.0)), 3000.0, &format!("{method} {conv} Cu doublet 3000 K"));
        }
    }
}

/// AO matrices that do not apply are refused as `not_applicable` (qcdata's code), an input
/// error, before anything runs: the Hamiltonian of g-xTB 2.0.1 in the program's conventions
/// (also as the composite's base; with max_iter 1 its SCF would fail) and both matrices of
/// GFN-FF. Their overlap is offered; the Hamiltonian brings the overlap; an analysis refuses both.
#[test]
fn ao_matrices_not_applicable() {
    let water = structure("water.xyz");
    let refused = |req: Request| {
        let e = run(&req).unwrap_err();
        assert!(e.kind() == "not_applicable" && e.is_input(), "{e:?}");
        assert!(e.to_string().starts_with("not applicable: "), "{e}");
    };
    let program = with(&water, "gxtb-v2.0.1", "gxtb");
    refused(Request { hamiltonian: true, max_iter: Some(1), ..program.clone() });
    refused(Request {
        hamiltonian: true,
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        ..program.clone()
    });
    refused(Request { hamiltonian: true, hessian: true, ..program.clone() });
    refused(Request { overlap: true, ..with(&water, "gfnff", "axtb") });
    refused(Request { hamiltonian: true, ..with(&water, "gfnff", "axtb") });
    let s = run(&Request { overlap: true, ..program }).unwrap().scf.unwrap();
    assert!(s.overlap.is_some() && s.hamiltonian.is_none());
    let s = run(&Request { hamiltonian: true, ..with(&water, "gfn2", "axtb") }).unwrap().scf.unwrap();
    let (o, h) = (s.overlap.unwrap(), s.hamiltonian.unwrap());
    assert_eq!((o.rows, o.cols, h.len(), h[0].rows), (6, 6, 1, 6));
    let analysis = Request {
        overlap: true,
        analyze_hessian: Some(vec![vec![0.0; 9]; 9]),
        ..Request::new(water.numbers.clone(), water.positions.clone())
    };
    assert_eq!(run(&analysis).unwrap_err().kind(), "input");
}

/// Every method's terms partition its energy: the pieces (every term but the subtotals the
/// response names) add up to it to 1e-12 relative. GFN2-xTB in each of its conventions, with
/// ALPB, GBSA, an open shell, 3000 K and without dispersion; g-xTB 2.0.1 and v1.1 in both
/// conventions, closed and open shells; the solvation composite; GFN-FF in each of its
/// conventions. The SCF methods have one subtotal, `electronic`; GFN-FF none.
#[test]
fn terms_partition_the_energy() {
    let caffeine = structure("caffeine.xyz");
    let water = structure("water.xyz");
    let o2 = Request {
        uhf: Some(2),
        ..structure("../../crates/gxtb/tests/reference_v2/reference/dioxygen_triplet/mol.xyz")
    };
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let check = |label: &str, req: Request, subtotals: &[&str]| {
        let r = run(&Request { gradient: true, ..req }).unwrap_or_else(|e| panic!("{label}: {e}"));
        assert_eq!(r.subtotals, subtotals, "{label}");
        let e = r.energy.unwrap();
        let pieces: f64 = r.terms.unwrap().0.iter().filter(|(n, _)| !r.subtotals.contains(n)).map(|(_, v)| v).sum();
        assert!((e - pieces).abs() <= 1e-12 * e.abs(), "{label}: energy - pieces {:e}", e - pieces);
    };
    let elec = &["electronic"][..];
    for conv in ["axtb", "xtb", "xtb-6.7.1", "tblite"] {
        for (name, mol) in [("water", &water), ("caffeine", &caffeine), ("O2 triplet", &o2), ("Cu doublet", &cu)] {
            check(&format!("gfn2 {conv} {name}"), with(mol, "gfn2", conv), elec);
        }
        let solvent = |model: &str| Request {
            alpb: (model == "alpb").then(|| "water".into()),
            gbsa: (model == "gbsa").then(|| "water".into()),
            ..with(&caffeine, "gfn2", conv)
        };
        check(&format!("gfn2 {conv} caffeine ALPB"), solvent("alpb"), elec);
        check(&format!("gfn2 {conv} caffeine GBSA"), solvent("gbsa"), elec);
        check(
            &format!("gfn2 {conv} caffeine 3000 K"),
            Request { etemp: Some(3000.0), ..with(&caffeine, "gfn2", conv) },
            elec,
        );
        check(
            &format!("gfn2 {conv} caffeine no D4"),
            Request { dispersion: Some(false), ..with(&caffeine, "gfn2", conv) },
            elec,
        );
    }
    for method in ["gxtb-v2.0.1", "gxtb-v1.1"] {
        for conv in ["axtb", "gxtb"] {
            for (name, mol) in [("caffeine", &caffeine), ("Cu doublet", &cu)] {
                check(&format!("{method} {conv} {name}"), with(mol, method, conv), elec);
            }
            let composite = Request {
                alpb: Some("water".into()),
                solvation_correction: Some("gfn2".into()),
                ..with(&caffeine, method, conv)
            };
            check(&format!("{method} {conv} composite"), composite, elec);
        }
    }
    for conv in ["axtb", "xtb", "xtb-6.7.1"] {
        check(&format!("gfnff {conv} caffeine"), with(&caffeine, "gfnff", conv), &[]);
    }
}
