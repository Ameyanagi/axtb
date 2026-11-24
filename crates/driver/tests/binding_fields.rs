// SPDX-License-Identifier: LGPL-3.0-or-later
//! The fields the Python binding reads from a response and the CLI's JSON leaves out: the
//! structured orbitals (`Task::orbitals`, no Molden text), the core electrons per atom, the shell
//! populations with their atoms, angular momenta and reference occupations, the traceless
//! molecular quadrupole in every conventions, and the thermochemistry's average moment and
//! imaginary frequencies. None of them changes another number or the CLI's JSON.

use axtb_driver::{run, Method, Request};
#[cfg(feature = "gxtb")]
use axtb_driver::{Calculator, Response, Task};
use std::path::Path;

fn structure(name: &str) -> Request {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/structures").join(name);
    let (numbers, positions) =
        axtb_driver::read_xyz(&std::fs::read_to_string(path).unwrap(), Method::Gfn2, None).unwrap();
    Request::new(numbers, positions)
}

#[cfg(feature = "gxtb")]
fn with(req: &Request, method: &str, conv: &str) -> Request {
    Request { method: Some(method.into()), conventions: Some(conv.into()), ..req.clone() }
}

#[cfg(feature = "gxtb")]
fn numbers(r: &Response) -> Vec<u64> {
    let s = r.scf.as_ref().unwrap();
    let mut v = vec![r.energy.unwrap()];
    v.extend(r.gradient.iter().flatten().flatten());
    v.extend(r.charges.iter().flatten());
    v.extend(&s.shell_charges);
    v.extend(&s.orbital_energies_alpha);
    v.extend(&s.occupations_beta);
    v.iter().map(|x| x.to_bits()).collect()
}

/// The structured orbitals come without Molden text and change no number; their coefficients
/// are S-orthonormal with the overlap of the same response (axtb's AO order in both).
#[cfg(feature = "gxtb")]
#[test]
fn structured_orbitals() {
    let caffeine = structure("caffeine.xyz");
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    for (label, req) in [
        ("gfn2 caffeine", with(&caffeine, "gfn2", "axtb")),
        ("gfn2 xtb Cu doublet", with(&cu, "gfn2", "xtb")),
        ("gxtb-v2.0.1 caffeine", with(&caffeine, "gxtb-v2.0.1", "axtb")),
        ("gxtb-v1.1 gxtb Cu doublet", with(&cu, "gxtb-v1.1", "gxtb")),
    ] {
        let plain = Task { gradient: true, ..Task::default() };
        let orbitals = Task { orbitals: true, overlap: true, ..plain.clone() };
        let a = Calculator::new(&req).unwrap().compute(&req.positions, &plain).unwrap();
        let b = Calculator::new(&req).unwrap().compute(&req.positions, &orbitals).unwrap();
        assert!(numbers(&a) == numbers(&b), "{label}: numbers");
        assert!(a.orbitals.is_none() && b.molden.is_none(), "{label}");
        let o = b.orbitals.as_ref().unwrap();
        let s = b.scf.as_ref().unwrap().overlap.as_ref().unwrap();
        let n = o.nao();
        assert_eq!((s.rows, o.sets[0].coefficients.len()), (n, n * n), "{label}");
        for set in &o.sets {
            let c = |k: usize, mu: usize| set.coefficients[k * n + mu];
            for k in 0..n {
                for l in 0..=k {
                    let x: f64 =
                        (0..n).map(|mu| c(k, mu) * (0..n).map(|nu| s.data[mu * n + nu] * c(l, nu)).sum::<f64>()).sum();
                    let want = if k == l { 1.0 } else { 0.0 };
                    assert!((x - want).abs() < 1e-11, "{label}: C^T S C ({k}, {l}) = {x}");
                }
            }
        }
        // the CLI's JSON has none of the binding's fields
        let json = serde_json::to_value(&b).unwrap();
        for key in ["orbitals", "core_electrons", "shell_populations", "quadrupole_traceless", "subtotals"] {
            assert!(json.get(key).is_none(), "{label}: {key} in the JSON");
        }
    }
}

/// Core electrons (Z less the valence reference occupations; Zn 28, Cl 10, C, N, O 2, H 0),
/// shell populations (reference occupation less shell charge, to round-off; their atoms and
/// angular momenta), the traceless quadrupole (the conventions' own where that is traceless).
#[cfg(feature = "gxtb")]
#[test]
fn electrons_shells_and_quadrupole() {
    let zncl4 = Request { charge: -2.0, ..structure("zncl4.xyz") };
    let caffeine = structure("caffeine.xyz");
    for method in ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"] {
        for conv in if method == "gfn2" { ["axtb", "tblite"] } else { ["axtb", "gxtb"] } {
            for (name, req) in [("zncl4", &zncl4), ("caffeine", &caffeine)] {
                let label = format!("{method} {conv} {name}");
                let r = run(&with(req, method, conv)).unwrap();
                let s = r.scf.as_ref().unwrap();
                let core: Vec<usize> = req
                    .numbers
                    .iter()
                    .map(|&z| match z {
                        30 => 28,
                        17 => 10,
                        6..=8 => 2,
                        1 => 0,
                        _ => unreachable!(),
                    })
                    .collect();
                assert_eq!(s.core_electrons, core, "{label}");
                let nsh = s.shell_charges.len();
                assert_eq!(
                    (s.shell_populations.len(), s.shell_atoms.len(), s.shell_angular_momenta.len()),
                    (nsh, nsh, nsh),
                    "{label}"
                );
                assert_eq!(s.reference_shell_occupations.len(), nsh, "{label}");
                assert!(
                    s.shell_atoms.windows(2).all(|w| w[0] <= w[1]) && s.shell_atoms[nsh - 1] == req.numbers.len() - 1
                );
                assert!(s.shell_angular_momenta.iter().all(|&l| l <= 3), "{label}");
                for k in 0..nsh {
                    let q = s.reference_shell_occupations[k] - s.shell_populations[k];
                    assert!(
                        (q - s.shell_charges[k]).abs() < 1e-12,
                        "{label}: shell {k}: {q} vs {}",
                        s.shell_charges[k]
                    );
                }
                let tblite_form = (method, conv) == ("gfn2", "tblite") || (method, conv) == ("gxtb-v2.0.1", "gxtb");
                if tblite_form {
                    assert!(s.quadrupole_traceless != s.quadrupole, "{label}");
                    let trace = s.quadrupole_traceless[0] + s.quadrupole_traceless[2] + s.quadrupole_traceless[5];
                    assert!(trace.abs() < 1e-10, "{label}: trace {trace}");
                } else {
                    let b = |q: &[f64; 6]| q.map(f64::to_bits);
                    assert_eq!(b(&s.quadrupole_traceless), b(&s.quadrupole), "{label}");
                }
            }
        }
    }
}

/// The thermochemistry's average moment of inertia (kg m²; the same within the constants of
/// either conventions) and the imaginary frequencies it inverts (|f| below the cutoff) and
/// leaves out.
#[test]
fn thermo_moment_and_imaginary_frequencies() {
    use axtb::vibrations::{thermochemistry, Conventions, ThermoOptions};
    let water = structure("water.xyz");
    let masses = [15.999, 1.008, 1.008];
    let freqs = [-10.0, -100.0, 1600.0];
    let mut moments = Vec::new();
    for c in [Conventions::Axtb, Conventions::Xtb, Conventions::Xtb671] {
        let t =
            thermochemistry(&water.numbers, &water.positions, &masses, &freqs, c, &ThermoOptions::default()).unwrap();
        assert_eq!((t.n_inverted, t.n_dropped), (1, 1), "{c:?}");
        assert_eq!(
            (t.inverted_frequencies.clone(), t.dropped_frequencies.clone()),
            (vec![-10.0], vec![-100.0]),
            "{c:?}"
        );
        assert!(t.average_moment > 1e-48 && t.average_moment < 1e-46, "{c:?}: {}", t.average_moment);
        moments.push(t.average_moment);
    }
    assert!((moments[0] - moments[1]).abs() < 1e-3 * moments[0], "{moments:?}");
    let r = run(&Request { method: Some("gfn2".into()), hessian: true, ..water }).unwrap();
    let thermo = r.vibrations.unwrap().thermo.unwrap();
    assert!(thermo.average_moment > 0.0 && thermo.inverted_frequencies.is_empty());
}
