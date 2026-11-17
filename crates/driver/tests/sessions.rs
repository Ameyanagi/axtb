// SPDX-License-Identifier: LGPL-3.0-or-later
//! The calls of a session over one calculator: the history of a lineage installed before each
//! call (`Calculator::set_history`) reproduces one calculator per lineage bit for bit, the
//! restart policy given apart from the start (`Request::restart`), the SCF accuracy of later
//! calls (`Calculator::set_scf_accuracy`), the errors' codes and the name and element helpers.

use axtb_driver::{canonical, run, supported_elements, Calculator, Error, Method, Request, Response, Task};
use std::path::Path;

/// Atomic numbers and positions (bohr) of a structure of the repository's `tests/structures`.
fn structure(name: &str) -> Request {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/structures").join(name);
    let (numbers, positions) =
        axtb_driver::read_xyz(&std::fs::read_to_string(path).unwrap(), Method::Gfn2, None).unwrap();
    Request::new(numbers, positions)
}

/// The structure and two displacements of it (golden.py's trajectories; `sign` -1: the other way).
fn trajectory(positions: &[[f64; 3]], sign: f64) -> Vec<Vec<[f64; 3]>> {
    (0..3)
        .map(|k| {
            positions
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    let mut r = *r;
                    for (c, x) in r.iter_mut().enumerate() {
                        *x += sign * 0.01 * k as f64 * (((7 * i + 3 * c) % 5) as f64 - 2.0);
                    }
                    r
                })
                .collect()
        })
        .collect()
}

/// The bits of a response's numbers (and its SCF's start and iterations).
fn bits(r: &Response) -> (Vec<u64>, Option<(&'static str, usize)>) {
    let mut v: Vec<f64> = vec![r.energy.unwrap()];
    v.extend(r.gradient.iter().flatten().flatten());
    v.extend(r.charges.iter().flatten());
    if let Some(s) = &r.scf {
        v.extend(&s.shell_charges);
        v.extend(&s.orbital_energies_alpha);
    }
    (v.iter().map(|x| x.to_bits()).collect(), r.scf.as_ref().map(|s| (s.guess, s.iterations)))
}

/// The trajectory requests: previous-geometry starts in every SCF method and conventions, an
/// open shell, the composite (golden.py's trajectory cases).
fn lineage_requests() -> Vec<(String, Request)> {
    let caffeine = structure("caffeine.xyz");
    let water = structure("water.xyz");
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let mut out = Vec::new();
    let with = |base: &Request, method: &str, conv: &str| Request {
        method: Some(method.into()),
        conventions: Some(conv.into()),
        guess: Some("previous".into()),
        ..base.clone()
    };
    for (method, convs) in
        [("gfn2", ["axtb", "xtb"]), ("gxtb-v2.0.1", ["axtb", "gxtb"]), ("gxtb-v1.1", ["axtb", "gxtb"])]
    {
        for conv in convs {
            out.push((format!("{method} {conv} caffeine"), with(&caffeine, method, conv)));
        }
    }
    out.push(("gxtb-v2.0.1 axtb Cu doublet".into(), with(&cu, "gxtb-v2.0.1", "axtb")));
    for conv in ["axtb", "gxtb"] {
        let composite = Request {
            method: Some("gxtb-v2.0.1".into()),
            conventions: Some(conv.into()),
            alpb: Some("water".into()),
            solvation_correction: Some("gfn2".into()),
            ..water.clone()
        };
        out.push((format!("composite {conv} water"), composite));
    }
    out
}

/// Two lineages interleaved on one calculator, each call after `set_history` with its lineage's
/// wavefunction, equal two calculators run one after the other, bit for bit; one lineage with its
/// own history installed before every call equals the plain calculator.
#[test]
fn interleaved_lineages_equal_sequential_calculators() {
    let task = Task { gradient: true, ..Task::default() };
    for (label, req) in lineage_requests() {
        let (a, b) = (trajectory(&req.positions, 1.0), trajectory(&req.positions, -1.0));
        let sequential = |traj: &[Vec<[f64; 3]>]| {
            let mut calc = Calculator::new(&req).unwrap();
            traj.iter().map(|x| bits(&calc.compute(x, &task).unwrap())).collect::<Vec<_>>()
        };
        let (want_a, want_b) = (sequential(&a), sequential(&b));
        let mut calc = Calculator::new(&req).unwrap();
        let (mut hist_a, mut hist_b) = (None, None);
        for k in 0..3 {
            calc.set_history(hist_a.as_ref()).unwrap();
            let ra = calc.compute(&a[k], &task).unwrap();
            hist_a = calc.wavefunction();
            calc.set_history(hist_b.as_ref()).unwrap();
            let rb = calc.compute(&b[k], &task).unwrap();
            hist_b = calc.wavefunction();
            assert!(bits(&ra) == want_a[k], "{label}: lineage a, step {k}");
            assert!(bits(&rb) == want_b[k], "{label}: lineage b, step {k}");
        }
        // one lineage, its history installed before every call
        let mut calc = Calculator::new(&req).unwrap();
        let mut hist = None;
        for (k, x) in a.iter().enumerate() {
            calc.set_history(hist.as_ref()).unwrap();
            assert!(bits(&calc.compute(x, &task).unwrap()) == want_a[k], "{label}: one lineage, step {k}");
            hist = calc.wavefunction();
        }
    }
}

/// `restart: fresh` starts every SCF as a new calculator would; `restart` and a named start go
/// together (GFN2-xTB from neutral atoms, then from the previous state); `guess` with `previous`
/// or `fresh` beside `restart` and `restart` for GFN-FF are refused.
#[test]
fn restart_apart_from_the_start() {
    let task = Task { gradient: true, ..Task::default() };
    let caffeine = structure("caffeine.xyz");
    let traj = trajectory(&caffeine.positions, 1.0);
    for method in ["gfn2", "gxtb-v2.0.1", "gxtb-v1.1"] {
        let req = Request { method: Some(method.into()), restart: Some("fresh".into()), ..caffeine.clone() };
        let mut calc = Calculator::new(&req).unwrap();
        calc.compute(&traj[0], &task).unwrap();
        let second = calc.compute(&traj[1], &task).unwrap();
        let fresh = Calculator::new(&req).unwrap().compute(&traj[1], &task).unwrap();
        assert!(bits(&second) == bits(&fresh), "{method}: restart fresh");
    }
    let sad = Request {
        method: Some("gfn2".into()),
        guess: Some("sad".into()),
        restart: Some("previous".into()),
        ..caffeine.clone()
    };
    let mut calc = Calculator::new(&sad).unwrap();
    assert_eq!(calc.compute(&traj[0], &task).unwrap().scf.unwrap().guess, "sad");
    assert_eq!(calc.compute(&traj[1], &task).unwrap().scf.unwrap().guess, "previous");
    // the same as guess "sad" alone at the first geometry
    let alone = run(&Request { gradient: true, guess: Some("sad".into()), restart: None, ..sad.clone() }).unwrap();
    let first = Calculator::new(&sad).unwrap().compute(&caffeine.positions, &task).unwrap();
    assert!(bits(&alone) == bits(&first));
    let refused = |req: Request| assert_eq!(Calculator::new(&req).err().map(|e| e.kind()), Some("input"));
    refused(Request { guess: Some("previous".into()), restart: Some("fresh".into()), ..sad.clone() });
    refused(Request { method: Some("gfnff".into()), restart: Some("previous".into()), ..caffeine.clone() });
    refused(Request { restart: Some("sometimes".into()), ..caffeine.clone() });
    let eeq_bc = Request {
        method: Some("gxtb-v2.0.1".into()),
        guess: Some("eeq-bc".into()),
        restart: Some("previous".into()),
        ..caffeine
    };
    let mut calc = Calculator::new(&eeq_bc).unwrap();
    calc.compute(&traj[0], &task).unwrap();
    assert_eq!(calc.compute(&traj[1], &task).unwrap().scf.unwrap().guess, "previous");
}

/// `set_scf_accuracy` changes the SCF thresholds of later calls: g-xTB (no screening by the
/// accuracy) then equals a calculator built with that accuracy, bit for bit; GFN2-xTB keeps the
/// integral screening it was built with (its energy differs from a calculator built with the
/// accuracy by the screening alone, the SCF converges as tightly). GFN-FF has none.
#[test]
fn scf_accuracy_of_later_calls() {
    let task = Task { gradient: true, ..Task::default() };
    let caffeine = structure("caffeine.xyz");
    let x = &caffeine.positions;
    for (method, conv, same_bits) in
        [("gxtb-v2.0.1", "axtb", true), ("gxtb-v1.1", "gxtb", true), ("gfn2", "axtb", false), ("gfn2", "xtb", false)]
    {
        let req = Request { method: Some(method.into()), conventions: Some(conv.into()), ..caffeine.clone() };
        let mut set = Calculator::new(&req).unwrap();
        let loose = set.compute(x, &task).unwrap();
        // the next SCF from the method's own start, as a new calculator's
        set.reset_wavefunction(false);
        set.set_scf_accuracy(1e-3).unwrap();
        let tight = set.compute(x, &task).unwrap();
        let built =
            Calculator::new(&Request { accuracy: Some(1e-3), ..req.clone() }).unwrap().compute(x, &task).unwrap();
        let label = format!("{method} {conv}");
        assert!(tight.scf.as_ref().unwrap().iterations > loose.scf.as_ref().unwrap().iterations, "{label}");
        if same_bits {
            assert!(bits(&tight) == bits(&built), "{label}: set vs built");
        } else {
            // the screening the calculator was built with stays: close to, not the bits of, one
            // built with the accuracy
            let de = (tight.energy.unwrap() - built.energy.unwrap()).abs();
            assert!(de < 1e-6, "{label}: {de:e}");
        }
    }
    let composite = Request {
        method: Some("gxtb-v2.0.1".into()),
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        ..structure("water.xyz")
    };
    let w = composite.positions.clone();
    let mut set = Calculator::new(&composite).unwrap();
    set.set_scf_accuracy(1e-3).unwrap();
    let built = Calculator::new(&Request { accuracy: Some(1e-3), ..composite.clone() }).unwrap();
    let mut built = built;
    assert!(bits(&set.compute(&w, &task).unwrap()) == bits(&built.compute(&w, &task).unwrap()), "composite");
    let mut ff = Calculator::new(&Request { method: Some("gfnff".into()), ..caffeine.clone() }).unwrap();
    assert!(ff.set_scf_accuracy(1e-3).is_err());
    assert!(set.set_scf_accuracy(-1.0).is_err());
}

/// Every error variant's code, the canonical names and the supported elements.
#[test]
fn codes_names_and_elements() {
    let codes = [
        (Error::Input(String::new()), "invalid_input"),
        (Error::UnsupportedElement { z: 0, message: String::new() }, "unsupported_element"),
        (Error::Topology(String::new()), "topology_failed"),
        (Error::ScfNotConverged { iterations: 0, energy_change: 0.0, residual: 0.0 }, "scf_not_converged"),
        (Error::Numerical(String::new()), "numerical_failure"),
        (Error::NotApplicable(String::new()), "not_applicable"),
    ];
    for (e, code) in codes {
        assert_eq!(e.code(), code);
    }
    let water = structure("water.xyz");
    let unsupported = run(&Request { numbers: vec![8, 1, 99], ..water.clone() }).unwrap_err();
    assert_eq!(unsupported.code(), "unsupported_element");
    assert_eq!(canonical("method", "GFN2-xTB", None).unwrap(), "gfn2");
    assert_eq!(canonical("method", "g-xTB 2.0", None).unwrap(), "gxtb-v2.0.1");
    assert_eq!(canonical("conventions", "xtb-6.8.0", Some("gfn2")).unwrap(), "xtb");
    assert_eq!(canonical("conventions", "XTB 6.7.1", Some("gfnff")).unwrap(), "xtb-6.7.1");
    assert_eq!(canonical("conventions", "default", Some("gxtb-v1.1")).unwrap(), "axtb");
    assert_eq!(canonical("solvent", "H2O", None).unwrap(), "water");
    assert_eq!(canonical("solvent", "Methanol", None).unwrap(), "methanol");
    assert!(canonical("solvent", "not a solvent", None).is_err());
    assert!(canonical("conventions", "axtb", None).is_err());
    assert!(canonical("color", "red", None).is_err());
    let range = |a: usize, b: usize| (a..=b).collect::<Vec<_>>();
    assert_eq!(supported_elements(Method::Gfn2), range(1, 86));
    assert_eq!(supported_elements(Method::Gfnff), range(1, 103));
    assert_eq!(supported_elements(Method::GxtbV2_0_1), range(1, 103));
    // v1.1's parameter file: H-Ce, Yb-Ac and U
    let v11 = supported_elements(Method::GxtbV1_1);
    let want: Vec<usize> = (1..=58).chain(70..=89).chain([92]).collect();
    assert_eq!(v11, want);
}

/// The valence electrons of every element of a method's tables are the electrons its
/// calculations' orbitals hold: Z less the core electrons the response reports.
#[test]
fn valence_electrons_are_the_responses() {
    use axtb_driver::valence_electrons;
    let hard = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/scf_hard/random-mbheavy-112.xyz");
    let (numbers, positions) =
        axtb_driver::read_xyz(&std::fs::read_to_string(hard).unwrap(), Method::Gfn2, None).unwrap();
    // La, Gd, Pt, Sb, Se, B: an odd valence count in GFN2-xTB (47), an even one in g-xTB 2.0.1 (54)
    let heavy = Request::new(numbers, positions);
    let zncl4 = Request { charge: -2.0, ..structure("zncl4.xyz") };
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    for (method, molecules) in [
        ("gfn2", vec![structure("caffeine.xyz"), zncl4.clone(), cu.clone(), heavy.clone()]),
        ("gxtb-v2.0.1", vec![structure("caffeine.xyz"), zncl4.clone(), cu.clone(), heavy]),
        ("gxtb-v1.1", vec![structure("caffeine.xyz"), zncl4, cu]),
    ] {
        let m = Method::parse(method).unwrap();
        for req in molecules {
            let req = Request { method: Some(method.into()), accuracy: Some(10.0), ..req };
            let r = run(&req).unwrap();
            let scf = r.scf.as_ref().unwrap();
            for (&z, &core) in req.numbers.iter().zip(&scf.core_electrons) {
                let v = valence_electrons(m, z).unwrap();
                assert!((z as f64 - v - core as f64).abs() < 1e-6, "{method} Z = {z}: {v} valence, {core} core");
            }
        }
    }
    assert_eq!(valence_electrons(Method::Gfnff, 30), Some(30.0));
    assert_eq!(valence_electrons(Method::Gfn2, 87), None);
    assert_eq!(valence_electrons(Method::GxtbV2_0_1, 104), None);
}

/// The densities C n C^T of a response's orbitals, per set (a closed shell's set: occupations up to
/// 2, the total density).
fn densities(r: &Response) -> Vec<Vec<f64>> {
    let o = r.orbitals.as_ref().unwrap();
    let n = o.nao();
    o.sets
        .iter()
        .map(|set| {
            let mut p = vec![0.0; n * n];
            for (k, &occ) in set.occupations.iter().enumerate() {
                let c = &set.coefficients[k * n..(k + 1) * n];
                for mu in 0..n {
                    for nu in 0..n {
                        p[mu * n + nu] += occ * c[mu] * c[nu];
                    }
                }
            }
            p
        })
        .collect()
}

/// A density start (R7): the converged orbitals' density given back starts g-xTB (2.0.1 from its
/// Mulliken state, v1.1 from the densities, the composite's g-xTB part) at the same result in a few
/// iterations; densities of the wrong channels or electrons are refused; GFN2-xTB and GFN-FF start
/// from no density.
#[test]
fn a_density_start() {
    let caffeine = structure("caffeine.xyz");
    let cu = Request { charge: 2.0, uhf: Some(1), ..structure("cu-nh3-4.xyz") };
    let water = structure("water.xyz");
    let composite = Request {
        method: Some("gxtb-v2.0.1".into()),
        alpb: Some("water".into()),
        solvation_correction: Some("gfn2".into()),
        ..water
    };
    let task = Task { gradient: true, orbitals: true, ..Task::default() };
    for (label, req) in [
        ("2.0.1", Request { method: Some("gxtb-v2.0.1".into()), ..caffeine.clone() }),
        (
            "2.0.1 program",
            Request { method: Some("gxtb-v2.0.1".into()), conventions: Some("gxtb".into()), ..caffeine.clone() },
        ),
        ("2.0.1 doublet", Request { method: Some("gxtb-v2.0.1".into()), ..cu.clone() }),
        ("v1.1", Request { method: Some("gxtb-v1.1".into()), ..caffeine.clone() }),
        ("composite", composite),
    ] {
        let first = Calculator::new(&req).unwrap().compute(&req.positions, &task).unwrap();
        let mut calc = Calculator::new(&req).unwrap();
        calc.set_density(densities(&first)).unwrap();
        let again = calc.compute(&req.positions, &task).unwrap();
        let (s0, s1) = (first.scf.as_ref().unwrap(), again.scf.as_ref().unwrap());
        assert!(s1.iterations <= 3, "{label}: {} iterations", s1.iterations);
        assert_eq!(s1.guess, "previous", "{label}");
        assert!((again.energy.unwrap() - first.energy.unwrap()).abs() < 1e-7, "{label}");
        assert!(s0.iterations > s1.iterations, "{label}");
    }
    // the wrong channels and the wrong electrons
    let req = Request { method: Some("gxtb-v2.0.1".into()), ..caffeine.clone() };
    let first = Calculator::new(&req).unwrap().compute(&req.positions, &task).unwrap();
    let p = densities(&first);
    let mut calc = Calculator::new(&req).unwrap();
    let two = vec![p[0].clone(), p[0].clone()];
    assert!(matches!(calc.set_density(two), Err(Error::Input(m)) if m.contains("matrices")));
    let dication = Request { charge: 2.0, ..req.clone() };
    let mut calc = Calculator::new(&dication).unwrap();
    assert!(matches!(calc.set_density(p.clone()), Err(Error::Input(m)) if m.contains("electrons")));
    for method in ["gfn2", "gfnff"] {
        let mut calc = Calculator::new(&Request { method: Some(method.into()), ..caffeine.clone() }).unwrap();
        assert!(calc.set_density(p.clone()).is_err(), "{method}");
    }
}
