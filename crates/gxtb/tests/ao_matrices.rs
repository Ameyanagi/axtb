// SPDX-License-Identifier: LGPL-3.0-or-later
//! The AO matrices of a single point (`Want::overlap`, `Want::hamiltonian`): the returned
//! orbitals are generalized eigenvectors of the returned Hamiltonian with the returned overlap,
//! H C = S C ε to 1e-12 relative, for every method and conventions that offer the Hamiltonian
//! (GFN2-xTB in each of its conventions, g-xTB 2.0.1 in axtb's, v1.1 in both, the solvation
//! composite's g-xTB part), closed and open shells; the Hamiltonian brings the overlap; g-xTB
//! 2.0.1 in the program's conventions refuses the Hamiltonian (it diagonalizes a DIIS-mixed
//! matrix) before anything runs, and offers the overlap.

use axtb::method::{Method, SinglePoint, Want};
use axtb::Error;
use faer::{Mat, MatRef};
use std::path::Path;

/// Atomic numbers and positions (bohr) of a structure of the repository's `tests/structures`.
fn structure(name: &str) -> (Vec<usize>, Vec<[f64; 3]>) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/structures").join(name);
    axtb::io::read_xyz(&std::fs::read_to_string(path).unwrap(), axtb::units::AATOAU).unwrap()
}

const BOTH: Want = Want { orbitals: true, hamiltonian: true, ..Want::ENERGY };

/// max_j |H c_j − ε_j S c_j| / ((|H| + |ε_j| |S|) |c_j|) (Frobenius norms of H and S, 2-norms of
/// the columns): the backward error of the eigenpairs of every channel of `sp`.
fn residual(sp: &SinglePoint, label: &str) -> f64 {
    let e = sp.electronic.as_ref().unwrap();
    let orbitals = e.orbitals.as_ref().unwrap_or_else(|| panic!("{label}: no orbitals"));
    let s = e.overlap.as_ref().unwrap_or_else(|| panic!("{label}: no overlap with the Hamiltonian"));
    let h = e.hamiltonian.as_ref().unwrap_or_else(|| panic!("{label}: no Hamiltonian"));
    let n = orbitals.nao();
    assert_eq!(s.len(), n * n, "{label}: overlap size");
    let sm = MatRef::from_row_major_slice(s, n, n);
    let norm = |m: MatRef<'_, f64>| m.norm_l2();
    let mut worst = 0.0f64;
    for (set_index, set) in orbitals.sets.iter().enumerate() {
        // one Hamiltonian for every set of GFN2-xTB (one spatial set), one per channel of g-xTB
        let hs = &h[if h.len() == 1 { 0 } else { set_index }];
        assert_eq!(hs.len(), n * n, "{label}: Hamiltonian size");
        let hm = MatRef::from_row_major_slice(hs, n, n);
        let c = MatRef::from_column_major_slice(&set.coefficients, n, n);
        let (hc, sc): (Mat<f64>, Mat<f64>) = (hm * c, sm * c);
        let (nh, ns) = (norm(hm), norm(sm));
        for j in 0..n {
            let r: f64 = (0..n).map(|i| (hc[(i, j)] - set.energies[j] * sc[(i, j)]).powi(2)).sum::<f64>().sqrt();
            let scale = (nh + set.energies[j].abs() * ns) * norm(c.col(j).as_mat());
            worst = worst.max(r / scale);
        }
    }
    worst
}

/// Single points of `method` with the orbitals and the Hamiltonian: H C = S C ε.
fn check(label: &str, method: &mut dyn Method, positions: &[[f64; 3]]) {
    let sp = method.singlepoint(positions, &BOTH).unwrap_or_else(|e| panic!("{label}: {e}"));
    let r = residual(&sp, label);
    eprintln!("{label}: H C = S C e to {r:.1e}");
    assert!(r < 1e-12, "{label}: |H C - S C e| {r:e}");
    // the overlap alone comes without the Hamiltonian
    let sp = method.singlepoint(positions, &Want { overlap: true, ..Want::ENERGY }).unwrap();
    let e = sp.electronic.unwrap();
    assert!(e.overlap.is_some() && e.hamiltonian.is_none(), "{label}: the overlap alone");
}

#[test]
fn gfn2_hamiltonian_and_overlap() {
    use axtb::{Calculator, Conventions, Hamiltonian, Options};
    let (caffeine, o2) = (structure("caffeine.xyz"), (vec![8, 8], vec![[0.0, 0.0, -1.14], [0.0, 0.0, 1.14]]));
    for conventions in [Conventions::Axtb, Conventions::Xtb, Conventions::Xtb671, Conventions::Tblite] {
        for ((z, x), uhf, name) in [(&caffeine, None, "caffeine"), (&o2, Some(2), "O2 triplet")] {
            let o = Options { conventions, ..Default::default() };
            let mut calc = Calculator::new(Hamiltonian::Gfn2, z, x, 0.0, uhf, o).unwrap();
            check(&format!("gfn2 {conventions:?} {name}"), &mut calc, x);
        }
    }
    // 600 rows: the reduction by sygst, the virtual orbitals after the SCF on it
    let (z, x) = structure("large/water100.xyz");
    let mut calc = Calculator::new(Hamiltonian::Gfn2, &z, &x, 0.0, None, Options::default()).unwrap();
    check("gfn2 Axtb (H2O)100", &mut calc, &x);
}

#[test]
fn gxtb_hamiltonian_and_overlap() {
    use axtb_gxtb::{Calculator, Conventions, Options, Version};
    let caffeine = structure("caffeine.xyz");
    let cu = structure("cu-nh3-4.xyz");
    let offered =
        [(Version::V2_0_1, Conventions::Axtb), (Version::V1_1, Conventions::Axtb), (Version::V1_1, Conventions::Gxtb)];
    for (version, conventions) in offered {
        for ((z, x), charge, uhf, name) in [(&caffeine, 0.0, None, "caffeine"), (&cu, 2.0, Some(1), "Cu(NH3)4 2+")] {
            let o = Options { conventions, ..Default::default() };
            let mut calc = Calculator::new(version, z, x, charge, uhf, o).unwrap();
            check(&format!("{version:?} {conventions:?} {name}"), &mut calc, x);
        }
    }
    // 762 rows (2.0.1 taxol)
    let (z, x) = structure("taxol.xyz");
    let mut calc = Calculator::new(Version::V2_0_1, &z, &x, 0.0, None, Options::default()).unwrap();
    check("V2_0_1 Axtb taxol", &mut calc, &x);
    // the solvation composite (2.0.1 + ALPB(water) correction): its g-xTB part's matrices and
    // orbitals
    let (z, x) = &caffeine;
    let o = axtb_gxtb::CompositeOptions::default();
    let mut composite = axtb_gxtb::CompositeCalculator::new(z, x, 0.0, None, o).unwrap();
    check("composite caffeine", &mut composite, x);
    // 2.0.1 in the program's conventions: the Hamiltonian is refused before anything runs, the
    // overlap is offered
    let o = Options { conventions: Conventions::Gxtb, ..Default::default() };
    let mut calc = Calculator::new(Version::V2_0_1, z, x, 0.0, None, o).unwrap();
    let want = Want { hamiltonian: true, ..Want::ENERGY };
    assert!(matches!(calc.check(&want, None), Err(Error::NotApplicable(_))));
    assert!(matches!(calc.singlepoint(x, &want), Err(Error::NotApplicable(_))));
    let sp = calc.singlepoint(x, &Want { overlap: true, ..Want::ENERGY }).unwrap();
    assert!(sp.electronic.unwrap().overlap.is_some());
}
