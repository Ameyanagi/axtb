// SPDX-License-Identifier: LGPL-3.0-or-later
//! Dense row-major matrices of both g-xTB versions on faer: square matrix products and
//! conversions between row-major vectors and faer matrices. (The dot products in the two
//! rounding orders the crate uses are `axtb_base::math`'s `dot8` and `dot_seq`.)
//!
//! g-xTB keeps its n x n AO matrices as row-major vectors; the response and Hessian code
//! (2.0.1) and the AO/MO transforms of `axtb::response::MoSpace` work on faer matrices
//! (column-major). Products run on faer's kernels with an explicit `Par` (never faer's global
//! setting).

use faer::{Mat, MatMut, MatRef, Par};

/// `a b` for row-major square matrices (faer on row-major views, no copies), on `par`.
pub fn matmul(n: usize, a: &[f64], b: &[f64], par: Par) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    faer::linalg::matmul::matmul(
        MatMut::from_row_major_slice_mut(&mut out, n, n),
        faer::Accum::Replace,
        MatRef::from_row_major_slice(a, n, n),
        MatRef::from_row_major_slice(b, n, n),
        1.0,
        par,
    );
    out
}

/// Row-major (or symmetric) n x n slice as a faer matrix.
pub fn to_mat(n: usize, v: &[f64]) -> Mat<f64> {
    MatRef::from_row_major_slice(v, n, n).to_owned()
}

/// faer matrix as a row-major vector.
pub fn from_mat(m: &Mat<f64>) -> Vec<f64> {
    let (r, c) = (m.nrows(), m.ncols());
    let mut out = vec![0.0; r * c];
    for j in 0..c {
        let col = m.col_as_slice(j);
        for i in 0..r {
            out[i * c + j] = col[i];
        }
    }
    out
}

/// A symmetric faer matrix as a row-major vector (its columns are its rows: plain copies).
pub fn from_mat_sym(m: &Mat<f64>) -> Vec<f64> {
    let n = m.nrows();
    let mut out = vec![0.0; n * n];
    for j in 0..n {
        out[j * n..(j + 1) * n].copy_from_slice(m.col_as_slice(j));
    }
    out
}

/// A (nearly) symmetric row-major n x n slice as a faer matrix by plain copies (the transpose
/// of [`to_mat`]: equal for a symmetric matrix, round-off apart for a nearly symmetric one).
pub fn to_mat_sym(n: usize, v: &[f64]) -> Mat<f64> {
    let mut m = Mat::<f64>::zeros(n, n);
    for j in 0..n {
        m.col_as_slice_mut(j).copy_from_slice(&v[j * n..(j + 1) * n]);
    }
    m
}
