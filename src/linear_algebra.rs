//! Tropical linear algebra: tropical eigenvalues, tropical determinant.
//!
//! Tropical matrices operate in the max-plus semiring.
//! Tropical determinant = max over all permutations of (sum of selected entries).

use crate::semiring::Tropical;
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::fmt;

/// A tropical matrix: entries in the max-plus semiring.
#[derive(Debug, Clone)]
pub struct TropicalMatrix {
    /// Underlying matrix of tropical entries.
    pub data: DMatrix<Tropical>,
}

impl TropicalMatrix {
    /// Create a tropical matrix from a Vec of rows.
    pub fn from_rows(rows: &[Vec<f64>]) -> Self {
        let ncols = rows.first().map(|r| r.len()).unwrap_or(0);
        let flat: Vec<Tropical> = rows
            .iter()
            .flat_map(|r| r.iter().map(|&v| Tropical::new(v)))
            .collect();
        let nrows = rows.len();
        TropicalMatrix {
            data: DMatrix::from_row_slice(nrows, ncols, &flat),
        }
    }

    /// Create from a nalgebra DMatrix of f64.
    pub fn from_dmatrix(m: &DMatrix<f64>) -> Self {
        TropicalMatrix {
            data: m.map(|v| Tropical::new(v)),
        }
    }

    /// Create an n×n identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut data = DMatrix::repeat(n, n, Tropical::ZERO);
        for i in 0..n {
            data[(i, i)] = Tropical::ONE;
        }
        TropicalMatrix { data }
    }

    /// Create a matrix filled with tropical zero.
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        TropicalMatrix {
            data: DMatrix::repeat(nrows, ncols, Tropical::ZERO),
        }
    }

    /// Number of rows.
    pub fn nrows(&self) -> usize {
        self.data.nrows()
    }

    /// Number of columns.
    pub fn ncols(&self) -> usize {
        self.data.ncols()
    }

    /// Get entry at (i, j).
    pub fn get(&self, i: usize, j: usize) -> Tropical {
        self.data[(i, j)]
    }

    /// Set entry at (i, j).
    pub fn set(&mut self, i: usize, j: usize, val: Tropical) {
        self.data[(i, j)] = val;
    }

    /// Tropical matrix addition (entrywise max).
    pub fn tropical_add(&self, other: &TropicalMatrix) -> TropicalMatrix {
        let mut result = self.data.clone();
        for i in 0..self.nrows() {
            for j in 0..self.ncols() {
                result[(i, j)] = result[(i, j)] + other.data[(i, j)];
            }
        }
        TropicalMatrix { data: result }
    }

    /// Tropical matrix multiplication.
    /// (AB)_{ij} = max_k (A_{ik} + B_{kj})
    pub fn tropical_mul(&self, other: &TropicalMatrix) -> TropicalMatrix {
        let n = self.nrows();
        let m = other.ncols();
        let k = self.ncols();
        let mut result = DMatrix::repeat(n, m, Tropical::ZERO);
        for i in 0..n {
            for j in 0..m {
                let mut best = Tropical::ZERO;
                for l in 0..k {
                    let val = self.data[(i, l)] * other.data[(l, j)];
                    best = best + val;
                }
                result[(i, j)] = best;
            }
        }
        TropicalMatrix { data: result }
    }

    /// Tropical scalar multiplication (add scalar to all entries).
    pub fn scalar_mul(&self, scalar: Tropical) -> TropicalMatrix {
        TropicalMatrix {
            data: self.data.map(|v| v * scalar),
        }
    }

    /// Tropical matrix power: A^k (tropical).
    pub fn tropical_pow(&self, k: u32) -> TropicalMatrix {
        let n = self.nrows();
        if k == 0 {
            return TropicalMatrix::identity(n);
        }
        let mut result = self.clone();
        for _ in 1..k {
            result = result.tropical_mul(self);
        }
        result
    }

    /// Compute the tropical determinant.
    /// det_trop(A) = max over all permutations σ of (sum_i A_{i,σ(i)}).
    pub fn tropical_determinant(&self) -> Tropical {
        let n = self.nrows();
        assert_eq!(n, self.ncols(), "Matrix must be square");
        if n == 0 {
            return Tropical::ONE;
        }
        if n == 1 {
            return self.data[(0, 0)];
        }

        // Generate all permutations and compute max
        let mut perm: Vec<usize> = (0..n).collect();
        let mut best = Tropical::ZERO;
        loop {
            let mut sum = Tropical::ONE;
            for i in 0..n {
                sum = sum * self.data[(i, perm[i])];
            }
            best = best + sum;
            if !next_permutation(&mut perm) {
                break;
            }
        }
        best
    }

    /// Compute tropical eigenvalue (max cycle mean).
    /// λ = max_i max_k (A^k)_{ii} / k
    /// Uses Karp's algorithm.
    pub fn tropical_eigenvalue(&self) -> Option<Tropical> {
        let n = self.nrows();
        assert_eq!(n, self.ncols(), "Matrix must be square");
        if n == 0 {
            return None;
        }

        // Karp's algorithm for max cycle mean
        // v_k[i] = max over all paths of length k ending at i
        let mut v: Vec<Tropical> = vec![Tropical::ONE; n]; // v_0 = 0 for all
        let mut v_history: Vec<Vec<Tropical>> = vec![v.clone()];

        // Compute v_1, ..., v_n
        for _ in 0..n {
            let mut new_v = vec![Tropical::ZERO; n];
            for j in 0..n {
                for i in 0..n {
                    // new_v[j] = max_i (v[i] + A[i][j])
                    let candidate = v[i] * self.data[(i, j)];
                    new_v[j] = new_v[j] + candidate;
                }
            }
            v = new_v;
            v_history.push(v.clone());
        }

        // Compute max cycle mean
        let mut best_mean = Tropical::ZERO;
        for j in 0..n {
            for k in 0..n {
                // (v_n[j] - v_k[j]) / (n - k)
                let numerator = v_history[n][j] - v_history[k][j];
                let denominator = (n - k) as f64;
                if denominator > 0.0 {
                    let mean = Tropical::new(numerator.value() / denominator);
                    best_mean = best_mean + mean;
                }
            }
        }

        if best_mean.is_zero() {
            None
        } else {
            Some(best_mean)
        }
    }

    /// Tropical eigenvectors corresponding to the eigenvalue.
    /// Solve (A ⊗ v) = λ ⊗ v in the tropical sense.
    pub fn tropical_eigenvectors(&self, eigenvalue: Tropical) -> Option<DVector<Tropical>> {
        let n = self.nrows();
        // v_j such that max_i(A_{ji} + v_i) = λ + v_j
        // => A_{ji} + v_i ≤ λ + v_j for all i,j
        // => v_j ≥ A_{ji} + v_i - λ for all i,j
        // This defines a system of tropical inequalities

        // One eigenvector: v_j = max_i(A_{ji} - λ + ... )
        // Use Kleene star approach
        let mut v = DVector::from_element(n, Tropical::ONE);

        for _ in 0..n {
            let mut new_v = DVector::from_element(n, Tropical::ZERO);
            for j in 0..n {
                let mut best = Tropical::ZERO;
                for i in 0..n {
                    // A_{ji} + v_i - λ
                    let candidate = self.data[(j, i)] * v[i] - eigenvalue;
                    best = best + candidate;
                }
                new_v[j] = best;
            }
            // Check convergence
            let mut converged = true;
            for j in 0..n {
                if new_v[j] != v[j] {
                    converged = false;
                    break;
                }
            }
            v = new_v;
            if converged {
                break;
            }
        }

        Some(v)
    }

    /// Compute the Kleene star A* = I ⊕ A ⊕ A² ⊕ ... ⊕ Aⁿ
    pub fn kleene_star(&self) -> TropicalMatrix {
        let n = self.nrows();
        let mut result = TropicalMatrix::identity(n);
        let mut power = self.clone();
        for _ in 0..n {
            result = result.tropical_add(&power);
            power = self.tropical_mul(&power);
        }
        result
    }

    /// Tropical trace: max of diagonal entries.
    pub fn tropical_trace(&self) -> Tropical {
        let n = self.nrows().min(self.ncols());
        (0..n).map(|i| self.data[(i, i)]).fold(Tropical::ZERO, |a, b| a + b)
    }

    /// Transpose.
    pub fn transpose(&self) -> TropicalMatrix {
        TropicalMatrix {
            data: self.data.transpose(),
        }
    }

    /// Check if the matrix is square.
    pub fn is_square(&self) -> bool {
        self.nrows() == self.ncols()
    }

    /// Tropical matrix-vector multiplication.
    pub fn mul_vector(&self, v: &DVector<Tropical>) -> DVector<Tropical> {
        let n = self.nrows();
        let mut result = DVector::from_element(n, Tropical::ZERO);
        for i in 0..n {
            let mut best = Tropical::ZERO;
            for j in 0..self.ncols() {
                best = best + (self.data[(i, j)] * v[j]);
            }
            result[i] = best;
        }
        result
    }
}

/// Generate the next permutation in lexicographic order.
fn next_permutation(perm: &mut Vec<usize>) -> bool {
    let n = perm.len();
    if n < 2 {
        return false;
    }
    // Find the largest k such that perm[k] < perm[k+1]
    let mut k = n - 2;
    while k > 0 && perm[k] >= perm[k + 1] {
        k -= 1;
    }
    if perm[k] >= perm[k + 1] {
        return false;
    }
    // Find the largest l > k such that perm[k] < perm[l]
    let mut l = n - 1;
    while perm[k] >= perm[l] {
        l -= 1;
    }
    perm.swap(k, l);
    perm[k + 1..].reverse();
    true
}

impl fmt::Display for TropicalMatrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.nrows() {
            for j in 0..self.ncols() {
                write!(f, "{}", self.data[(i, j)])?;
                if j + 1 < self.ncols() {
                    write!(f, " ")?;
                }
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tropical_matrix_creation() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert_eq!(m.nrows(), 2);
        assert_eq!(m.ncols(), 2);
        assert_eq!(m.get(0, 0), Tropical::new(1.0));
    }

    #[test]
    fn test_tropical_matrix_identity() {
        let m = TropicalMatrix::identity(3);
        assert_eq!(m.get(0, 0), Tropical::ONE);
        assert_eq!(m.get(0, 1), Tropical::ZERO);
        assert_eq!(m.get(1, 1), Tropical::ONE);
    }

    #[test]
    fn test_tropical_matrix_add() {
        let a = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let b = TropicalMatrix::from_rows(&[vec![5.0, 1.0], vec![2.0, 6.0]]);
        let c = a.tropical_add(&b);
        assert_eq!(c.get(0, 0), Tropical::new(5.0));
        assert_eq!(c.get(0, 1), Tropical::new(2.0));
        assert_eq!(c.get(1, 0), Tropical::new(3.0));
        assert_eq!(c.get(1, 1), Tropical::new(6.0));
    }

    #[test]
    fn test_tropical_matrix_mul() {
        // [1 2] * [5 1] = [max(1+5,2+2), max(1+1,2+6)] = [max(6,4), max(2,8)] = [6, 8]
        // [3 4]   [2 6]   [max(3+5,4+2), max(3+1,4+6)] = [max(8,6), max(4,10)] = [8, 10]
        let a = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let b = TropicalMatrix::from_rows(&[vec![5.0, 1.0], vec![2.0, 6.0]]);
        let c = a.tropical_mul(&b);
        assert_eq!(c.get(0, 0), Tropical::new(6.0));
        assert_eq!(c.get(0, 1), Tropical::new(8.0));
        assert_eq!(c.get(1, 0), Tropical::new(8.0));
        assert_eq!(c.get(1, 1), Tropical::new(10.0));
    }

    #[test]
    fn test_tropical_determinant_2x2() {
        // det = max(a11+a22, a12+a21)
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let det = m.tropical_determinant();
        // max(1+4, 2+3) = max(5, 5) = 5
        assert_eq!(det, Tropical::new(5.0));
    }

    #[test]
    fn test_tropical_determinant_1x1() {
        let m = TropicalMatrix::from_rows(&[vec![7.0]]);
        assert_eq!(m.tropical_determinant(), Tropical::new(7.0));
    }

    #[test]
    fn test_tropical_determinant_3x3() {
        let m = TropicalMatrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ]);
        let det = m.tropical_determinant();
        // Check: max over 6 permutations
        // σ=(0,1,2): 1+5+9=15
        // σ=(0,2,1): 1+6+8=15
        // σ=(1,0,2): 2+4+9=15
        // σ=(1,2,0): 2+6+7=15
        // σ=(2,0,1): 3+4+8=15
        // σ=(2,1,0): 3+5+7=15
        assert_eq!(det, Tropical::new(15.0));
    }

    #[test]
    fn test_tropical_eigenvalue() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let eig = m.tropical_eigenvalue();
        assert!(eig.is_some());
        // Eigenvalue should be ≥ max diagonal = max(1,4) = 4
        assert!(eig.unwrap() >= Tropical::new(4.0));
    }

    #[test]
    fn test_tropical_trace() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert_eq!(m.tropical_trace(), Tropical::new(4.0)); // max(1,4)
    }

    #[test]
    fn test_tropical_transpose() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let t = m.transpose();
        assert_eq!(t.get(0, 0), Tropical::new(1.0));
        assert_eq!(t.get(0, 1), Tropical::new(3.0));
        assert_eq!(t.get(1, 0), Tropical::new(2.0));
    }

    #[test]
    fn test_tropical_pow() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let m2 = m.tropical_pow(2);
        // A² = [max(1+1,0+0), max(1+0,0+1)] = [2, 1]
        //      [max(0+1,1+0), max(0+0,1+1)]   [1, 2]
        assert_eq!(m2.get(0, 0), Tropical::new(2.0));
    }

    #[test]
    fn test_scalar_mul() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0]]);
        let s = m.scalar_mul(Tropical::new(3.0));
        assert_eq!(s.get(0, 0), Tropical::new(4.0));
        assert_eq!(s.get(0, 1), Tropical::new(5.0));
    }

    #[test]
    fn test_mul_vector() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let v = DVector::from_column_slice(&[Tropical::new(5.0), Tropical::new(6.0)]);
        let r = m.mul_vector(&v);
        // [max(1+5, 2+6), max(3+5, 4+6)] = [8, 10]
        assert_eq!(r[0], Tropical::new(8.0));
        assert_eq!(r[1], Tropical::new(10.0));
    }

    #[test]
    fn test_kleene_star() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let star = m.kleene_star();
        // Should include identity + all powers
        assert_eq!(star.nrows(), 2);
    }

    #[test]
    fn test_is_square() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0]]);
        assert!(!m.is_square());
        let m2 = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert!(m2.is_square());
    }

    #[test]
    fn test_eigenvectors() {
        let m = TropicalMatrix::from_rows(&[vec![2.0, 1.0], vec![1.0, 2.0]]);
        if let Some(eig) = m.tropical_eigenvalue() {
            if let Some(_vec) = m.tropical_eigenvectors(eig) {
                // Eigenvector computed successfully
            }
        }
    }

    #[test]
    fn test_display() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0]]);
        let s = format!("{}", m);
        assert!(s.contains("1"));
    }
}
