//! Tropical optimization: solve min-max problems via tropical algebra.
//!
//! Tropical optimization uses the max-plus (or min-plus) semiring to solve
//! optimization problems that can be expressed as tropical polynomial systems.

use crate::convexity::TropicalPolytope;
use crate::linear_algebra::TropicalMatrix;
use crate::semiring::Tropical;
use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// A tropical optimization solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalOptimizer {
    /// Problem dimension.
    pub dim: usize,
}

impl TropicalOptimizer {
    pub fn new(dim: usize) -> Self {
        TropicalOptimizer { dim }
    }

    /// Solve a tropical linear assignment problem.
    /// Given cost matrix C, find assignment x minimizing max_j(C_j + x_j).
    /// In min-plus: min over x of max_j(C_j + x_j)
    pub fn solve_minmax(&self, cost_matrix: &TropicalMatrix) -> DVector<Tropical> {
        let n = cost_matrix.nrows();
        let m = cost_matrix.ncols();
        // For each row i, find the minimum cost column
        let mut result = DVector::from_element(m, Tropical::ONE);
        for j in 0..m {
            let mut best = Tropical::ONE; // 0.0 in tropical
            for i in 0..n {
                let val = cost_matrix.get(i, j);
                best = if val < best { val } else { best };
            }
            result[j] = best;
        }
        result
    }

    /// Solve a tropical eigenvalue optimization.
    /// Find x maximizing (A ⊗ x) / x in tropical sense.
    pub fn eigenvalue_optimization(&self, matrix: &TropicalMatrix) -> Option<Tropical> {
        matrix.tropical_eigenvalue()
    }

    /// Solve the tropical linear programming problem:
    /// minimize c^T ⊗ x subject to A ⊗ x ≤ b
    /// where ≤ is in the tropical sense (componentwise max).
    pub fn tropical_lp(
        &self,
        cost: &[Tropical],
        a: &TropicalMatrix,
        b: &[Tropical],
    ) -> Option<DVector<Tropical>> {
        let n = a.ncols();
        let m = a.nrows();

        // Tropical LP: min max_j(c_j + x_j) s.t. max_j(a_ij + x_j) ≤ b_i
        // This means: a_ij + x_j ≤ b_i for all i,j achieving the max
        // => x_j ≤ b_i - a_ij for all i
        // => x_j ≤ min_i(b_i - a_ij)

        let mut x = DVector::from_element(n, Tropical::ONE);
        for j in 0..n {
            let mut upper_bound = Tropical::new(f64::INFINITY);
            for i in 0..m {
                let bound = b[i] - a.get(i, j);
                if bound < upper_bound {
                    upper_bound = bound;
                }
            }
            x[j] = upper_bound;
        }

        // Verify feasibility
        for i in 0..m {
            let mut row_max = Tropical::ZERO;
            for j in 0..n {
                let val = a.get(i, j) * x[j];
                row_max = row_max + val;
            }
            if row_max > b[i] {
                return None; // Infeasible
            }
        }

        Some(x)
    }

    /// Solve a Chebyshev approximation problem in tropical arithmetic.
    /// Find x minimizing max_i |f_i - (A ⊗ x)_i| tropically.
    pub fn chebyshev_approx(
        &self,
        a: &TropicalMatrix,
        target: &[Tropical],
    ) -> (DVector<Tropical>, Tropical) {
        let n = a.ncols();
        let m = a.nrows();

        // Minimize max_i max(target_i - (Ax)_i, (Ax)_i - target_i)
        // This reduces to finding x such that (Ax)_i is close to target_i

        // Approximate solution: for each i, (Ax)_i ≈ target_i
        // max_j(a_ij + x_j) ≈ target_i => choose x_j ≈ target_i - a_ij

        let mut x = DVector::from_element(n, Tropical::ONE);
        let mut residual = Tropical::ZERO;

        for j in 0..n {
            let mut sum = Tropical::ONE;
            let mut count = 0;
            for i in 0..m {
                let xj_candidate = target[i] - a.get(i, j);
                sum = sum * xj_candidate; // tropical sum = regular addition
                count += 1;
            }
            // Average
            x[j] = Tropical::new(sum.value() / count as f64);
        }

        // Compute residual
        for i in 0..m {
            let mut ax_i = Tropical::ZERO;
            for j in 0..n {
                ax_i = ax_i + (a.get(i, j) * x[j]);
            }
            let diff = (target[i] - ax_i).value().abs();
            residual = residual + Tropical::new(diff); // max of diffs
        }

        (x, residual)
    }

    /// Find the tropical nearest point in a polytope to a given point.
    pub fn tropical_nearest_point(
        &self,
        polytope: &TropicalPolytope,
        point: &[Tropical],
    ) -> Option<Vec<Tropical>> {
        if polytope.generators.is_empty() {
            return None;
        }
        // Find the generator closest in tropical distance
        let mut best_dist = f64::INFINITY;
        let mut best_gen = 0;
        for (i, gen) in polytope.generators.iter().enumerate() {
            let dist = TropicalPolytope::tropical_distance(point, gen);
            if dist < best_dist {
                best_dist = dist;
                best_gen = i;
            }
        }
        Some(polytope.generators[best_gen].clone())
    }

    /// Solve a tropical two-player game.
    /// Player 1 chooses row, Player 2 chooses column.
    /// Payoff matrix in tropical arithmetic.
    pub fn solve_tropical_game(&self, payoff: &TropicalMatrix) -> (Tropical, Tropical) {
        // Player 1 wants to maximize, Player 2 wants to minimize
        // Tropical saddle point value
        let n = payoff.nrows();
        let m = payoff.ncols();

        // For each row, find the min (tropical min = regular min)
        let mut row_mins: Vec<Tropical> = Vec::with_capacity(n);
        for i in 0..n {
            let mut row_min = Tropical::new(f64::INFINITY);
            for j in 0..m {
                let val = payoff.get(i, j);
                if val < row_min {
                    row_min = val;
                }
            }
            row_mins.push(row_min);
        }

        // Player 1 picks the row with the max of row-mins
        let lower = row_mins.into_iter().fold(Tropical::ZERO, |a, b| a + b);

        // For each column, find the max (tropical max)
        let mut col_maxes: Vec<Tropical> = Vec::with_capacity(m);
        for j in 0..m {
            let mut col_max = Tropical::ZERO;
            for i in 0..n {
                let val = payoff.get(i, j);
                col_max = col_max + val;
            }
            col_maxes.push(col_max);
        }

        // Player 2 picks the column with the min of col-maxes
        let upper = col_maxes
            .into_iter()
            .fold(Tropical::new(f64::INFINITY), |a, b| if a < b { a } else { b });

        (lower, upper)
    }

    /// Compute the tropical bottleneck assignment.
    /// Minimize the maximum cost in an assignment.
    pub fn bottleneck_assignment(&self, cost_matrix: &TropicalMatrix) -> Vec<(usize, usize)> {
        let n = cost_matrix.nrows();
        let m = cost_matrix.ncols();
        let mut assignment = Vec::new();

        // Greedy: for each row, pick the column with minimum cost
        let mut used_cols = vec![false; m];
        for i in 0..n {
            let mut best_j = 0;
            let mut best_cost = Tropical::new(f64::INFINITY);
            for j in 0..m {
                if !used_cols[j] {
                    let c = cost_matrix.get(i, j);
                    if c < best_cost {
                        best_cost = c;
                        best_j = j;
                    }
                }
            }
            used_cols[best_j] = true;
            assignment.push((i, best_j));
        }

        assignment
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_optimizer_creation() {
        let opt = TropicalOptimizer::new(3);
        assert_eq!(opt.dim, 3);
    }

    #[test]
    fn test_minmax_solve() {
        let cost = TropicalMatrix::from_rows(&[vec![3.0, 1.0], vec![2.0, 4.0]]);
        let opt = TropicalOptimizer::new(2);
        let result = opt.solve_minmax(&cost);
        assert_eq!(result.nrows(), 2);
    }

    #[test]
    fn test_eigenvalue_opt() {
        let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
        let opt = TropicalOptimizer::new(2);
        let eig = opt.eigenvalue_optimization(&m);
        assert!(eig.is_some());
    }

    #[test]
    fn test_tropical_lp_feasible() {
        let opt = TropicalOptimizer::new(2);
        let cost = vec![Tropical::new(1.0), Tropical::new(2.0)];
        let a = TropicalMatrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let b = vec![Tropical::new(5.0), Tropical::new(3.0)];
        let result = opt.tropical_lp(&cost, &a, &b);
        assert!(result.is_some());
        let x = result.unwrap();
        // x_1 ≤ 5-1=4, x_2 ≤ 3-0=3
        assert!(x[0] <= Tropical::new(4.0));
        assert!(x[1] <= Tropical::new(3.0));
    }

    #[test]
    fn test_chebyshev_approx() {
        let opt = TropicalOptimizer::new(2);
        let a = TropicalMatrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1.0]]);
        let target = [Tropical::new(3.0), Tropical::new(4.0)];
        let (x, residual) = opt.chebyshev_approx(&a, &target);
        assert_eq!(x.nrows(), 2);
        assert!(residual.value() >= 0.0);
    }

    #[test]
    fn test_tropical_nearest_point() {
        let opt = TropicalOptimizer::new(2);
        let tp = TropicalPolytope::new(vec![vec![1.0, 2.0], vec![5.0, 6.0]]);
        let point = [Tropical::new(1.5), Tropical::new(2.5)];
        let nearest = opt.tropical_nearest_point(&tp, &point);
        assert!(nearest.is_some());
        let n = nearest.unwrap();
        // Should be closer to (1,2) than (5,6)
        assert_eq!(n[0], Tropical::new(1.0));
        assert_eq!(n[1], Tropical::new(2.0));
    }

    #[test]
    fn test_tropical_game() {
        let opt = TropicalOptimizer::new(2);
        let payoff = TropicalMatrix::from_rows(&[vec![1.0, 3.0], vec![2.0, 4.0]]);
        let (lower, upper) = opt.solve_tropical_game(&payoff);
        assert!(lower <= upper);
    }

    #[test]
    fn test_bottleneck_assignment() {
        let opt = TropicalOptimizer::new(2);
        let cost = TropicalMatrix::from_rows(&[vec![3.0, 1.0], vec![2.0, 4.0]]);
        let assignment = opt.bottleneck_assignment(&cost);
        assert_eq!(assignment.len(), 2);
        // Each row assigned to a different column
        assert_ne!(assignment[0].1, assignment[1].1);
    }
}
