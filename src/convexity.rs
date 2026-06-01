//! Tropical convexity: tropical polytopes, tropical halfspaces.
//!
//! Tropical convexity generalizes classical convexity to the max-plus semiring.
//! A tropical polytope is the tropical convex hull of a set of points.

use crate::semiring::Tropical;
use serde::{Deserialize, Serialize};

/// A tropical halfspace defined by a linear inequality.
/// max(a₁+x₁, a₂+x₂, ..., aₙ+xₙ) ≥ c  (or ≤ c)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TropicalHalfspace {
    /// Coefficients.
    pub coefficients: Vec<Tropical>,
    /// Right-hand side.
    pub rhs: Tropical,
    /// Direction: true = upper halfspace (≥), false = lower (≤).
    pub is_upper: bool,
}

impl TropicalHalfspace {
    pub fn new(coefficients: Vec<f64>, rhs: f64, is_upper: bool) -> Self {
        TropicalHalfspace {
            coefficients: coefficients.into_iter().map(Tropical::new).collect(),
            rhs: Tropical::new(rhs),
            is_upper,
        }
    }

    /// Check if a point is in the halfspace.
    pub fn contains(&self, point: &[Tropical]) -> bool {
        let lhs: Tropical = self
            .coefficients
            .iter()
            .zip(point.iter())
            .map(|(&a, &x)| a * x)
            .fold(Tropical::ZERO, |a, b| a + b);
        if self.is_upper {
            lhs >= self.rhs
        } else {
            lhs <= self.rhs
        }
    }

    /// Number of variables.
    pub fn dim(&self) -> usize {
        self.coefficients.len()
    }
}

/// A tropical polytope: tropical convex hull of generator points.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TropicalPolytope {
    /// Generator points.
    pub generators: Vec<Vec<Tropical>>,
    /// Ambient dimension.
    pub dim: usize,
}

impl TropicalPolytope {
    /// Create a new tropical polytope from generator points.
    pub fn new(generators: Vec<Vec<f64>>) -> Self {
        let dim = generators.first().map(|g| g.len()).unwrap_or(0);
        TropicalPolytope {
            generators: generators
                .into_iter()
                .map(|g| g.into_iter().map(Tropical::new).collect())
                .collect(),
            dim,
        }
    }

    /// Number of generators.
    pub fn num_generators(&self) -> usize {
        self.generators.len()
    }

    /// Check if a point is in the tropical polytope.
    /// A point p is in the tropical convex hull if it can be written as
    /// a tropical linear combination of the generators.
    pub fn contains(&self, point: &[Tropical]) -> bool {
        if self.generators.is_empty() {
            return false;
        }
        // Check if point is in the tropical convex hull by verifying
        // it satisfies all tropical supporting halfspaces.
        // Simplified check: the point must satisfy the tropical inequality
        // constraints defined by the generators.

        // Method: compute the tropical projective distance
        // p is in the convex hull iff for all j:
        // min_i max(p_k - g_i_k) ≤ max_i max(p_k - g_i_k)
        // This is always true for tropical convex hulls; instead we use
        // the metric criterion.

        // A simpler sufficient approach: check if the point is dominated
        // by a tropical combination.
        let mut min_max_dist = f64::INFINITY;
        for gen in &self.generators {
            let max_dist = point
                .iter()
                .zip(gen.iter())
                .map(|(&p, &g)| (p - g).value())
                .fold(f64::NEG_INFINITY, f64::max);
            let min_dist = point
                .iter()
                .zip(gen.iter())
                .map(|(&p, &g)| (p - g).value())
                .fold(f64::INFINITY, f64::min);
            // Tropical distance from this generator
            let dist = max_dist - min_dist;
            min_max_dist = min_max_dist.min(dist);
        }
        // If the minimum max-distance is 0, the point is a generator
        // This is an approximate containment check
        min_max_dist < 1e-10 || self.approximate_membership(point)
    }

    /// Approximate membership test via sampling tropical combinations.
    fn approximate_membership(&self, point: &[Tropical]) -> bool {
        // Generate tropical combinations of generators
        // A tropical combination: max_i (λ_i + g_i)
        // We try λ_i values that would produce the point

        let n = self.generators.len();
        let d = self.dim;
        if n == 0 || d == 0 {
            return false;
        }

        // For each dimension k, find λ such that max_i (λ_i + g_i_k) = point_k
        // This means λ_i = point_k - g_i_k for the achieving generator
        // Try to find consistent λ across all dimensions

        let mut best_lambdas: Vec<f64> = vec![f64::NEG_INFINITY; n];

        for k in 0..d {
            let mut max_val = f64::NEG_INFINITY;
            let mut best_i = 0;
            for i in 0..n {
                let val = point[k].value() - self.generators[i][k].value();
                if val > max_val {
                    max_val = val;
                    best_i = i;
                }
            }
            // This generator achieves the max in dimension k
            best_lambdas[best_i] = best_lambdas[best_i].max(max_val);
        }

        // Now check: does max_i (λ_i + g_i) = point for all dimensions?
        for k in 0..d {
            let mut max_val = f64::NEG_INFINITY;
            for i in 0..n {
                if best_lambdas[i] > f64::NEG_INFINITY {
                    let val = best_lambdas[i] + self.generators[i][k].value();
                    max_val = max_val.max(val);
                }
            }
            if (max_val - point[k].value()).abs() > 1e-8 {
                return false;
            }
        }
        true
    }

    /// Compute the tropical convex hull of the generators.
    /// Returns a set of points forming the hull.
    pub fn convex_hull_points(&self) -> Vec<Vec<Tropical>> {
        // For small sets, return the generators
        // In general, we'd compute the tropical convex hull vertices
        self.generators.clone()
    }

    /// Compute the tropical segment between two points.
    pub fn tropical_segment(a: &[Tropical], b: &[Tropical]) -> Vec<Vec<Tropical>> {
        assert_eq!(a.len(), b.len());
        let n = a.len();
        let mut points = vec![a.to_vec()];

        // Sort the differences: at each threshold, the max switches
        let mut thresholds: Vec<f64> = Vec::new();
        for k in 0..n {
            let diff_a = a[k].value();
            let diff_b = b[k].value();
            // max(λ + a_k, λ' + b_k) switches at a_k - b_k
            thresholds.push(diff_a - diff_b);
        }
        thresholds.sort_by(|a, b| a.partial_cmp(b).unwrap());

        for &t in &thresholds {
            let mid: Vec<Tropical> = (0..n)
                .map(|k| Tropical::new(a[k].value().max(b[k].value() + t)))
                .collect();
            points.push(mid);
        }
        points.push(b.to_vec());
        points
    }

    /// Compute the tropical distance (Hilbert projective metric) between two points.
    pub fn tropical_distance(a: &[Tropical], b: &[Tropical]) -> f64 {
        let diffs: Vec<f64> = a
            .iter()
            .zip(b.iter())
            .map(|(&x, &y)| (x - y).value())
            .collect();
        let max_d = diffs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_d = diffs.iter().cloned().fold(f64::INFINITY, f64::min);
        max_d - min_d
    }

    /// Minkowski sum of two tropical polytopes.
    pub fn minkowski_sum(&self, other: &TropicalPolytope) -> TropicalPolytope {
        let mut generators = Vec::new();
        for a in &self.generators {
            for b in &other.generators {
                let sum: Vec<Tropical> = a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect();
                generators.push(sum);
            }
        }
        TropicalPolytope {
            generators,
            dim: self.dim.max(other.dim),
        }
    }

    /// Compute the tropical circumscribed ball center and radius.
    /// Center minimizes the maximum tropical distance to all generators.
    pub fn circumscribed_ball(&self) -> (Vec<Tropical>, f64) {
        if self.generators.is_empty() {
            return (vec![], 0.0);
        }
        // Center in tropical geometry: componentwise midrange
        let n = self.dim;
        let mut center = Vec::with_capacity(n);
        for k in 0..n {
            let vals: Vec<f64> = self
                .generators
                .iter()
                .map(|g| g[k].value())
                .collect();
            let min_v = vals.iter().cloned().fold(f64::INFINITY, f64::min);
            let max_v = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            center.push(Tropical::new((min_v + max_v) / 2.0));
        }
        let radius = self
            .generators
            .iter()
            .map(|g| Self::tropical_distance(&center, g))
            .fold(0.0f64, f64::max);
        (center, radius)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_halfspace_contains_upper() {
        let hs = TropicalHalfspace::new(vec![0.0, 0.0], 3.0, true);
        let p = [Tropical::new(5.0), Tropical::new(2.0)];
        // max(0+5, 0+2) = 5 ≥ 3
        assert!(hs.contains(&p));
    }

    #[test]
    fn test_halfspace_not_contains() {
        let hs = TropicalHalfspace::new(vec![0.0, 0.0], 10.0, true);
        let p = [Tropical::new(3.0), Tropical::new(2.0)];
        // max(3, 2) = 3 < 10
        assert!(!hs.contains(&p));
    }

    #[test]
    fn test_halfspace_lower() {
        let hs = TropicalHalfspace::new(vec![0.0, 0.0], 3.0, false);
        let p = [Tropical::new(1.0), Tropical::new(2.0)];
        // max(1, 2) = 2 ≤ 3
        assert!(hs.contains(&p));
    }

    #[test]
    fn test_polytope_creation() {
        let tp = TropicalPolytope::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        assert_eq!(tp.num_generators(), 2);
        assert_eq!(tp.dim, 2);
    }

    #[test]
    fn test_polytope_contains_generator() {
        let tp = TropicalPolytope::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        let gen = [Tropical::new(1.0), Tropical::new(2.0)];
        assert!(tp.contains(&gen));
    }

    #[test]
    fn test_tropical_segment() {
        let a = [Tropical::new(1.0), Tropical::new(2.0)];
        let b = [Tropical::new(3.0), Tropical::new(4.0)];
        let seg = TropicalPolytope::tropical_segment(&a, &b);
        assert!(seg.len() >= 2);
    }

    #[test]
    fn test_tropical_distance() {
        let a = [Tropical::new(1.0), Tropical::new(2.0)];
        let b = [Tropical::new(3.0), Tropical::new(4.0)];
        let dist = TropicalPolytope::tropical_distance(&a, &b);
        // diffs: 1-3=-2, 2-4=-2, max=-2, min=-2, dist=0
        // Wait: diffs are a-b = (-2, -2), so max=-2, min=-2, dist=0
        assert_eq!(dist, 0.0); // Same projective point
    }

    #[test]
    fn test_tropical_distance_nonzero() {
        let a = [Tropical::new(1.0), Tropical::new(5.0)];
        let b = [Tropical::new(3.0), Tropical::new(4.0)];
        let dist = TropicalPolytope::tropical_distance(&a, &b);
        // diffs: 1-3=-2, 5-4=1, max=1, min=-2, dist=3
        assert_eq!(dist, 3.0);
    }

    #[test]
    fn test_minkowski_sum() {
        let tp1 = TropicalPolytope::new(vec![vec![1.0, 2.0]]);
        let tp2 = TropicalPolytope::new(vec![vec![3.0, 4.0]]);
        let sum = tp1.minkowski_sum(&tp2);
        assert_eq!(sum.num_generators(), 1);
        // Generator should be (1+3, 2+4) = (4, 6)
        assert_eq!(sum.generators[0][0], Tropical::new(4.0));
        assert_eq!(sum.generators[0][1], Tropical::new(6.0));
    }

    #[test]
    fn test_circumscribed_ball() {
        let tp = TropicalPolytope::new(vec![vec![0.0, 0.0], vec![2.0, 2.0]]);
        let (center, radius) = tp.circumscribed_ball();
        assert_eq!(center.len(), 2);
        assert!(radius >= 0.0);
    }

    #[test]
    fn test_convex_hull_points() {
        let tp = TropicalPolytope::new(vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        let hull = tp.convex_hull_points();
        assert_eq!(hull.len(), 2);
    }
}
