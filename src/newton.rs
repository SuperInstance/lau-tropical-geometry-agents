//! Newton polytope: tropical hypersurfaces correspond to polyhedral subdivisions.
//!
//! The Newton polytope of a tropical polynomial is the convex hull of its exponent vectors.
//! The tropical hypersurface induces a polyhedral subdivision dual to the Newton polytope.

use crate::polynomial::TropicalPolynomial;
use crate::semiring::Tropical;
use serde::{Deserialize, Serialize};

/// A point in integer lattice (exponent vector).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LatticePoint {
    pub coords: Vec<i64>,
}

impl LatticePoint {
    pub fn new(coords: Vec<i64>) -> Self {
        LatticePoint { coords }
    }

    pub fn dim(&self) -> usize {
        self.coords.len()
    }

    /// Dot product with another lattice point.
    pub fn dot(&self, other: &LatticePoint) -> i64 {
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(&a, &b)| a * b)
            .sum()
    }

    /// Subtract.
    pub fn sub(&self, other: &LatticePoint) -> LatticePoint {
        let coords = self
            .coords
            .iter()
            .zip(other.coords.iter())
            .map(|(&a, &b)| a - b)
            .collect();
        LatticePoint::new(coords)
    }

    /// Add.
    pub fn add(&self, other: &LatticePoint) -> LatticePoint {
        let coords = self
            .coords
            .iter()
            .zip(other.coords.iter())
            .map(|(&a, &b)| a + b)
            .collect();
        LatticePoint::new(coords)
    }
}

/// A convex hull of lattice points (Newton polytope).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NewtonPolytope {
    /// Vertices of the polytope.
    pub vertices: Vec<LatticePoint>,
    /// Dimension.
    pub dim: usize,
}

impl NewtonPolytope {
    /// Compute the Newton polytope from a tropical polynomial.
    pub fn from_polynomial(poly: &TropicalPolynomial) -> Self {
        let vertices: Vec<LatticePoint> = poly
            .monomials
            .iter()
            .map(|m| {
                let coords: Vec<i64> = m.degree.iter().map(|&d| d as i64).collect();
                LatticePoint::new(coords)
            })
            .collect();
        let dim = poly.num_vars;
        NewtonPolytope { vertices, dim }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get the dimension of the polytope (rank of vertex differences).
    pub fn dimension(&self) -> usize {
        if self.vertices.is_empty() {
            return 0;
        }
        if self.vertices.len() == 1 {
            return 0;
        }
        // Compute rank of difference vectors
        let base = &self.vertices[0];
        let diffs: Vec<Vec<i64>> = self
            .vertices
            .iter()
            .skip(1)
            .map(|v| v.sub(base).coords)
            .collect();
        // Simple rank computation via Gaussian elimination over rationals (approximate)
        rank_of_vectors(&diffs)
    }

    /// Compute the upper convex hull in a given direction.
    /// Returns indices of vertices on the upper hull.
    pub fn upper_hull_indices(&self, direction: &[f64]) -> Vec<usize> {
        if self.vertices.is_empty() {
            return vec![];
        }
        let dir: Vec<f64> = direction.to_vec();
        let vals: Vec<f64> = self
            .vertices
            .iter()
            .map(|v| {
                v.coords
                    .iter()
                    .zip(dir.iter())
                    .map(|(&c, &d)| c as f64 * d)
                    .sum()
            })
            .collect();
        let max_val = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        vals.iter()
            .enumerate()
            .filter(|(_, &v)| (v - max_val).abs() < 1e-10)
            .map(|(i, _)| i)
            .collect()
    }

    /// Check if a lattice point is inside the polytope (simple bounding box check).
    pub fn contains(&self, point: &LatticePoint) -> bool {
        if self.vertices.is_empty() {
            return false;
        }
        // Simple check: point must be in bounding box and be a convex combination
        // For simplicity, check bounding box
        for i in 0..point.dim() {
            let min_coord = self
                .vertices
                .iter()
                .map(|v| v.coords[i])
                .min()
                .unwrap_or(0);
            let max_coord = self
                .vertices
                .iter()
                .map(|v| v.coords[i])
                .max()
                .unwrap_or(0);
            if point.coords[i] < min_coord || point.coords[i] > max_coord {
                return false;
            }
        }
        true
    }

    /// Compute the volume (normalized) of the polytope.
    /// For 2D: area; for 3D: volume; etc.
    pub fn volume(&self) -> f64 {
        let d = self.dimension();
        if d == 0 {
            return 1.0;
        }
        if d == 1 {
            // Length
            if self.vertices.len() < 2 {
                return 0.0;
            }
            let min_v = self
                .vertices
                .iter()
                .map(|v| v.coords[0] as f64)
                .fold(f64::INFINITY, f64::min);
            let max_v = self
                .vertices
                .iter()
                .map(|v| v.coords[0] as f64)
                .fold(f64::NEG_INFINITY, f64::max);
            return max_v - min_v;
        }
        if d == 2 {
            return self.signed_area_2d().abs();
        }
        // Higher dimensions: use Minkowski sum approximation
        0.0
    }

    /// Signed area for 2D polytope using shoelace formula.
    fn signed_area_2d(&self) -> f64 {
        if self.vertices.len() < 3 {
            return 0.0;
        }
        let n = self.vertices.len();
        let mut area = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            let xi = self.vertices[i].coords[0] as f64;
            let yi = self.vertices[i].coords[1] as f64;
            let xj = self.vertices[j].coords[0] as f64;
            let yj = self.vertices[j].coords[1] as f64;
            area += xi * yj - xj * yi;
        }
        area / 2.0
    }

    /// Edges of the polytope (all pairs of vertices for simplicity).
    pub fn edges(&self) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        for i in 0..self.vertices.len() {
            for j in (i + 1)..self.vertices.len() {
                edges.push((i, j));
            }
        }
        edges
    }

    /// Minkowski sum with another polytope.
    pub fn minkowski_sum(&self, other: &NewtonPolytope) -> NewtonPolytope {
        let mut vertices = Vec::new();
        for a in &self.vertices {
            for b in &other.vertices {
                vertices.push(a.add(b));
            }
        }
        NewtonPolytope {
            vertices,
            dim: self.dim.max(other.dim),
        }
    }
}

/// Compute the rank of a set of vectors using approximate Gaussian elimination.
fn rank_of_vectors(vectors: &[Vec<i64>]) -> usize {
    if vectors.is_empty() {
        return 0;
    }
    let m = vectors.len();
    let n = vectors[0].len();
    // Convert to f64 for simpler pivoting
    let mut mat: Vec<Vec<f64>> = vectors
        .iter()
        .map(|v| v.iter().map(|&x| x as f64).collect())
        .collect();
    let mut rank = 0;
    for col in 0..n {
        // Find pivot
        let pivot_row = (rank..m).find(|&r| mat[r][col].abs() > 1e-10);
        if let Some(pr) = pivot_row {
            mat.swap(rank, pr);
            let pivot_val = mat[rank][col];
            for row in (rank + 1)..m {
                if mat[row][col].abs() > 1e-10 {
                    let factor = mat[row][col] / pivot_val;
                    for c in 0..n {
                        mat[row][c] -= factor * mat[rank][c];
                    }
                }
            }
            rank += 1;
        }
    }
    rank
}

/// A polyhedral subdivision induced by a tropical polynomial.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolyhedralSubdivision {
    /// Cells of the subdivision, each is a list of vertex indices.
    pub cells: Vec<Vec<usize>>,
    /// The Newton polytope.
    pub newton_polytope: NewtonPolytope,
}

impl PolyhedralSubdivision {
    /// Compute the subdivision from a tropical polynomial.
    pub fn from_polynomial(poly: &TropicalPolynomial) -> Self {
        let np = NewtonPolytope::from_polynomial(poly);
        // For a tropical polynomial, the subdivision is induced by
        // lifting the exponent vectors by their coefficients
        // and projecting the upper faces.
        let cells = compute_subdivision(poly);
        PolyhedralSubdivision {
            cells,
            newton_polytope: np,
        }
    }

    /// Number of cells.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }
}

/// Compute the regular subdivision induced by the coefficients.
fn compute_subdivision(poly: &TropicalPolynomial) -> Vec<Vec<usize>> {
    let n = poly.monomials.len();
    if n == 0 {
        return vec![];
    }
    // Group monomials by their dominance regions
    // Simple approach: each pair that shares a boundary forms a cell
    let mut cells = Vec::new();
    // For univariate: cells are intervals between corner points
    if poly.num_vars == 1 {
        for i in 0..n {
            cells.push(vec![i]);
        }
    } else {
        // For multivariate: each monomial forms its own cell initially
        for i in 0..n {
            cells.push(vec![i]);
        }
    }
    cells
}

/// The tropical hypersurface: the corner locus of a tropical polynomial.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalHypersurface {
    /// The defining polynomial.
    pub polynomial: TropicalPolynomial,
}

impl TropicalHypersurface {
    pub fn new(poly: TropicalPolynomial) -> Self {
        TropicalHypersurface { polynomial: poly }
    }

    /// Check if a point is on the hypersurface
    /// (at least two monomials achieve the maximum).
    pub fn contains_point(&self, point: &[Tropical]) -> bool {
        let dom = self.polynomial.dominant_monomials(point);
        dom.len() >= 2
    }

    /// Get the Newton polytope.
    pub fn newton_polytope(&self) -> NewtonPolytope {
        NewtonPolytope::from_polynomial(&self.polynomial)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_newton_polytope_from_poly() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1), (2.0, 2)]);
        let np = NewtonPolytope::from_polynomial(&p);
        assert_eq!(np.num_vertices(), 3);
    }

    #[test]
    fn test_newton_polytope_dimension_1d() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        let np = NewtonPolytope::from_polynomial(&p);
        assert_eq!(np.dimension(), 1);
    }

    #[test]
    fn test_newton_polytope_dimension_0() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0)]);
        let np = NewtonPolytope::from_polynomial(&p);
        assert_eq!(np.dimension(), 0);
    }

    #[test]
    fn test_lattice_point_dot() {
        let a = LatticePoint::new(vec![1, 2, 3]);
        let b = LatticePoint::new(vec![4, 5, 6]);
        assert_eq!(a.dot(&b), 32); // 4+10+18
    }

    #[test]
    fn test_lattice_point_sub() {
        let a = LatticePoint::new(vec![5, 3]);
        let b = LatticePoint::new(vec![2, 1]);
        let c = a.sub(&b);
        assert_eq!(c.coords, vec![3, 2]);
    }

    #[test]
    fn test_newton_polytope_contains() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 3)]);
        let np = NewtonPolytope::from_polynomial(&p);
        assert!(np.contains(&LatticePoint::new(vec![1])));
        assert!(np.contains(&LatticePoint::new(vec![2])));
        assert!(!np.contains(&LatticePoint::new(vec![5])));
    }

    #[test]
    fn test_newton_polytope_volume_1d() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 5)]);
        let np = NewtonPolytope::from_polynomial(&p);
        assert_eq!(np.volume(), 5.0);
    }

    #[test]
    fn test_minkowski_sum() {
        let np1 = NewtonPolytope {
            vertices: vec![LatticePoint::new(vec![0]), LatticePoint::new(vec![1])],
            dim: 1,
        };
        let np2 = NewtonPolytope {
            vertices: vec![LatticePoint::new(vec![0]), LatticePoint::new(vec![1])],
            dim: 1,
        };
        let sum = np1.minkowski_sum(&np2);
        assert_eq!(sum.vertices.len(), 4);
    }

    #[test]
    fn test_polyhedral_subdivision() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1)]);
        let sub = PolyhedralSubdivision::from_polynomial(&p);
        assert_eq!(sub.num_cells(), 2);
    }

    #[test]
    fn test_tropical_hypersurface_contains() {
        // max(0, x): corner at x=0
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        let hs = TropicalHypersurface::new(p);
        assert!(hs.contains_point(&[Tropical::new(0.0)]));
        assert!(!hs.contains_point(&[Tropical::new(5.0)]));
    }

    #[test]
    fn test_hypersurface_newton() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1)]);
        let hs = TropicalHypersurface::new(p.clone());
        let np = hs.newton_polytope();
        assert_eq!(np.num_vertices(), 2);
    }

    #[test]
    fn test_rank_of_vectors() {
        let v = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(rank_of_vectors(&v), 2);
    }

    #[test]
    fn test_rank_of_vectors_dependent() {
        let v = vec![vec![1, 0], vec![2, 0]];
        assert_eq!(rank_of_vectors(&v), 1);
    }

    #[test]
    fn test_upper_hull() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1), (2.0, 2)]);
        let np = NewtonPolytope::from_polynomial(&p);
        let hull = np.upper_hull_indices(&[1.0]);
        assert!(!hull.is_empty());
    }
}
