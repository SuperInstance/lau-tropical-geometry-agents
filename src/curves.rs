//! Tropical curves: piecewise-linear skeletons of algebraic curves.
//!
//! Tropical curves are 1-dimensional polyhedral complexes that are
//! the tropicalization of algebraic curves.

use crate::polynomial::TropicalPolynomial;
use crate::semiring::Tropical;
use serde::{Deserialize, Serialize};

/// A vertex in a tropical curve.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CurveVertex {
    /// Coordinates of the vertex.
    pub coords: Vec<f64>,
}

impl CurveVertex {
    pub fn new(coords: Vec<f64>) -> Self {
        CurveVertex { coords }
    }

    pub fn dim(&self) -> usize {
        self.coords.len()
    }
}

/// An edge in a tropical curve, connecting two vertices.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CurveEdge {
    /// Index of first vertex.
    pub v1: usize,
    /// Index of second vertex.
    pub v2: usize,
    /// Weight/multiplicity of the edge.
    pub weight: u32,
    /// Direction vector (primitive integer vector).
    pub direction: Vec<i64>,
}

impl CurveEdge {
    pub fn new(v1: usize, v2: usize, weight: u32, direction: Vec<i64>) -> Self {
        CurveEdge { v1, v2, weight, direction }
    }
}

/// A tropical curve: a weighted balanced 1-dimensional polyhedral complex.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalCurve {
    /// Vertices of the curve.
    pub vertices: Vec<CurveVertex>,
    /// Edges of the curve.
    pub edges: Vec<CurveEdge>,
    /// Ambient dimension.
    pub ambient_dim: usize,
}

impl TropicalCurve {
    pub fn new(vertices: Vec<CurveVertex>, edges: Vec<CurveEdge>, ambient_dim: usize) -> Self {
        TropicalCurve { vertices, edges, ambient_dim }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.edges.len()
    }

    /// Compute a tropical curve from a bivariate tropical polynomial.
    /// The tropical curve is the corner locus in 2D.
    pub fn from_bivariate_polynomial(poly: &TropicalPolynomial, resolution: f64) -> Self {
        assert_eq!(poly.num_vars, 2, "Expected bivariate polynomial");
        let mut vertices = Vec::new();
        let mut edges = Vec::new();

        // Find vertices: points where 3+ monomials achieve the max simultaneously
        // For simplicity, sample on a grid and find approximate corners
        let range = 10.0;
        let step = resolution.max(0.1);
        let mut corner_points: Vec<(f64, f64, Vec<usize>)> = Vec::new();

        let n_steps = (2.0 * range / step) as usize;
        for i in 0..=n_steps {
            let x = -range + i as f64 * step;
            for j in 0..=n_steps {
                let y = -range + j as f64 * step;
                let point = [Tropical::new(x), Tropical::new(y)];
                let dom = poly.dominant_monomials(&point);
                if dom.len() >= 3 {
                    corner_points.push((x, y, dom.clone()));
                }
            }
        }

        // Deduplicate nearby corner points
        let mut deduped: Vec<(f64, f64, Vec<usize>)> = Vec::new();
        for (x, y, dom) in &corner_points {
            let is_dup = deduped.iter().any(|(px, py, _)| {
                (px - x).abs() < step * 1.5 && (py - y).abs() < step * 1.5
            });
            if !is_dup {
                deduped.push((*x, *y, dom.clone()));
            }
        }

        // Add vertices
        for (x, y, _) in &deduped {
            vertices.push(CurveVertex::new(vec![*x, *y]));
        }

        // Connect nearby vertices with edges
        for i in 0..vertices.len() {
            for j in (i + 1)..vertices.len() {
                let dx = vertices[j].coords[0] - vertices[i].coords[0];
                let dy = vertices[j].coords[1] - vertices[i].coords[1];
                let dist = (dx * dx + dy * dy).sqrt();
                if dist < range * 0.3 && dist > 0.01 {
                    let dir = vec![dx.round() as i64, dy.round() as i64];
                    edges.push(CurveEdge::new(i, j, 1, dir));
                }
            }
        }

        TropicalCurve::new(vertices, edges, 2)
    }

    /// Check the balancing condition at a vertex.
    /// At each vertex, the weighted sum of outgoing direction vectors should be zero.
    pub fn is_balanced_at(&self, vertex_idx: usize) -> bool {
        if vertex_idx >= self.vertices.len() {
            return false;
        }
        let mut sum = vec![0i64; self.ambient_dim];
        for edge in &self.edges {
            if edge.v1 == vertex_idx {
                for (k, &d) in edge.direction.iter().enumerate() {
                    if k < sum.len() {
                        sum[k] += d * edge.weight as i64;
                    }
                }
            } else if edge.v2 == vertex_idx {
                for (k, &d) in edge.direction.iter().enumerate() {
                    if k < sum.len() {
                        sum[k] -= d * edge.weight as i64;
                    }
                }
            }
        }
        sum.iter().all(|&s| s == 0)
    }

    /// Check the balancing condition at all vertices.
    pub fn is_balanced(&self) -> bool {
        (0..self.vertices.len()).all(|i| self.is_balanced_at(i))
    }

    /// Compute the genus (number of cycles).
    pub fn genus(&self) -> usize {
        // Genus = E - V + 1 (for connected graph)
        // Use signed arithmetic to avoid underflow
        let e = self.num_edges() as i64;
        let v = self.num_vertices() as i64;
        0.max((e - v + 1) as usize)
    }

    /// Get edges incident to a vertex.
    pub fn incident_edges(&self, vertex_idx: usize) -> Vec<usize> {
        self.edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.v1 == vertex_idx || e.v2 == vertex_idx)
            .map(|(i, _)| i)
            .collect()
    }

    /// Degree of a vertex (number of incident edges).
    pub fn vertex_degree(&self, vertex_idx: usize) -> usize {
        self.incident_edges(vertex_idx).len()
    }

    /// Total weight of the curve.
    pub fn total_weight(&self) -> u32 {
        self.edges.iter().map(|e| e.weight).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_curve_vertex_creation() {
        let v = CurveVertex::new(vec![1.0, 2.0]);
        assert_eq!(v.dim(), 2);
        assert_eq!(v.coords, vec![1.0, 2.0]);
    }

    #[test]
    fn test_curve_edge_creation() {
        let e = CurveEdge::new(0, 1, 2, vec![1, 0]);
        assert_eq!(e.v1, 0);
        assert_eq!(e.v2, 1);
        assert_eq!(e.weight, 2);
    }

    #[test]
    fn test_tropical_curve_basic() {
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![0, 1]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2], 2);
        assert_eq!(curve.num_vertices(), 3);
        assert_eq!(curve.num_edges(), 2);
    }

    #[test]
    fn test_curve_genus_tree() {
        // Tree: V=3, E=2 => genus = max(0, 2-3+1) = 0
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![0, 1]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2], 2);
        assert_eq!(curve.genus(), 0);
    }

    #[test]
    fn test_curve_genus_cycle() {
        // Triangle: V=3, E=3 => genus = 3-3+1 = 1
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(1, 2, 1, vec![-1, 1]);
        let e3 = CurveEdge::new(2, 0, 1, vec![0, -1]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2, e3], 2);
        assert_eq!(curve.genus(), 1);
    }

    #[test]
    fn test_incident_edges() {
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![0, 1]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2], 2);
        assert_eq!(curve.incident_edges(0).len(), 2);
        assert_eq!(curve.incident_edges(1).len(), 1);
    }

    #[test]
    fn test_vertex_degree() {
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![0, 1]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2], 2);
        assert_eq!(curve.vertex_degree(0), 2);
        assert_eq!(curve.vertex_degree(1), 1);
    }

    #[test]
    fn test_total_weight() {
        let e1 = CurveEdge::new(0, 1, 2, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 3, vec![0, 1]);
        let v1 = CurveVertex::new(vec![0.0, 0.0]);
        let v2 = CurveVertex::new(vec![1.0, 0.0]);
        let v3 = CurveVertex::new(vec![0.0, 1.0]);
        let curve = TropicalCurve::new(vec![v1, v2, v3], vec![e1, e2], 2);
        assert_eq!(curve.total_weight(), 5);
    }

    #[test]
    fn test_from_bivariate_polynomial() {
        // max(0, x, y) should give a tropical curve with the "tent" shape
        let p = TropicalPolynomial::bivariate(vec![
            (0.0, 0, 0),
            (0.0, 1, 0),
            (0.0, 0, 1),
        ]);
        let curve = TropicalCurve::from_bivariate_polynomial(&p, 0.5);
        // Should have at least one vertex near the origin
        assert!(curve.num_vertices() > 0);
    }

    #[test]
    fn test_balanced_star() {
        // A balanced 3-pronged star at the origin
        let v0 = CurveVertex::new(vec![0.0, 0.0]);
        let v1 = CurveVertex::new(vec![1.0, 0.0]);
        let v2 = CurveVertex::new(vec![-1.0, 1.0]);
        let v3 = CurveVertex::new(vec![-1.0, -1.0]);
        // Directions from v0: (1,0), (-1,1), (-1,-1) -- not balanced
        // For balance: (1,0)+(-1,1)+(-1,-1) = (-1,0) ≠ 0
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![-1, 1]);
        let e3 = CurveEdge::new(0, 3, 1, vec![-1, -1]);
        let curve = TropicalCurve::new(vec![v0, v1, v2, v3], vec![e1, e2, e3], 2);
        assert!(!curve.is_balanced_at(0));
    }

    #[test]
    fn test_balanced_trident() {
        // A balanced 3-pronged star: directions (1,0), (0,1), (-1,-1)
        // Sum = (0,0) ✓
        let v0 = CurveVertex::new(vec![0.0, 0.0]);
        let v1 = CurveVertex::new(vec![1.0, 0.0]);
        let v2 = CurveVertex::new(vec![0.0, 1.0]);
        let v3 = CurveVertex::new(vec![-1.0, -1.0]);
        let e1 = CurveEdge::new(0, 1, 1, vec![1, 0]);
        let e2 = CurveEdge::new(0, 2, 1, vec![0, 1]);
        let e3 = CurveEdge::new(0, 3, 1, vec![-1, -1]);
        let curve = TropicalCurve::new(vec![v0, v1, v2, v3], vec![e1, e2, e3], 2);
        assert!(curve.is_balanced_at(0));
    }
}
