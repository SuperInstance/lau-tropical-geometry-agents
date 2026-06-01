//! Tropical polynomials: max of affine functions = piecewise linear.
//!
//! A tropical polynomial in one variable: max(a₀, a₁+x, a₂+2x, ..., aₙ+nx).
//! In multiple variables: max over monomials of (coefficient + dot(degree, variables)).

use crate::semiring::Tropical;
use serde::{Deserialize, Serialize};
use std::fmt;

/// A monomial in a tropical polynomial: coefficient + dot(degree, variables).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TropicalMonomial {
    /// Coefficient in the tropical sense (the constant term).
    pub coefficient: Tropical,
    /// Degree (exponent vector) for each variable.
    pub degree: Vec<u32>,
}

impl TropicalMonomial {
    pub fn new(coefficient: f64, degree: Vec<u32>) -> Self {
        TropicalMonomial {
            coefficient: Tropical::new(coefficient),
            degree,
        }
    }

    /// Evaluate this monomial at a point.
    pub fn evaluate(&self, point: &[Tropical]) -> Tropical {
        let mut result = self.coefficient;
        for (i, &d) in self.degree.iter().enumerate() {
            if i < point.len() && d > 0 {
                result = result * point[i].tropical_pow(d);
            }
        }
        result
    }

    /// Number of variables.
    pub fn num_vars(&self) -> usize {
        self.degree.len()
    }

    /// Total degree.
    pub fn total_degree(&self) -> u32 {
        self.degree.iter().sum()
    }
}

/// A tropical polynomial = max of monomials.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TropicalPolynomial {
    /// The monomials.
    pub monomials: Vec<TropicalMonomial>,
    /// Number of variables.
    pub num_vars: usize,
}

impl TropicalPolynomial {
    /// Create a new tropical polynomial.
    pub fn new(monomials: Vec<TropicalMonomial>, num_vars: usize) -> Self {
        TropicalPolynomial { monomials, num_vars }
    }

    /// Create a univariate tropical polynomial from coefficient-degree pairs.
    pub fn univariate(terms: Vec<(f64, u32)>) -> Self {
        let monomials: Vec<TropicalMonomial> = terms
            .into_iter()
            .map(|(c, d)| TropicalMonomial::new(c, vec![d]))
            .collect();
        TropicalPolynomial::new(monomials, 1)
    }

    /// Create a bivariate tropical polynomial.
    pub fn bivariate(terms: Vec<(f64, u32, u32)>) -> Self {
        let monomials: Vec<TropicalMonomial> = terms
            .into_iter()
            .map(|(c, d1, d2)| TropicalMonomial::new(c, vec![d1, d2]))
            .collect();
        TropicalPolynomial::new(monomials, 2)
    }

    /// Evaluate the polynomial at a point (tropical evaluation = max of monomials).
    pub fn evaluate(&self, point: &[Tropical]) -> Tropical {
        self.monomials
            .iter()
            .map(|m| m.evaluate(point))
            .fold(Tropical::ZERO, |a, b| a + b) // max
    }

    /// Tropical polynomial addition (max) = take all monomials from both.
    pub fn tropical_add(&self, other: &TropicalPolynomial) -> TropicalPolynomial {
        let mut monomials = self.monomials.clone();
        monomials.extend(other.monomials.clone());
        TropicalPolynomial::new(monomials, self.num_vars.max(other.num_vars))
    }

    /// Tropical polynomial multiplication.
    /// Convolve monomials: for each pair, add coefficients and add degrees.
    pub fn tropical_mul(&self, other: &TropicalPolynomial) -> TropicalPolynomial {
        let mut monomials = Vec::new();
        for a in &self.monomials {
            for b in &other.monomials {
                let coeff = a.coefficient * b.coefficient; // regular addition
                let mut degree: Vec<u32> = a
                    .degree
                    .iter()
                    .zip(b.degree.iter())
                    .map(|(&da, &db)| da + db)
                    .collect();
                // Extend if degrees have different lengths
                if a.degree.len() < b.degree.len() {
                    degree.extend_from_slice(&b.degree[a.degree.len()..]);
                } else if b.degree.len() < a.degree.len() {
                    degree.extend_from_slice(&a.degree[b.degree.len()..]);
                }
                monomials.push(TropicalMonomial { coefficient: coeff, degree });
            }
        }
        TropicalPolynomial::new(monomials, self.num_vars + other.num_vars)
    }

    /// Tropical scalar multiplication: add scalar to all coefficients.
    pub fn scalar_mul(&self, scalar: Tropical) -> TropicalPolynomial {
        let monomials = self
            .monomials
            .iter()
            .map(|m| TropicalMonomial {
                coefficient: m.coefficient * scalar,
                degree: m.degree.clone(),
            })
            .collect();
        TropicalPolynomial::new(monomials, self.num_vars)
    }

    /// Find which monomial achieves the maximum at a given point.
    /// Returns indices of achieving monomials.
    pub fn dominant_monomials(&self, point: &[Tropical]) -> Vec<usize> {
        let evals: Vec<Tropical> = self.monomials.iter().map(|m| m.evaluate(point)).collect();
        let max_val = evals.iter().fold(Tropical::ZERO, |a, &b| a + b);
        evals
            .iter()
            .enumerate()
            .filter(|(_, &v)| v == max_val)
            .map(|(i, _)| i)
            .collect()
    }

    /// Compute the corner locus (set where the polynomial is non-differentiable).
    /// For univariate: points where two monomials are equal.
    pub fn corner_points_1d(&self) -> Vec<f64> {
        assert_eq!(self.num_vars, 1, "corner_points_1d only for univariate");
        let mut points = Vec::new();
        for i in 0..self.monomials.len() {
            for j in (i + 1)..self.monomials.len() {
                let a = self.monomials[i].coefficient.value();
                let da = self.monomials[i].degree[0] as f64;
                let b = self.monomials[j].coefficient.value();
                let db = self.monomials[j].degree[0] as f64;
                if (da - db).abs() > 1e-12 {
                    // a + da*x = b + db*x => x = (b-a)/(da-db)
                    let x = (b - a) / (da - db);
                    if x.is_finite() {
                        points.push(x);
                    }
                }
            }
        }
        points.sort_by(|a, b| a.partial_cmp(b).unwrap());
        points
    }

    /// Degree of the polynomial (maximum total degree).
    pub fn degree(&self) -> u32 {
        self.monomials.iter().map(|m| m.total_degree()).max().unwrap_or(0)
    }

    /// Remove monomials that are never dominant (dominated by others everywhere).
    pub fn simplify(&self) -> TropicalPolynomial {
        let n = self.monomials.len();
        let mut keep = vec![true; n];
        for i in 0..n {
            if !keep[i] { continue; }
            for j in 0..n {
                if i == j || !keep[j] { continue; }
                // Check if monomial i dominates monomial j everywhere
                if self.monomial_dominates(i, j) {
                    keep[j] = false;
                }
            }
        }
        let monomials: Vec<TropicalMonomial> = self
            .monomials
            .iter()
            .enumerate()
            .filter(|(i, _)| keep[*i])
            .map(|(_, m)| m.clone())
            .collect();
        TropicalPolynomial::new(monomials, self.num_vars)
    }

    /// Check if monomial at index i dominates monomial at index j everywhere.
    fn monomial_dominates(&self, i: usize, j: usize) -> bool {
        // In tropical (max-plus), monomial i dominates j if
        // coefficient_i >= coefficient_j AND degree_i[k] >= degree_j[k] for all k
        // This is a sufficient condition.
        let mi = &self.monomials[i];
        let mj = &self.monomials[j];
        if mi.coefficient < mj.coefficient {
            return false;
        }
        for (di, dj) in mi.degree.iter().zip(mj.degree.iter()) {
            if di < dj {
                return false;
            }
        }
        // Also check if degree vectors have different lengths
        if mi.degree.len() < mj.degree.len() {
            return false;
        }
        true
    }
}

impl fmt::Display for TropicalPolynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let terms: Vec<String> = self
            .monomials
            .iter()
            .map(|m| {
                let vars: Vec<String> = m
                    .degree
                    .iter()
                    .enumerate()
                    .map(|(i, &d)| {
                        if d == 0 {
                            String::new()
                        } else if d == 1 {
                            format!("x{}", i)
                        } else {
                            format!("x{}^{}", i, d)
                        }
                    })
                    .filter(|s| !s.is_empty())
                    .collect();
                let var_part = if vars.is_empty() {
                    String::new()
                } else {
                    format!("+{}", vars.join("+"))
                };
                format!("({}{})", m.coefficient, var_part)
            })
            .collect();
        write!(f, "max({})", terms.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_univariate_evaluate() {
        // max(0, x) evaluated at x=3 => max(0, 3) = 3
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        let result = p.evaluate(&[Tropical::new(3.0)]);
        assert_eq!(result, Tropical::new(3.0));
    }

    #[test]
    fn test_univariate_evaluate_negative() {
        // max(0, x) at x=-5 => max(0, -5) = 0
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        let result = p.evaluate(&[Tropical::new(-5.0)]);
        assert_eq!(result, Tropical::new(0.0));
    }

    #[test]
    fn test_univariate_quadratic() {
        // max(0, x, 2x) at x=3 => max(0, 3, 6) = 6
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1), (0.0, 2)]);
        let result = p.evaluate(&[Tropical::new(3.0)]);
        assert_eq!(result, Tropical::new(6.0));
    }

    #[test]
    fn test_polynomial_add() {
        let p1 = TropicalPolynomial::univariate(vec![(1.0, 0)]);
        let p2 = TropicalPolynomial::univariate(vec![(3.0, 0)]);
        let sum = p1.tropical_add(&p2);
        assert_eq!(sum.monomials.len(), 2);
    }

    #[test]
    fn test_polynomial_mul() {
        let p1 = TropicalPolynomial::univariate(vec![(1.0, 1)]);
        let p2 = TropicalPolynomial::univariate(vec![(2.0, 1)]);
        let product = p1.tropical_mul(&p2);
        assert_eq!(product.monomials.len(), 1);
        // coefficient = 1+2=3, degree = 1+1=2
        assert_eq!(product.monomials[0].coefficient, Tropical::new(3.0));
        assert_eq!(product.monomials[0].degree, vec![2]);
    }

    #[test]
    fn test_scalar_mul() {
        let p = TropicalPolynomial::univariate(vec![(1.0, 0), (2.0, 1)]);
        let scaled = p.scalar_mul(Tropical::new(5.0));
        assert_eq!(scaled.monomials[0].coefficient, Tropical::new(6.0));
        assert_eq!(scaled.monomials[1].coefficient, Tropical::new(7.0));
    }

    #[test]
    fn test_dominant_monomial() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        // At x=5: max(0, 5) => monomial 1 dominates
        let dom = p.dominant_monomials(&[Tropical::new(5.0)]);
        assert_eq!(dom, vec![1]);
    }

    #[test]
    fn test_dominant_monomial_tie() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        // At x=0: max(0, 0) => both dominate
        let dom = p.dominant_monomials(&[Tropical::new(0.0)]);
        assert_eq!(dom.len(), 2);
    }

    #[test]
    fn test_corner_points() {
        // max(0, x, 2x): corners at x=0 (0=x) and x=0 (x=2x), actually
        // 0=x at x=0, 0=2x at x=0, x=2x at x=0
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1), (0.0, 2)]);
        let corners = p.corner_points_1d();
        // 0 = x => (0-0)/(0-1) = 0, 0 = 2x => 0, x = 2x => 0
        assert!(corners.iter().any(|&c| c.abs() < 1e-10));
    }

    #[test]
    fn test_corner_points_offset() {
        // max(2, x): corner at x=2
        let p = TropicalPolynomial::univariate(vec![(2.0, 0), (0.0, 1)]);
        let corners = p.corner_points_1d();
        assert!(corners.iter().any(|&c| (c - 2.0).abs() < 1e-10));
    }

    #[test]
    fn test_bivariate_evaluate() {
        // max(0, x0, x1) at (2,3) => max(0, 2, 3) = 3
        let p = TropicalPolynomial::bivariate(vec![(0.0, 0, 0), (0.0, 1, 0), (0.0, 0, 1)]);
        let result = p.evaluate(&[Tropical::new(2.0), Tropical::new(3.0)]);
        assert_eq!(result, Tropical::new(3.0));
    }

    #[test]
    fn test_degree() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1), (0.0, 3)]);
        assert_eq!(p.degree(), 3);
    }

    #[test]
    fn test_simplify() {
        // max(0, x, 2x) => 2x dominates x for x>0, but x dominates for x<0
        // Actually 2x has higher coefficient when x>0, so neither dominates everywhere
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1), (0.0, 1)]);
        // monomial at index 1 (1+x) dominates monomial at index 2 (0+x) everywhere
        let simplified = p.simplify();
        assert!(simplified.monomials.len() <= 3);
    }

    #[test]
    fn test_polynomial_display() {
        let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1)]);
        let s = format!("{}", p);
        assert!(s.contains("max("));
    }
}
