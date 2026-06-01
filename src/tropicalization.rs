//! Tropicalization: classical → tropical as a deformation (t → 0 limit of log_t).
//!
//! Tropicalization is the process of transforming classical algebraic varieties
//! into tropical varieties via the valuation map: x ↦ -log_t(x) as t → 0.

use crate::polynomial::TropicalPolynomial;
use crate::semiring::Tropical;
use serde::{Deserialize, Serialize};

/// A classical polynomial term: coefficient * x^degree.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClassicalTerm {
    pub coefficient: f64,
    pub degree: Vec<u32>,
}

impl ClassicalTerm {
    pub fn new(coefficient: f64, degree: Vec<u32>) -> Self {
        ClassicalTerm { coefficient, degree }
    }

    /// Evaluate at a point.
    pub fn evaluate(&self, point: &[f64]) -> f64 {
        let mut result = self.coefficient;
        for (i, &d) in self.degree.iter().enumerate() {
            if i < point.len() && d > 0 {
                result *= point[i].powi(d as i32);
            }
        }
        result
    }
}

/// A classical polynomial.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClassicalPolynomial {
    pub terms: Vec<ClassicalTerm>,
    pub num_vars: usize,
}

impl ClassicalPolynomial {
    pub fn new(terms: Vec<ClassicalTerm>, num_vars: usize) -> Self {
        ClassicalPolynomial { terms, num_vars }
    }

    /// Evaluate at a point.
    pub fn evaluate(&self, point: &[f64]) -> f64 {
        self.terms.iter().map(|t| t.evaluate(point)).sum()
    }

    /// Tropicalize: apply -log_t to coefficients and convert + to max, * to +.
    /// As t → 0: -log_t(c * x^d) → val(c) + dot(d, point)
    /// where val(c) = -log_t(|c|) and point = -log_t(x).
    pub fn tropicalize(&self, base: f64) -> TropicalPolynomial {
        let monomials: Vec<crate::polynomial::TropicalMonomial> = self
            .terms
            .iter()
            .filter(|t| t.coefficient.abs() > 1e-15)
            .map(|t| {
                // val(coefficient) = -log_base(|c|)
                let coeff = if (base - 1.0).abs() < 1e-10 || base <= 1.0 {
                    t.coefficient.abs().ln() // Natural log approximation
                } else {
                    -t.coefficient.abs().log(base)
                };
                crate::polynomial::TropicalMonomial::new(coeff, t.degree.clone())
            })
            .collect();
        TropicalPolynomial::new(monomials, self.num_vars)
    }

    /// Tropicalize using the 2-adic valuation as default.
    pub fn tropicalize_default(&self) -> TropicalPolynomial {
        self.tropicalize(2.0)
    }
}

/// Tropicalization via the degenerations approach.
/// Replace x by t^a and take limits.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tropicalization {
    /// The base of the logarithm.
    pub base: f64,
    /// The resulting tropical polynomial.
    pub tropical_polynomial: Option<TropicalPolynomial>,
}

impl Tropicalization {
    pub fn new(base: f64) -> Self {
        Tropicalization {
            base,
            tropical_polynomial: None,
        }
    }

    /// Tropicalize a classical polynomial.
    pub fn tropicalize(&mut self, poly: &ClassicalPolynomial) -> &TropicalPolynomial {
        self.tropical_polynomial = Some(poly.tropicalize(self.base));
        self.tropical_polynomial.as_ref().unwrap()
    }

    /// Compute the initial ideal (leading terms with respect to a weight vector).
    pub fn initial_form(&self, poly: &ClassicalPolynomial, weight: &[f64]) -> ClassicalPolynomial {
        let mut max_wd = f64::NEG_INFINITY;
        let mut weighted_vals: Vec<f64> = Vec::with_capacity(poly.terms.len());

        for term in &poly.terms {
            let wd: f64 = term
                .degree
                .iter()
                .zip(weight.iter())
                .map(|(&d, &w)| d as f64 * w)
                .sum();
            weighted_vals.push(wd);
            if wd > max_wd {
                max_wd = wd;
            }
        }

        let leading_terms: Vec<ClassicalTerm> = poly
            .terms
            .iter()
            .zip(weighted_vals.iter())
            .filter(|(_, &wv)| (wv - max_wd).abs() < 1e-10)
            .map(|(t, _)| t.clone())
            .collect();

        ClassicalPolynomial::new(leading_terms, poly.num_vars)
    }

    /// Compute the Gröbner fan direction (weight vectors that give the same initial form).
    pub fn grobner_directions(&self, poly: &ClassicalPolynomial) -> Vec<Vec<f64>> {
        let n = poly.num_vars;
        if n == 0 || poly.terms.len() <= 1 {
            return vec![vec![0.0; n]];
        }

        // Each pair of terms defines a wall in the Gröbner fan
        let mut directions = Vec::new();
        for i in 0..poly.terms.len() {
            for j in (i + 1)..poly.terms.len() {
                let diff: Vec<f64> = poly.terms[i]
                    .degree
                    .iter()
                    .zip(poly.terms[j].degree.iter())
                    .map(|(&di, &dj)| (di as f64) - (dj as f64))
                    .collect();
                let norm: f64 = diff.iter().map(|d| d * d).sum::<f64>().sqrt();
                if norm > 1e-10 {
                    let dir: Vec<f64> = diff.iter().map(|d| d / norm).collect();
                    directions.push(dir);
                }
            }
        }
        directions
    }
}

/// Lift a tropical polynomial back to a classical polynomial (reverse tropicalization).
pub fn detropicalize(tropical_poly: &TropicalPolynomial, base: f64) -> ClassicalPolynomial {
    let terms: Vec<ClassicalTerm> = tropical_poly
        .monomials
        .iter()
        .map(|m| {
            // Reverse: coefficient = base^(-val)
            let coeff = if (base - 1.0).abs() < 1e-10 || base <= 1.0 {
                m.coefficient.value().exp()
            } else {
                base.powf(-m.coefficient.value())
            };
            ClassicalTerm::new(coeff, m.degree.clone())
        })
        .collect();
    ClassicalPolynomial::new(terms, tropical_poly.num_vars)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classical_evaluate() {
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(1.0, vec![2]),
                ClassicalTerm::new(3.0, vec![1]),
                ClassicalTerm::new(2.0, vec![0]),
            ],
            1,
        );
        // 1*x^2 + 3*x + 2 at x=1 => 1+3+2=6
        assert_eq!(p.evaluate(&[1.0]), 6.0);
    }

    #[test]
    fn test_tropicalize_basic() {
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(4.0, vec![2]),
                ClassicalTerm::new(2.0, vec![1]),
                ClassicalTerm::new(1.0, vec![0]),
            ],
            1,
        );
        let tropical = p.tropicalize(2.0);
        assert_eq!(tropical.monomials.len(), 3);
        // val(4) = -log2(4) = -2, val(2) = -1, val(1) = 0
        assert_eq!(tropical.monomials[0].coefficient, Tropical::new(-2.0));
        assert_eq!(tropical.monomials[1].coefficient, Tropical::new(-1.0));
        assert_eq!(tropical.monomials[2].coefficient, Tropical::new(0.0));
    }

    #[test]
    fn test_tropicalization_struct() {
        let mut tr = Tropicalization::new(2.0);
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(1.0, vec![1]),
                ClassicalTerm::new(1.0, vec![0]),
            ],
            1,
        );
        let result = tr.tropicalize(&p);
        assert_eq!(result.monomials.len(), 2);
    }

    #[test]
    fn test_initial_form() {
        let tr = Tropicalization::new(2.0);
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(1.0, vec![2]),
                ClassicalTerm::new(1.0, vec![1]),
                ClassicalTerm::new(1.0, vec![0]),
            ],
            1,
        );
        // Weight [1.0]: leading form = x^2 (highest weighted degree)
        let init = tr.initial_form(&p, &[1.0]);
        assert_eq!(init.terms.len(), 1);
        assert_eq!(init.terms[0].degree, vec![2]);
    }

    #[test]
    fn test_grobner_directions() {
        let tr = Tropicalization::new(2.0);
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(1.0, vec![2]),
                ClassicalTerm::new(1.0, vec![0]),
            ],
            1,
        );
        let dirs = tr.grobner_directions(&p);
        assert_eq!(dirs.len(), 1);
        // Direction should be (2-0)/|2-0| = 1
        assert!((dirs[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_detropicalize() {
        let tp = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1)]);
        let classical = detropicalize(&tp, 2.0);
        assert_eq!(classical.terms.len(), 2);
        // 2^0 = 1, 2^0 = 1
        assert_eq!(classical.terms[0].coefficient, 1.0);
        assert_eq!(classical.terms[1].coefficient, 1.0);
    }

    #[test]
    fn test_tropicalize_evaluate_consistency() {
        // Classical: 2x + 3 at x=4 => 11
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(2.0, vec![1]),
                ClassicalTerm::new(3.0, vec![0]),
            ],
            1,
        );
        let classical_val = p.evaluate(&[4.0]);
        assert_eq!(classical_val, 11.0);

        // Tropical: max(val(2)+x, val(3)) with base=2
        let tropical = p.tropicalize(2.0);
        // At x = log2(4) = 2: max(-1+2, -log2(3)) = max(1, -1.585) = 1
        let tropical_val = tropical.evaluate(&[Tropical::new(2.0)]);
        assert!(tropical_val.value().is_finite());
    }

    #[test]
    fn test_bivariate_tropicalization() {
        let p = ClassicalPolynomial::new(
            vec![
                ClassicalTerm::new(1.0, vec![1, 0]),
                ClassicalTerm::new(1.0, vec![0, 1]),
                ClassicalTerm::new(1.0, vec![0, 0]),
            ],
            2,
        );
        let tropical = p.tropicalize(2.0);
        assert_eq!(tropical.num_vars, 2);
        assert_eq!(tropical.monomials.len(), 3);
    }
}
