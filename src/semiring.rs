//! Tropical semiring: (ℝ ∪ {-∞}, max, +)
//!
//! The fundamental algebraic structure. Tropical addition is max,
//! tropical multiplication is ordinary addition.

use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fmt;

/// The tropical number type wrapping f64.
/// Represents elements of the max-plus semiring (ℝ ∪ {-∞}, max, +).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Tropical(pub f64);

impl Tropical {
    /// Tropical zero element: -∞ (additive identity for max)
    pub const NEG_INF: Tropical = Tropical(f64::NEG_INFINITY);

    /// Tropical one (multiplicative identity): 0.0
    pub const ONE: Tropical = Tropical(0.0);

    /// Tropical zero (additive identity)
    pub const ZERO: Tropical = Tropical(f64::NEG_INFINITY);

    /// Create a tropical number from a regular f64.
    pub fn new(val: f64) -> Self {
        Tropical(val)
    }

    /// Create from an integer.
    pub fn from_int(val: i64) -> Self {
        Tropical(val as f64)
    }

    /// Tropical addition = max.
    #[inline]
    pub fn tropical_add(self, other: Self) -> Self {
        Tropical(self.0.max(other.0))
    }

    /// Tropical multiplication = ordinary addition.
    #[inline]
    pub fn tropical_mul(self, other: Self) -> Self {
        Tropical(self.0 + other.0)
    }

    /// Tropical subtraction (inverse of tropical multiplication).
    #[inline]
    pub fn tropical_sub(self, other: Self) -> Self {
        Tropical(self.0 - other.0)
    }

    /// Tropical division (inverse of tropical multiplication).
    pub fn tropical_div(self, other: Self) -> Option<Self> {
        if other.0.is_infinite() && other.0.is_sign_negative() {
            None // Division by -∞ undefined
        } else {
            Some(Tropical(self.0 - other.0))
        }
    }

    /// Tropical exponentiation: repeated tropical multiplication = scalar multiplication.
    pub fn tropical_pow(self, n: u32) -> Self {
        if n == 0 {
            Tropical::ONE
        } else {
            Tropical(self.0 * n as f64)
        }
    }

    /// Check if this is the zero element (-∞).
    pub fn is_zero(&self) -> bool {
        self.0.is_infinite() && self.0.is_sign_negative()
    }

    /// Check if this is the one element (0.0).
    pub fn is_one(&self) -> bool {
        self.0 == 0.0
    }

    /// Get the inner f64 value.
    pub fn value(&self) -> f64 {
        self.0
    }

    /// Tropical absolute value (same as regular absolute value in max-plus).
    pub fn tropical_abs(self) -> Self {
        // In the tropical semiring, "absolute value" maps to regular abs
        Tropical(self.0.abs())
    }

    /// Tropical negation: returns -self (relevant for min-plus).
    pub fn tropical_neg(self) -> Self {
        Tropical(-self.0)
    }
}

impl Default for Tropical {
    fn default() -> Self {
        Tropical::ZERO
    }
}

impl Eq for Tropical {}

impl PartialOrd for Tropical {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Tropical {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.partial_cmp(&other.0).unwrap_or(Ordering::Equal)
    }
}

impl std::ops::Add for Tropical {
    type Output = Tropical;
    /// Tropical addition = max.
    fn add(self, rhs: Self) -> Self::Output {
        self.tropical_add(rhs)
    }
}

impl std::ops::Mul for Tropical {
    type Output = Tropical;
    /// Tropical multiplication = ordinary addition.
    fn mul(self, rhs: Self) -> Self::Output {
        self.tropical_mul(rhs)
    }

}

impl std::ops::Sub for Tropical {
    type Output = Tropical;
    fn sub(self, rhs: Self) -> Self::Output {
        self.tropical_sub(rhs)
    }
}

impl std::iter::Sum for Tropical {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Tropical::ZERO, |a, b| a + b)
    }
}

impl std::iter::Product for Tropical {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Tropical::ONE, |a, b| a * b)
    }
}

impl fmt::Display for Tropical {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            write!(f, "-∞")
        } else {
            write!(f, "{}", self.0)
        }
    }
}

/// Trait for types that support tropical arithmetic.
pub trait TropicalOps: Clone {
    fn tropical_add(&self, other: &Self) -> Self;
    fn tropical_mul(&self, other: &Self) -> Self;
    fn tropical_zero() -> Self;
    fn tropical_one() -> Self;
}

impl TropicalOps for Tropical {
    fn tropical_add(&self, other: &Self) -> Self {
        *self + *other
    }
    fn tropical_mul(&self, other: &Self) -> Self {
        *self * *other
    }
    fn tropical_zero() -> Self {
        Tropical::ZERO
    }
    fn tropical_one() -> Self {
        Tropical::ONE
    }
}

/// Min-plus semiring: (ℝ ∪ {+∞}, min, +).
/// Dual of the max-plus semiring.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MinPlus(pub f64);

impl MinPlus {
    /// Min-plus zero element: +∞ (additive identity for min)
    pub const POS_INF: MinPlus = MinPlus(f64::INFINITY);
    /// Min-plus one (multiplicative identity): 0.0
    pub const ONE: MinPlus = MinPlus(0.0);

    pub fn new(val: f64) -> Self {
        MinPlus(val)
    }

    /// Min-plus addition = min.
    pub fn min_add(self, other: Self) -> Self {
        MinPlus(self.0.min(other.0))
    }

    /// Min-plus multiplication = ordinary addition.
    pub fn min_mul(self, other: Self) -> Self {
        MinPlus(self.0 + other.0)
    }

    pub fn value(&self) -> f64 {
        self.0
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_infinite() && self.0.is_sign_positive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tropical_add_basic() {
        let a = Tropical::new(3.0);
        let b = Tropical::new(5.0);
        assert_eq!(a + b, Tropical::new(5.0));
    }

    #[test]
    fn test_tropical_add_commutative() {
        let a = Tropical::new(3.0);
        let b = Tropical::new(7.0);
        assert_eq!(a + b, b + a);
    }

    #[test]
    fn test_tropical_add_associative() {
        let a = Tropical::new(1.0);
        let b = Tropical::new(2.0);
        let c = Tropical::new(3.0);
        assert_eq!((a + b) + c, a + (b + c));
    }

    #[test]
    fn test_tropical_mul_basic() {
        let a = Tropical::new(3.0);
        let b = Tropical::new(5.0);
        assert_eq!(a * b, Tropical::new(8.0));
    }

    #[test]
    fn test_tropical_mul_commutative() {
        let a = Tropical::new(3.0);
        let b = Tropical::new(7.0);
        assert_eq!(a * b, b * a);
    }

    #[test]
    fn test_tropical_mul_associative() {
        let a = Tropical::new(1.0);
        let b = Tropical::new(2.0);
        let c = Tropical::new(3.0);
        assert_eq!((a * b) * c, a * (b * c));
    }

    #[test]
    fn test_tropical_distributive() {
        let a = Tropical::new(1.0);
        let b = Tropical::new(2.0);
        let c = Tropical::new(3.0);
        // a * (b + c) = a * b + a * c (tropical distributivity)
        let left = a * (b + c);
        let right = (a * b) + (a * c);
        assert_eq!(left, right);
    }

    #[test]
    fn test_tropical_zero_identity_add() {
        let a = Tropical::new(5.0);
        assert_eq!(a + Tropical::ZERO, a);
        assert_eq!(Tropical::ZERO + a, a);
    }

    #[test]
    fn test_tropical_one_identity_mul() {
        let a = Tropical::new(5.0);
        assert_eq!(a * Tropical::ONE, a);
        assert_eq!(Tropical::ONE * a, a);
    }

    #[test]
    fn test_tropical_zero_annihilates_mul() {
        let a = Tropical::new(5.0);
        assert_eq!(a * Tropical::ZERO, Tropical::ZERO);
    }

    #[test]
    fn test_tropical_pow() {
        let a = Tropical::new(3.0);
        assert_eq!(a.tropical_pow(0), Tropical::ONE);
        assert_eq!(a.tropical_pow(1), Tropical::new(3.0));
        assert_eq!(a.tropical_pow(2), Tropical::new(6.0));
        assert_eq!(a.tropical_pow(3), Tropical::new(9.0));
    }

    #[test]
    fn test_tropical_idempotent() {
        let a = Tropical::new(5.0);
        assert_eq!(a + a, a);
    }

    #[test]
    fn test_tropical_sub() {
        let a = Tropical::new(8.0);
        let b = Tropical::new(3.0);
        assert_eq!(a - b, Tropical::new(5.0));
    }

    #[test]
    fn test_tropical_div() {
        let a = Tropical::new(8.0);
        let b = Tropical::new(3.0);
        assert_eq!(a.tropical_div(b), Some(Tropical::new(5.0)));
    }

    #[test]
    fn test_tropical_div_by_zero() {
        let a = Tropical::new(5.0);
        assert_eq!(a.tropical_div(Tropical::ZERO), None);
    }

    #[test]
    fn test_tropical_ordering() {
        let a = Tropical::new(3.0);
        let b = Tropical::new(5.0);
        assert!(a < b);
        assert!(b > a);
    }

    #[test]
    fn test_tropical_neg() {
        let a = Tropical::new(3.0);
        assert_eq!(a.tropical_neg(), Tropical::new(-3.0));
    }

    #[test]
    fn test_tropical_is_zero() {
        assert!(Tropical::ZERO.is_zero());
        assert!(!Tropical::new(5.0).is_zero());
    }

    #[test]
    fn test_tropical_is_one() {
        assert!(Tropical::ONE.is_one());
        assert!(!Tropical::new(5.0).is_one());
    }

    #[test]
    fn test_tropical_sum_iterator() {
        let vals = vec![Tropical::new(1.0), Tropical::new(3.0), Tropical::new(2.0)];
        let sum: Tropical = vals.into_iter().sum();
        assert_eq!(sum, Tropical::new(3.0)); // max(1, 3, 2)
    }

    #[test]
    fn test_tropical_product_iterator() {
        let vals = vec![Tropical::new(1.0), Tropical::new(3.0), Tropical::new(2.0)];
        let product: Tropical = vals.into_iter().product();
        assert_eq!(product, Tropical::new(6.0)); // 1 + 3 + 2
    }

    #[test]
    fn test_tropical_display() {
        assert_eq!(format!("{}", Tropical::new(5.0)), "5");
        assert_eq!(format!("{}", Tropical::ZERO), "-∞");
    }

    #[test]
    fn test_tropical_serde() {
        let a = Tropical::new(3.14);
        let json = serde_json::to_string(&a).unwrap();
        let b: Tropical = serde_json::from_str(&json).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn test_min_plus_basic() {
        let a = MinPlus::new(3.0);
        let b = MinPlus::new(5.0);
        assert_eq!(a.min_add(b), MinPlus::new(3.0)); // min(3, 5)
        assert_eq!(a.min_mul(b), MinPlus::new(8.0)); // 3 + 5
    }

    #[test]
    fn test_min_plus_zero() {
        let a = MinPlus::new(5.0);
        assert_eq!(a.min_add(MinPlus::POS_INF), a);
    }
}
