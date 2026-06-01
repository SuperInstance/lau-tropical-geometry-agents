//! # lau-tropical-geometry-agents
//!
//! Tropical geometry for agents — where + becomes max and × becomes +.
//!
//! Tropical geometry replaces usual arithmetic with the **max-plus semiring**:
//! - Tropical addition = max
//! - Tropical multiplication = +
//!
//! This simplifies polynomial systems and reveals combinatorial structure.
//! For agents, tropical optimization = optimal decisions via max-plus algebra.

pub mod semiring;
pub mod polynomial;
pub mod newton;
pub mod curves;
pub mod linear_algebra;
pub mod convexity;
pub mod optimization;
pub mod tropicalization;
pub mod intersection;
pub mod scheduling;

pub use semiring::Tropical;
pub use polynomial::TropicalPolynomial;
pub use linear_algebra::TropicalMatrix;
pub use convexity::TropicalPolytope;
pub use optimization::TropicalOptimizer;
pub use scheduling::AgentScheduler;
