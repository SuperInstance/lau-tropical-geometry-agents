# lau-tropical-geometry-agents

Tropical geometry for agents — where + becomes max and × becomes +.

Tropical geometry replaces the usual arithmetic with the **max-plus semiring**:
- **Tropical addition** = `max`
- **Tropical multiplication** = `+`

This simplifies polynomial systems and reveals combinatorial structure. For agents, tropical optimization = optimal decisions via max-plus algebra.

## Features

- **Tropical semiring**: `(ℝ ∪ {-∞}, max, +)` — tropical addition = max, tropical multiplication = +
- **Tropical polynomials**: max of affine functions = piecewise linear
- **Newton polytope**: tropical hypersurfaces correspond to polyhedral subdivisions
- **Tropical curves**: piecewise-linear skeletons of algebraic curves
- **Tropical linear algebra**: tropical eigenvalues, tropical determinant
- **Tropical convexity**: tropical polytopes, tropical halfspaces
- **Tropical optimization**: solve min-max problems via tropical algebra
- **Tropicalization**: classical → tropical as a deformation (t → 0 limit of log_t)
- **Tropical intersection theory**: stable intersections, tropical Bézout
- **Agent scheduling**: optimal scheduling and resource allocation via tropical optimization

## Usage

```rust
use lau_tropical_geometry_agents::{Tropical, TropicalPolynomial, TropicalMatrix, AgentScheduler};

// Tropical arithmetic
let a = Tropical::new(3.0);
let b = Tropical::new(5.0);
assert_eq!(a + b, Tropical::new(5.0)); // max(3, 5)
assert_eq!(a * b, Tropical::new(8.0)); // 3 + 5

// Tropical polynomial evaluation
let p = TropicalPolynomial::univariate(vec![(0.0, 0), (0.0, 1), (0.0, 2)]);
let val = p.evaluate(&[Tropical::new(3.0)]); // max(0, 3, 6) = 6

// Tropical matrix operations
let m = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);
let det = m.tropical_determinant(); // max(1+4, 2+3) = 5
let eig = m.tropical_eigenvalue(); // max cycle mean

// Agent scheduling
let scheduler = AgentScheduler::new(agents, tasks);
let assignments = scheduler.schedule();
```

## License

MIT
