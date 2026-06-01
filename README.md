# lau-tropical-geometry-agents

**Tropical geometry for agents — where + becomes max and × becomes +.**

A Rust library implementing the **max-plus semiring** (ℝ ∪ {−∞}, max, +) and the machinery of tropical geometry: tropical polynomials, tropical linear algebra, Newton polytopes, tropical curves, convexity, intersection theory, optimization, and agent scheduling.

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

---

## What This Does

Tropical geometry replaces ordinary arithmetic with tropical arithmetic:

| Operation | Usual | Tropical |
|-----------|-------|----------|
| Addition  | a + b | max(a, b) |
| Multiplication | a × b | a + b |
| Zero (additive identity) | 0 | −∞ |
| One (multiplicative identity) | 1 | 0 |

This deceptively simple change transforms polynomial systems into **piecewise-linear** objects with rich combinatorial structure. Tropical polynomials become convex piecewise-linear functions; their zero sets are **tropical hypersurfaces** — polyhedral complexes that encode combinatorial type.

This crate uses that machinery to:

- Solve **tropical optimization** problems (shortest paths, scheduling, optimal decisions)
- Compute with **tropical polynomials** and their Newton polytopes
- Build and analyze **tropical curves** (skeletons of classical algebraic curves)
- Perform **tropical linear algebra** (determinants, eigenvalues, Kleene star)
- Schedule agents via max-plus matrix algebra

---

## Key Idea

In the tropical world, a polynomial like `max(a₁+x₁, a₂+x₂, …, aₙ+xₙ)` encodes a shortest-path or optimal-decision problem. The tropical semiring's idempotent addition (max(a, a) = a) means there is no "cancelation" — every term matters combinatorially. This makes tropical geometry a natural language for:

- **Optimization**: Tropical matrix multiplication is the Floyd–Warshall algorithm
- **Scheduling**: The tropical eigenvalue of a dependency matrix gives the critical path length
- **Agent systems**: Agent scheduling reduces to tropical matrix powers

---

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-tropical-geometry-agents = "0.1"
```

Requires Rust 2021 edition. Dependencies: `nalgebra` and `serde`.

---

## Quick Start

```rust
use lau_tropical_geometry_agents::{Tropical, TropicalPolynomial, TropicalMatrix, AgentScheduler};

// --- Tropical arithmetic ---
let a = Tropical::new(3.0);
let b = Tropical::new(5.0);
assert_eq!(a.tropical_add(b), Tropical::new(5.0));  // max(3, 5) = 5
assert_eq!(a.tropical_mul(b), Tropical::new(8.0));   // 3 + 5 = 8

// --- Tropical polynomials ---
let poly = TropicalPolynomial::bivariate(
    vec![(1.0, 1, 0), (2.0, 0, 1), (0.0, 0, 0)],  // max(x + 1, y + 2, 0)
);
let val = poly.evaluate(&[3.0, 1.0]); // max(4, 3, 0) = 4

// --- Tropical matrix (shortest paths) ---
let dist = TropicalMatrix::from_rows(&[
    &[0.0, 5.0, f64::NEG_INFINITY],
    &[f64::NEG_INFINITY, 0.0, 3.0],
    &[2.0, f64::NEG_INFINITY, 0.0],
]);
let closure = dist.kleene_star(); // All-pairs shortest paths

// --- Agent scheduling ---
let mut scheduler = AgentScheduler::new(4);
scheduler.set_duration(0, 1, 3.0); // Task 0 → Task 1 takes 3 units
scheduler.set_duration(1, 2, 2.0);
scheduler.set_duration(0, 2, 6.0);
let makespan = scheduler.makespan(); // Critical path length
```

---

## API Reference

### `Tropical` — The Max-Plus Number Type

```rust
pub struct Tropical(pub f64);
```

| Method | Description |
|--------|-------------|
| `Tropical::NEG_INF` | Additive identity (−∞) |
| `Tropical::ONE` | Multiplicative identity (0.0) |
| `tropical_add(self, other)` | max(self, other) |
| `tropical_mul(self, other)` | self + other |
| `tropical_sub(self, other)` | self − other |
| `tropical_div(self, other)` | self − other (same as sub) |
| `tropical_pow(self, n)` | Repeated tropical multiplication = n·self |
| `to_f64(self) → f64` | Unwrap the inner value |

Implements `Add`, `Mul`, `Display`, `Serialize`, `Deserialize`, `PartialOrd`.

### `TropicalPolynomial` — Tropical Polynomial Functions

```rust
pub struct TropicalPolynomial { /* monomials with coefficients and exponents */ }
```

| Method | Description |
|--------|-------------|
| `from_monomials(coeffs, exponents)` | Build from raw data |
| `bivariate(terms)` | Convenience for 2-variable polynomials |
| `evaluate(&[x₁, …, xₙ])` | Evaluate at a point |
| `degree() → u32` | Total degree |
| `num_variables() → usize` | Number of variables |
| `tropical_add(&self, other)` | Pointwise max |
| `tropical_mul(&self, other)` | Tropical convolution |
| `corner_locus() → Vec<Vec<f64>>` | Points where two monomials tie |

### `NewtonPolytope` — Newton Polytope of a Polynomial

| Method | Description |
|--------|-------------|
| `from_polynomial(&poly)` | Extract the Newton polytope |
| `upper_hull_indices(direction)` | Upper hull facets |
| `contains(point)` | Point-in-polytope test |
| `volume() → f64` | Volume (2D → area) |
| `edges()` | Edge list |
| `minkowski_sum(other)` | Minkowski sum of polytopes |

### `TropicalCurve` — Skeleton Graph of a Tropical Hypersurface

| Method | Description |
|--------|-------------|
| `from_bivariate_polynomial(poly, resolution)` | Extract curve from polynomial |
| `is_balanced_at(vertex)` | Check balancing condition at vertex |
| `is_balanced()` | Check global balancing |
| `genus() → usize` | First Betti number (number of loops) |

### `TropicalMatrix` — Max-Plus Matrix Algebra

| Method | Description |
|--------|-------------|
| `from_rows(&[rows])` | Construct from row data |
| `identity(n)` | Tropical identity (0 on diagonal, −∞ elsewhere) |
| `tropical_mul(&self, other)` | Max-plus matrix multiply |
| `tropical_pow(k)` | k-th tropical power |
| `tropical_determinant() → Tropical` | Max of permutation products |
| `tropical_eigenvalue() → Option<Tropical>` | Cycle mean (critical path weight) |
| `kleene_star() → TropicalMatrix` | Transitive closure (all-pairs shortest paths) |
| `tropical_trace() → Tropical` | max of diagonal entries |

### `TropicalPolytope` — Tropical Convexity

| Method | Description |
|--------|-------------|
| `new(generators)` | Create from generator points |
| `contains(point)` | Membership test |
| `tropical_distance(a, b)` | Tropical (Chebyshev-like) distance |
| `minkowski_sum(other)` | Minkowski sum |
| `circumscribed_ball()` | Center and radius of enclosing ball |

### `TropicalOptimizer` — Tropical Optimization

| Method | Description |
|--------|-------------|
| `shortest_path(matrix, from, to)` | Solve via tropical matrix algebra |
| `all_pairs(&matrix)` | Kleene star (all-pairs shortest paths) |
| `assignment_problem(cost_matrix)` | Optimal assignment via tropical determinant |

### `AgentScheduler` — Agent/Task Scheduling via Tropical Algebra

| Method | Description |
|--------|-------------|
| `new(n_tasks)` | Create scheduler for n tasks |
| `set_duration(i, j, d)` | Set duration for task i → j |
| `makespan() → f64` | Critical path length (tropical eigenvalue) |
| `schedule() → Vec<f64>` | Start times for each task |

---

## How It Works

### Architecture

```
Tropical (number type)
  └─ TropicalPolynomial (piecewise-linear functions)
       ├─ NewtonPolytope (combinatorial type)
       ├─ TropicalCurve (skeleton graph)
       └─ TropicalHypersurface (zero set)
  └─ TropicalMatrix (max-plus linear algebra)
       ├─ TropicalOptimizer (shortest paths, assignment)
       └─ AgentScheduler (critical-path scheduling)
```

### Module Map

| Module | Contents |
|--------|----------|
| `semiring` | `Tropical` number type with max-plus arithmetic |
| `polynomial` | `TropicalPolynomial` — evaluation, addition, multiplication |
| `newton` | `NewtonPolytope`, `PolyhedralSubdivision`, `TropicalHypersurface` |
| `curves` | `TropicalCurve` — vertices, edges, balancing, genus |
| `linear_algebra` | `TropicalMatrix` — determinant, eigenvalue, Kleene star |
| `convexity` | `TropicalPolytope`, `TropicalHalfspace`, distance metrics |
| `intersection` | Tropical Bézout, stable intersection, mixed volume |
| `tropicalization` | Degeneration from classical to tropical curves |
| `optimization` | `TropicalOptimizer` — shortest paths, assignment |
| `scheduling` | `AgentScheduler` — critical-path task scheduling |

---

## The Math

### The Max-Plus Semiring

The **tropical semiring** is (ℝ ∪ {−∞}, ⊕, ⊗) where:

- a ⊕ b = max(a, b)
- a ⊗ b = a + b
- The additive identity is −∞ (denoted ⊥ or ε)
- The multiplicative identity is 0

This is also called the **max-plus algebra**. It is **idempotent**: a ⊕ a = a, meaning there is no additive inverse. This idempotence is what makes tropical geometry combinatorial — addition becomes a selection operation.

### Tropical Polynomials

A tropical polynomial in n variables:

> ⊕ᵢ (cᵢ ⊗ x₁^αᵢ¹ ⊗ … ⊗ xₙ^αᵢⁿ) = maxᵢ (cᵢ + αᵢ¹·x₁ + … + αᵢⁿ·xₙ)

This is a **convex piecewise-linear** function. Its "zero set" (where the maximum is achieved by at least two monomials) is the **tropical hypersurface** — a polyhedral complex dual to the Newton polytope's subdivision.

### Newton Polytopes

The **Newton polytope** of a tropical polynomial is the convex hull of its exponent vectors. The upper convex hull, projected down, gives the tropical hypersurface's combinatorial structure. Each facet of the Newton polytope corresponds to a region where a particular monomial dominates.

### Tropical Curves

For a bivariate tropical polynomial, the tropical curve is a **metric graph** embedded in ℝ². It satisfies the **balancing condition**: at every vertex, the sum of primitive edge directions (weighted by edge weights) is zero. This is the tropical analog of the classical curve's smoothness.

### Tropical Linear Algebra

A **tropical matrix** encodes a weighted directed graph. Key results:

- **Tropical determinant** = max-weight perfect matching (assignment problem)
- **Tropical eigenvalue** = maximum cycle mean = critical path length
- **Kleene star** A* = I ⊕ A ⊕ A² ⊕ … = transitive closure (all-pairs shortest paths)
- **Tropical matrix power** Aᵏ = k-step optimal paths

### Tropical Bézout's Theorem

Two tropical curves of degrees d₁ and d₂ intersect in exactly d₁ · d₂ points (counted with multiplicity), mirroring the classical Bézout theorem. This is verified by computing the **mixed volume** of their Newton polytopes.

### Agent Scheduling

Given a task dependency graph with durations d(i,j), the **makespan** (total completion time) equals the tropical eigenvalue of the duration matrix — the maximum cycle mean. Individual start times come from the tropical eigenvector, computed via the Kleene star.

---

## Testing

```bash
cargo test
```

The test suite includes **136 tests** covering:

- Tropical arithmetic properties (idempotence, associativity, distributivity)
- Polynomial evaluation, addition, and multiplication
- Newton polytope geometry (volume, edges, containment)
- Tropical curve balancing and genus computation
- Tropical matrix operations (determinant, eigenvalue, Kleene star)
- Convexity (polytope containment, distances, Minkowski sums)
- Intersection theory (Bézout verification, stable intersection)
- Optimization (shortest paths, assignment problems)
- Agent scheduling (makespan, start times)

---

## License

MIT
