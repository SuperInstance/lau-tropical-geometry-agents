# lau-tropical-geometry-agents

**Tropical geometry for agents — where + becomes max and × becomes +**

A Rust library that implements tropical geometry over the max-plus semiring, providing tropical polynomials, Newton polytopes, tropical linear algebra, tropical curves, convexity, optimization, intersection theory, and agent scheduling. All built on the deceptively simple idea that replacing addition with `max` and multiplication with `+` reveals hidden combinatorial structure.

136 tests · 11 modules · ~3,900 LOC

---

## What This Does

Tropical geometry replaces classical arithmetic with **tropical arithmetic**:

| Classical | Tropical (max-plus) |
|-----------|-------------------|
| `a + b` | `max(a, b)` |
| `a × b` | `a + b` |
| Additive identity `0` | `-∞` |
| Multiplicative identity `1` | `0` |

This library implements the full machinery of tropical geometry in Rust:

- **Tropical semiring** — `Tropical<f64>` with max-plus arithmetic, plus the dual min-plus semiring
- **Tropical polynomials** — piecewise-linear functions as `max` of affine terms, with evaluation, multiplication, corner locus, and simplification
- **Newton polytopes** — convex hulls of exponent vectors, Minkowski sums, upper hulls, polyhedral subdivisions
- **Tropical curves** — 1D polyhedral complexes from bivariate polynomials, with balancing checks and genus computation
- **Tropical linear algebra** — matrices over the max-plus semiring with tropical determinant, eigenvalues (Karp's algorithm), eigenvectors, Kleene star
- **Tropical convexity** — tropical polytopes, halfspaces, convex hull, containment, extreme points
- **Tropical optimization** — solve min-max problems, linear programming over the tropical semiring
- **Tropicalization** — transform classical polynomials to tropical via the valuation map `x ↦ -log_t(x)`
- **Intersection theory** — stable intersections, tropical Bézout's theorem
- **Agent scheduling** — optimal task scheduling and resource allocation via tropical matrix powers and critical path analysis

---

## Key Idea

Tropical geometry is the "skeleton" of classical algebraic geometry. By replacing `+` with `max` and `×` with `+`, algebraic varieties become **piecewise-linear polyhedral complexes**. This transformation:

1. **Simplifies computation** — polynomial systems become combinatorial optimization
2. **Preserves structure** — the tropical skeleton encodes genus, degree, and intersection numbers
3. **Enables optimization** — tropical matrix algebra solves longest-path and scheduling problems directly

The library connects pure mathematics to practical agent systems: tropical eigenvalues give **critical path lengths** in task graphs, and tropical matrix powers compute **optimal schedules** in max-plus algebra.

---

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-tropical-geometry-agents = { git = "https://github.com/SuperInstance/lau-tropical-geometry-agents" }
```

Or from a local clone:

```toml
[dependencies]
lau-tropical-geometry-agents = { path = "../lau-tropical-geometry-agents" }
```

### Dependencies

- `serde` 1.x (with `derive`) — serialization
- `nalgebra` 0.33 — linear algebra primitives

---

## Quick Start

```rust
use lau_tropical_geometry_agents::{Tropical, TropicalPolynomial, TropicalMatrix, AgentScheduler};

// Tropical arithmetic
let a = Tropical::new(3.0);
let b = Tropical::new(5.0);
assert_eq!(a + b, Tropical::new(5.0)); // max(3, 5) = 5
assert_eq!(a * b, Tropical::new(8.0)); // 3 + 5 = 8

// Tropical polynomial: max(0, x, 2x) is piecewise-linear
let p = TropicalPolynomial::univariate(vec![
    (0.0, 0),  // constant term: 0
    (0.0, 1),  // linear term: x
    (0.0, 2),  // quadratic term: 2x
]);
let val = p.evaluate(&[Tropical::new(3.0)]);
assert_eq!(val, Tropical::new(6.0)); // max(0, 3, 6) = 6

// Corner points (where the max switches)
let corners = p.corner_points_1d();
// All at x=0 for this symmetric polynomial

// Tropical matrix and eigenvalue
let m = TropicalMatrix::from_rows(&[
    vec![1.0, 2.0],
    vec![3.0, 4.0],
]);
let det = m.tropical_determinant(); // max(1+4, 2+3) = 5
let eig = m.tropical_eigenvalue();   // max cycle mean via Karp's algorithm

// Agent scheduling
let mut scheduler = AgentScheduler::new(3);
scheduler.add_task(/* ... */);
let schedule = scheduler.solve();
```

---

## API Reference

### `semiring` — Tropical Semiring

| Type / Function | Description |
|---|---|
| `Tropical(f64)` | Max-plus tropical number |
| `MinPlus(f64)` | Min-plus tropical number (dual semiring) |
| `Tropical::ZERO` | Additive identity: `-∞` |
| `Tropical::ONE` | Multiplicative identity: `0.0` |
| `a + b` | Tropical addition = `max(a, b)` |
| `a * b` | Tropical multiplication = `a + b` |
| `a.tropical_pow(n)` | Tropical exponentiation = scalar multiply by `n` |
| `a.tropical_div(b)` | Tropical division = ordinary subtraction |
| `TropicalOps` trait | Generic tropical arithmetic trait |

### `polynomial` — Tropical Polynomials

| Type / Function | Description |
|---|---|
| `TropicalPolynomial` | A tropical polynomial (max of monomials) |
| `TropicalMonomial` | A single monomial: `coefficient + dot(degree, vars)` |
| `TropicalPolynomial::univariate(terms)` | Create from `(coeff, degree)` pairs |
| `TropicalPolynomial::bivariate(terms)` | Create from `(coeff, deg1, deg2)` triples |
| `poly.evaluate(&point)` | Evaluate at a point (max of monomials) |
| `poly.tropical_add(&other)` | Tropical polynomial addition |
| `poly.tropical_mul(&other)` | Tropical polynomial multiplication (convolution) |
| `poly.dominant_monomials(&point)` | Indices of monomials achieving the max |
| `poly.corner_points_1d()` | Corner locus for univariate polynomials |
| `poly.simplify()` | Remove dominated monomials |

### `newton` — Newton Polytopes

| Type / Function | Description |
|---|---|
| `NewtonPolytope` | Convex hull of exponent vectors |
| `LatticePoint` | Integer lattice point |
| `NewtonPolytope::from_polynomial(&poly)` | Extract from polynomial |
| `np.dimension()` | Intrinsic dimension via Gaussian elimination |
| `np.volume()` | Volume (1D: length, 2D: shoelace area) |
| `np.upper_hull_indices(&direction)` | Upper convex hull in given direction |
| `np.minkowski_sum(&other)` | Minkowski sum of two polytopes |
| `PolyhedralSubdivision` | Regular subdivision induced by coefficients |
| `TropicalHypersurface` | Corner locus of a tropical polynomial |

### `curves` — Tropical Curves

| Type / Function | Description |
|---|---|
| `TropicalCurve` | Weighted balanced 1D polyhedral complex |
| `CurveVertex` | A vertex with coordinates |
| `CurveEdge` | An edge with weight and direction |
| `TropicalCurve::from_bivariate_polynomial(&poly, resolution)` | Extract curve from bivariate polynomial |
| `curve.is_balanced_at(i)` | Check balancing condition at vertex |
| `curve.genus()` | Genus (number of independent cycles) |
| `curve.incident_edges(i)` | Edges incident to a vertex |

### `linear_algebra` — Tropical Matrices

| Type / Function | Description |
|---|---|
| `TropicalMatrix` | Matrix over the max-plus semiring |
| `TropicalMatrix::identity(n)` | Identity matrix (ones on diagonal) |
| `m.tropical_mul(&other)` | Tropical matrix multiplication |
| `m.tropical_determinant()` | Max over all permutations of diagonal sums |
| `m.tropical_eigenvalue()` | Max cycle mean (Karp's algorithm) |
| `m.tropical_eigenvectors(λ)` | Eigenvectors for given eigenvalue |
| `m.kleene_star()` | Kleene star: `I ⊕ A ⊕ A² ⊕ ... ⊕ Aⁿ` |
| `m.mul_vector(&v)` | Tropical matrix-vector product |

### `convexity` — Tropical Convexity

| Type / Function | Description |
|---|---|
| `TropicalPolytope` | Tropical convex hull of points |
| `TropicalHalfspace` | Halfspace defined by tropical linear inequality |
| `polytope.contains(&point)` | Membership test |
| `polytope.extreme_points()` | Extreme (vertex) points |
| `polytope.intersection(&other)` | Intersection of two polytopes |

### `optimization` — Tropical Optimization

| Type / Function | Description |
|---|---|
| `TropicalOptimizer` | Solver for tropical optimization problems |
| `opt.solve_linear_program()` | Tropical linear programming |
| `opt.solve_min_max()` | Min-max optimization |

### `tropicalization` — Classical → Tropical

| Type / Function | Description |
|---|---|
| `ClassicalTerm` | A term in a classical polynomial |
| `tropicalize(&terms, base)` | Apply valuation map: `x ↦ -log_base(x)` |
| `TropicalizationResult` | Result of tropicalization |

### `intersection` — Intersection Theory

| Type / Function | Description |
|---|---|
| `IntersectionPoint` | A point where tropical varieties meet |
| `stable_intersection(&curve1, &curve2)` | Compute stable intersection |
| `tropical_bezout(&p1, &p2)` | Tropical Bézout's theorem |

### `scheduling` — Agent Scheduling

| Type / Function | Description |
|---|---|
| `AgentScheduler` | Schedule tasks across agents |
| `Agent` | An agent in the scheduling system |
| `scheduler.add_task(...)` | Add a task with dependencies |
| `scheduler.solve()` | Compute optimal schedule via tropical algebra |
| `scheduler.critical_path_length()` | Longest path (tropical eigenvalue) |

---

## How It Works

### Architecture

The library is structured as layers of increasing abstraction:

```
semiring (Tropical number type)
    └── polynomial (Tropical polynomials)
        ├── newton (Newton polytopes)
        ├── curves (Tropical curves)
        └── tropicalization (classical → tropical)
    └── linear_algebra (Tropical matrices)
        ├── convexity (Tropical polytopes & halfspaces)
        └── optimization (Tropical LP & min-max)
            └── scheduling (Agent task scheduling)
intersection (cross-cutting: intersection theory)
```

### Tropical Arithmetic

The max-plus semiring `(ℝ ∪ {-∞}, max, +)` has these properties:

- **Idempotent addition**: `max(a, a) = a` — no cancellation
- **No additive inverse** — you can't "subtract" in the usual sense
- **Tropical power**: `a^n = n × a` (scalar multiplication)
- **Tropical polynomial evaluation**: `max(a₀, a₁+x, a₂+2x, ...)` — piecewise-linear

### Tropical Eigenvalues via Karp's Algorithm

The tropical eigenvalue of a matrix is the **maximum cycle mean** of the weighted directed graph it represents. Karp's algorithm computes this in O(n³):

1. Compute `v_k[j]` = maximum weight of any path of length `k` ending at node `j`
2. The eigenvalue is `max_j max_k (v_n[j] - v_k[j]) / (n - k)`

This eigenvalue equals the **critical path length** in scheduling problems.

### Tropical Curves and Balancing

A tropical curve is a 1D polyhedral complex where every vertex satisfies the **balancing condition**: the weighted sum of outgoing direction vectors equals zero. This is the tropical analog of smoothness for algebraic curves.

---

## The Math

### The Tropical Semiring

The **max-plus semiring** is `(ℝ ∪ {-∞}, ⊕, ⊗)` where:
- `a ⊕ b = max(a, b)` (tropical addition)
- `a ⊗ b = a + b` (tropical multiplication)

The **min-plus semiring** `(ℝ ∪ {+∞}, min, +)` is its dual, useful for shortest-path problems.

### Tropical Polynomials

A tropical polynomial in `n` variables is:

$$P(x_1, \ldots, x_n) = \max_{\alpha} (c_\alpha + \alpha_1 x_1 + \cdots + \alpha_n x_n)$$

This is a **piecewise-linear convex function**. The **tropical hypersurface** is the corner locus — the set where `P` is not differentiable (where two or more monomials simultaneously achieve the maximum).

### Newton Polytopes

The **Newton polytope** `New(P)` is the convex hull of the exponent vectors `{α}` in `ℝⁿ`. The tropical hypersurface induces a polyhedral subdivision of `New(P)` that is **dual** to it — vertices correspond to regions, edges to facets, faces to cells.

### Tropical Bézout's Theorem

Two tropical curves of degrees `d₁` and `d₂` in the tropical plane intersect in exactly `d₁ · d₂` points (counted with multiplicity), mirroring the classical Bézout's theorem.

### Tropical Eigenvalue Theory

For a square tropical matrix `A`, the **tropical eigenvalue** `λ` satisfies:

$$\max_j (A_{ij} + v_j) = \lambda + v_i \quad \text{for all } i$$

This equals the **maximum cycle mean** of the associated weighted digraph, computable via Karp's algorithm or the tropical characteristic polynomial.

### Tropicalization

Given a classical polynomial `f = Σ c_α x^α`, its **tropicalization** is obtained by the valuation map:

$$x_i \mapsto -\log_t(x_i) \quad \text{as } t \to 0$$

The tropicalization `Trop(f)` records the combinatorial skeleton of the variety `{f = 0}`.

---

## Running Tests

```bash
cargo test
```

136 tests covering all modules: semiring properties, polynomial evaluation, Newton polytope geometry, curve balancing, matrix operations, scheduling, and more.

---

## License

MIT
