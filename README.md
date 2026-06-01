# lau-tropical-geometry-agents

**Tropical geometry for agents — where + becomes max and × becomes +.**

A Rust library implementing tropical (max-plus) algebra: the tropical semiring, tropical polynomials and matrices, Newton polytopes, tropical curves, convexity, optimization, tropicalization, intersection theory, and agent scheduling. 136 tests, all passing.

[![136 tests passing](https://img.shields.io/badge/tests-136%20passing-brightgreen)]()

---

## Table of Contents

- [What This Does](#what-this-does)
- [Key Idea](#key-idea)
- [Install](#install)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
  - [Tropical (Semiring)](#tropical-semiring)
  - [MinPlus (Dual Semiring)](#minplus-dual-semiring)
  - [TropicalPolynomial](#tropicalpolynomial)
  - [TropicalMatrix](#tropicalmatrix)
  - [NewtonPolytope](#newtonpolytope)
  - [TropicalCurve](#tropicalcurve)
  - [TropicalPolytope](#tropicalpolytope)
  - [TropicalHalfspace](#tropicalhalfspace)
  - [TropicalOptimizer](#tropicaloptimizer)
  - [Tropicalization](#tropicalization)
  - [Intersection Theory](#intersection-theory-functions)
  - [AgentScheduler](#agentscheduler)
- [How It Works](#how-it-works)
- [The Math](#the-math)
  - [The Max-Plus Semiring](#the-max-plus-semiring)
  - [Tropical Polynomials as Piecewise-Linear Functions](#tropical-polynomials-as-piecewise-linear-functions)
  - [Tropical Linear Algebra](#tropical-linear-algebra)
  - [Newton Polytopes and Duality](#newton-polytopes-and-duality)
  - [Tropical Curves and the Balancing Condition](#tropical-curves-and-the-balancing-condition)
  - [Tropical Convexity](#tropical-convexity)
  - [Tropicalization: Classical → Tropical](#tropicalization-classical--tropical)
  - [Tropical Bézout and Intersection Theory](#tropical-bézout-and-intersection-theory)
  - [Agent Scheduling as Tropical Optimization](#agent-scheduling-as-tropical-optimization)
- [License](#license)

---

## What This Does

This library implements **tropical geometry** — a piecewise-linear shadow of algebraic geometry where classical arithmetic is replaced by the max-plus semiring:

- Tropical addition: `a ⊕ b = max(a, b)`
- Tropical multiplication: `a ⊗ b = a + b`

Under this transformation, polynomials become piecewise-linear functions, algebraic varieties become polyhedral complexes, and deep algebraic theorems become combinatorial ones.

Modules:

- **semiring** — The `Tropical` type (max-plus) and `MinPlus` type (min-plus), with full operator overloading
- **polynomial** — Tropical polynomials as max of affine functions, evaluation, dominant monomials, corner locus
- **linear_algebra** — Tropical matrices, tropical determinant, eigenvalue (via Karp's algorithm), eigenvectors, Kleene star
- **newton** — Newton polytopes, lattice points, polyhedral subdivisions, tropical hypersurfaces, Minkowski sum
- **curves** — Tropical curves as weighted balanced graphs, construction from polynomials, genus, balancing
- **convexity** — Tropical halfspaces, tropical polytopes, tropical convex hull, Hilbert projective metric
- **optimization** — Min-max optimization, tropical LP, Chebyshev approximation, tropical games, bottleneck assignment
- **tropicalization** — Classical → tropical via valuation map, initial forms, Gröbner fan directions, detropicalization
- **intersection** — Stable intersection, tropical Bézout theorem, intersection multiplicity, mixed volume
- **scheduling** — Agent-task scheduling via tropical optimization, critical path, resource allocation, bottleneck analysis

---

## Key Idea

Tropical geometry replaces the smooth curves and surfaces of classical algebraic geometry with **piecewise-linear skeletons**. The key insight: as you take the limit `t → 0` of `log_t(V)` where `V` is an algebraic variety, you get a polyhedral complex that retains combinatorial information about `V`.

This library makes that transformation computational:

1. **Tropical polynomials** are `max(a₀, a₁+x, a₂+2x, ...)` — piecewise-linear, concave
2. **Tropical varieties** are corner loci where these functions are non-differentiable
3. **Tropical linear algebra** replaces `+` with `max` and `×` with `+` in matrix operations
4. **Tropical eigenvalues** = max cycle mean of a weighted graph (Karp's algorithm)
5. **Tropical optimization** solves min-max problems that arise naturally in scheduling and resource allocation
6. **Tropical Bézout** counts intersections: two curves of degree `d₁` and `d₂` meet in `d₁·d₂` points (counted with multiplicity)

---

## Install

```toml
[dependencies]
lau-tropical-geometry-agents = "0.1.0"
```

Or:

```bash
cargo add lau-tropical-geometry-agents
```

Requires Rust 2021+. Dependencies: `nalgebra` (with serde), `serde`.

---

## Quick Start

### Tropical arithmetic

```rust
use lau_tropical_geometry_agents::Tropical;

let a = Tropical::new(3.0);
let b = Tropical::new(5.0);

// Tropical addition = max
assert_eq!(a + b, Tropical::new(5.0));

// Tropical multiplication = ordinary addition
assert_eq!(a * b, Tropical::new(8.0));

// Idempotent: a ⊕ a = a
assert_eq!(a + a, a);

// Power: a^3 = a + a + a = 9
assert_eq!(a.tropical_pow(3), Tropical::new(9.0));
```

### Tropical polynomials

```rust
use lau_tropical_geometry_agents::TropicalPolynomial;
use lau_tropical_geometry_agents::Tropical;

// max(2, x) — univariate, corner at x=2
let p = TropicalPolynomial::univariate(vec![(2.0, 0), (0.0, 1)]);
assert_eq!(p.evaluate(&[Tropical::new(5.0)]), Tropical::new(5.0)); // max(2, 5) = 5
assert_eq!(p.evaluate(&[Tropical::new(1.0)]), Tropical::new(2.0)); // max(2, 1) = 2

// Corner points (non-differentiable locus)
let corners = p.corner_points_1d();
assert!(corners.iter().any(|&c| (c - 2.0).abs() < 1e-10));
```

### Tropical matrix algebra

```rust
use lau_tropical_geometry_agents::TropicalMatrix;

let a = TropicalMatrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]);

// Tropical determinant = max over all permutations of sum of selected entries
assert_eq!(a.tropical_determinant(), Tropical::new(5.0)); // max(1+4, 2+3)

// Tropical eigenvalue (max cycle mean via Karp's algorithm)
let eig = a.tropical_eigenvalue();
assert!(eig.unwrap() >= Tropical::new(4.0));
```

### Newton polytope

```rust
use lau_tropical_geometry_agents::TropicalPolynomial;
use lau_tropical_geometry_agents::newton::NewtonPolytope;

let p = TropicalPolynomial::univariate(vec![(0.0, 0), (1.0, 1), (2.0, 2)]);
let np = NewtonPolytope::from_polynomial(&p);
assert_eq!(np.num_vertices(), 3);
assert_eq!(np.dimension(), 1);
```

### Agent scheduling

```rust
use lau_tropical_geometry_agents::{AgentScheduler, Agent, Task};

let agents = vec![
    Agent::new("alice", vec![5.0, 3.0]),
    Agent::new("bob", vec![2.0, 7.0]),
];
let tasks = vec![
    Task::new("deploy", vec![4.0, 1.0], 3.0),
    Task::new("test", vec![1.0, 4.0], 2.0).with_dependencies(vec![0]),
];

let scheduler = AgentScheduler::new(agents, tasks);
let assignments = scheduler.schedule();
let makespan = scheduler.makespan(&assignments);
let (bottleneck, impact) = scheduler.bottleneck_task();
```

---

## API Reference

### Tropical (Semiring)

```rust
pub struct Tropical(pub f64);
```

Elements of the max-plus semiring `(ℝ ∪ {-∞}, max, +)`.

| Constant | Value | Meaning |
|----------|-------|---------|
| `Tropical::ZERO` | `-∞` | Additive identity (max identity) |
| `Tropical::ONE` | `0.0` | Multiplicative identity |
| `Tropical::NEG_INF` | `-∞` | Same as ZERO |

| Method / Op | Description |
|-------------|-------------|
| `new(val)` | Wrap an `f64` |
| `from_int(val)` | From `i64` |
| `a + b` | Tropical addition = `max(a, b)` |
| `a * b` | Tropical multiplication = `a + b` |
| `a - b` | Tropical subtraction = `a - b` |
| `a.tropical_div(b)` | `Some(a - b)` or `None` if `b = -∞` |
| `a.tropical_pow(n)` | `a * n` (scalar multiplication) |
| `a.tropical_neg()` | `-a` |
| `a.tropical_abs()` | `|a|` |
| `a.is_zero()` | `a == -∞` |
| `a.is_one()` | `a == 0.0` |
| `a.value()` | Inner `f64` |

Implements `Ord`, `Sum`, `Product`, `Display`, `Serialize`/`Deserialize`. Idempotent: `a + a = a`.

### MinPlus (Dual Semiring)

```rust
pub struct MinPlus(pub f64);
```

Elements of the min-plus semiring `(ℝ ∪ {+∞}, min, +)`.

| Method | Description |
|--------|-------------|
| `new(val)` | Wrap an `f64` |
| `a.min_add(b)` | `min(a, b)` |
| `a.min_mul(b)` | `a + b` |
| `MinPlus::POS_INF` | `+∞` (additive identity) |
| `MinPlus::ONE` | `0.0` |

### TropicalPolynomial

```rust
pub struct TropicalPolynomial {
    pub monomials: Vec<TropicalMonomial>,
    pub num_vars: usize,
}
```

A tropical polynomial = `max` of monomials. Each `TropicalMonomial` has a `coefficient: Tropical` and `degree: Vec<u32>`.

| Method | Description |
|--------|-------------|
| `new(monomials, num_vars)` | Create from monomials |
| `univariate(terms)` | From `Vec<(coeff, degree)>` |
| `bivariate(terms)` | From `Vec<(coeff, deg1, deg2)>` |
| `evaluate(point)` | `max` of monomial evaluations |
| `tropical_add(other)` | Combine all monomials |
| `tropical_mul(other)` | Convolve (add coeffs, add degrees) |
| `scalar_mul(s)` | Add scalar to all coefficients |
| `dominant_monomials(point)` | Indices achieving the max |
| `corner_points_1d()` | Non-differentiable points (univariate) |
| `degree()` | Maximum total degree |
| `simplify()` | Remove dominated monomials |

### TropicalMatrix

```rust
pub struct TropicalMatrix {
    pub data: DMatrix<Tropical>,
}
```

Tropical matrix in the max-plus semiring.

| Method | Description |
|--------|-------------|
| `from_rows(rows)` | From `Vec<Vec<f64>>` |
| `from_dmatrix(m)` | From `DMatrix<f64>` |
| `identity(n)` | `n×n` identity (diagonal = 0, off-diag = -∞) |
| `zeros(nrows, ncols)` | All entries = -∞ |
| `tropical_add(other)` | Entrywise max |
| `tropical_mul(other)` | `(AB)_{ij} = max_k(A_{ik} + B_{kj})` |
| `tropical_pow(k)` | `k`-th tropical power |
| `tropical_determinant()` | Max over all permutations of sum |
| `tropical_eigenvalue()` | Max cycle mean (Karp's algorithm) |
| `tropical_eigenvectors(λ)` | Solve `A ⊗ v = λ ⊗ v` |
| `kleene_star()` | `I ⊕ A ⊕ A² ⊕ ... ⊕ Aⁿ` |
| `tropical_trace()` | Max of diagonal |
| `mul_vector(v)` | Tropical matrix-vector product |
| `transpose()` | Transpose |
| `scalar_mul(s)` | Add scalar to all entries |

### NewtonPolytope

```rust
pub struct NewtonPolytope {
    pub vertices: Vec<LatticePoint>,
    pub dim: usize,
}
```

Convex hull of exponent vectors from a tropical polynomial.

| Method | Description |
|--------|-------------|
| `from_polynomial(poly)` | Build from tropical polynomial |
| `num_vertices()` | Vertex count |
| `dimension()` | Affine dimension via rank computation |
| `upper_hull_indices(direction)` | Vertices maximizing a direction |
| `contains(point)` | Bounding-box containment check |
| `volume()` | Length (1D), shoelace area (2D) |
| `edges()` | All vertex pairs |
| `minkowski_sum(other)` | Minkowski sum of polytopes |

Supporting type: `LatticePoint { coords: Vec<i64> }` with `dot`, `sub`, `add`.

### TropicalCurve

```rust
pub struct TropicalCurve {
    pub vertices: Vec<CurveVertex>,
    pub edges: Vec<CurveEdge>,
    pub ambient_dim: usize,
}
```

A weighted balanced 1-dimensional polyhedral complex (tropical curve).

| Method | Description |
|--------|-------------|
| `new(vertices, edges, ambient_dim)` | Construct |
| `from_bivariate_polynomial(poly, resolution)` | Build from bivariate polynomial via grid sampling |
| `is_balanced_at(idx)` | Check balancing at vertex |
| `is_balanced()` | Check balancing at all vertices |
| `genus()` | `max(0, E - V + 1)` |
| `incident_edges(idx)` | Edges touching a vertex |
| `vertex_degree(idx)` | Number of incident edges |
| `total_weight()` | Sum of edge weights |

### TropicalPolytope

```rust
pub struct TropicalPolytope {
    pub generators: Vec<Vec<Tropical>>,
    pub dim: usize,
}
```

Tropical convex hull of generator points.

| Method | Description |
|--------|-------------|
| `new(generators)` | From `Vec<Vec<f64>>` |
| `contains(point)` | Approximate membership via tropical combination check |
| `convex_hull_points()` | Current generators |
| `tropical_segment(a, b)` | Points on the tropical segment between `a` and `b` |
| `tropical_distance(a, b)` | Hilbert projective metric |
| `minkowski_sum(other)` | Minkowski sum |
| `circumscribed_ball()` | Center and radius of circumscribed ball |

### TropicalHalfspace

```rust
pub struct TropicalHalfspace {
    pub coefficients: Vec<Tropical>,
    pub rhs: Tropical,
    pub is_upper: bool,
}
```

| Method | Description |
|--------|-------------|
| `new(coefficients, rhs, is_upper)` | Construct |
| `contains(point)` | Check `max(aᵢ+xᵢ) ≥/≤ c` |

### TropicalOptimizer

```rust
pub struct TropicalOptimizer { pub dim: usize }
```

| Method | Description |
|--------|-------------|
| `solve_minmax(cost)` | Min-max assignment |
| `eigenvalue_optimization(matrix)` | Tropical eigenvalue |
| `tropical_lp(cost, A, b)` | Tropical linear program |
| `chebyshev_approx(A, target)` | Tropical Chebyshev approximation |
| `tropical_nearest_point(polytope, point)` | Nearest generator |
| `solve_tropical_game(payoff)` | Tropical saddle point (lower, upper value) |
| `bottleneck_assignment(cost)` | Minimize maximum cost assignment |

### Tropicalization

```rust
pub struct Tropicalization { pub base: f64, ... }
pub struct ClassicalPolynomial { pub terms: Vec<ClassicalTerm>, pub num_vars: usize }
```

| Method | Description |
|--------|-------------|
| `ClassicalPolynomial::evaluate(point)` | Classical polynomial evaluation |
| `ClassicalPolynomial::tropicalize(base)` | Classical → tropical via valuation map |
| `Tropicalization::new(base)` | Create tropicalization context |
| `tr.tropicalize(poly)` | Tropicalize and store |
| `tr.initial_form(poly, weight)` | Leading terms w.r.t. weight vector |
| `tr.grobner_directions(poly)` | Gröbner fan wall normals |
| `detropicalize(tropical, base)` | Reverse: tropical → classical |

### Intersection Theory Functions

| Function | Description |
|----------|-------------|
| `tropical_bezout_number(d1, d2)` | `d₁ · d₂` |
| `stable_intersection(p1, p2, tol)` | Find intersection points |
| `intersection_multiplicity(p1, p2, point)` | Multiplicity at a point |
| `verify_bezout(p1, p2)` | `(actual, expected)` intersection count |
| `mixed_volume(np1, np2)` | Bernstein's mixed volume |
| `check_transverse(c1, c2)` | Check transverse intersection of curves |

### AgentScheduler

```rust
pub struct AgentScheduler {
    pub agents: Vec<Agent>,
    pub tasks: Vec<Task>,
}
```

| Method | Description |
|--------|-------------|
| `new(agents, tasks)` | Construct |
| `cost_matrix()` | Agent-task tropical cost matrix |
| `dependency_matrix()` | Tropical dependency graph |
| `critical_path_length()` | Tropical eigenvalue of dependency matrix |
| `earliest_start_times()` | Via tropical matrix iteration |
| `schedule()` | Produce assignments |
| `makespan(assignments)` | Total duration |
| `resource_utilization(assignments)` | Fraction of capacity used |
| `allocate_resources(budget)` | Proportional allocation |
| `bottleneck_task()` | `(index, impact)` of critical task |

Supporting types: `Agent { id, capabilities }`, `Task { id, requirements, duration, dependencies, priority }`, `Assignment { task_idx, agent_idx, start_time, end_time }`.

---

## How It Works

1. **Define the semiring**: `Tropical` wraps `f64` with `+` = `max` and `*` = `+`. The zero element is `-∞` and the unit is `0.0`. This is idempotent (`a ⊕ a = a`).

2. **Build polynomials**: A tropical polynomial is `max(c₀, c₁+x, c₂+2x, ...)` — a piecewise-linear concave function. The **corner locus** (where it's non-differentiable) is the tropical variety.

3. **Do linear algebra**: Tropical matrix multiplication is `(AB)_{ij} = max_k(A_{ik} + B_{kj})`. The tropical determinant maximizes over permutations. Tropical eigenvalues (max cycle mean) are computed via **Karp's algorithm**.

4. **Construct Newton polytopes**: The exponent vectors of a polynomial form a lattice polytope whose geometry controls the tropical hypersurface via **duality**.

5. **Build tropical curves**: From bivariate polynomials, sample a grid to find vertices where 3+ monomials achieve the max simultaneously. Connect nearby vertices. Verify the **balancing condition** at each vertex.

6. **Optimize**: Tropical LP, Chebyshev approximation, and game theory all reduce to max-plus algebra. The tropical eigenvalue gives the critical path length in scheduling.

7. **Tropicalize**: Apply the valuation map `x ↦ -log_t(|x|)` to convert classical polynomials into tropical ones. The `t → 0` limit extracts the combinatorial skeleton.

8. **Intersect**: Stable intersection is the limit of perturbed intersections. **Tropical Bézout** says two curves of degree `d₁` and `d₂` meet in `d₁·d₂` points.

---

## The Math

### The Max-Plus Semiring

The **max-plus semiring** `(ℝ ∪ {-∞}, ⊕, ⊗)` has:
- Addition: `a ⊕ b = max(a, b)`
- Multiplication: `a ⊗ b = a + b`

Properties:
- `⊕` is **idempotent**: `a ⊕ a = a`
- `⊕` is **commutative** and **associative**
- `⊗` distributes over `⊕`: `a ⊗ (b ⊕ c) = (a ⊗ b) ⊕ (a ⊗ c)`
- Zero element: `-∞` (since `max(a, -∞) = a`)
- Unit element: `0` (since `a + 0 = a`)
- **No additive inverse**: you can't undo `max`

The dual **min-plus semiring** `(ℝ ∪ {+∞}, min, +)` has the same structure with `min` replacing `max`.

### Tropical Polynomials as Piecewise-Linear Functions

A univariate tropical polynomial:

```
f(x) = max(a₀, a₁ + x, a₂ + 2x, ..., aₙ + nx)
```

This is a **piecewise-linear concave function**. Each monomial `aᵢ + ix` is an affine function with slope `i`, and the tropical polynomial takes their maximum.

The **corner locus** (tropical variety) consists of points where two or more monomials achieve the maximum simultaneously — these are the non-differentiable points. For monomials `aᵢ + ix` and `aⱼ + jx`, they're equal at `x = (aᵢ - aⱼ) / (j - i)`.

For multivariate polynomials, `f(x₁,...,xₙ) = maxₖ(cₖ + dₖ₁x₁ + ... + dₖₙxₙ)` where each term is a hyperplane.

### Tropical Linear Algebra

Tropical matrix operations replace `+` with `max` and `×` with `+`:

- **Matrix multiplication**: `(AB)_{ij} = max_k(A_{ik} + B_{kj})` — this is the longest-path weight in a weighted graph.
- **Determinant**: `det_⊕(A) = max_σ Σᵢ Aᵢ,σᵢ)` — the maximum weight matching.
- **Eigenvalue**: The tropical eigenvalue is the **maximum cycle mean**: `λ = max_i max_k (Aᵏ)ᵢᵢ / k`, computed efficiently via **Karp's algorithm**.
- **Kleene star**: `A* = I ⊕ A ⊕ A² ⊕ ... ⊕ Aⁿ` gives the closure (all-path best weights).

### Newton Polytopes and Duality

The **Newton polytope** of a polynomial is the convex hull of its exponent vectors. Tropical geometry reveals a beautiful duality:

- The tropical hypersurface is **dual** to a polyhedral subdivision of the Newton polytope
- Each region where a monomial dominates corresponds to a vertex of the Newton polytope
- Each corner point corresponds to an edge of the Newton polytope
- The coefficients of the polynomial define a **lifting function** that induces the subdivision

### Tropical Curves and the Balancing Condition

A **tropical curve** is a weighted 1-dimensional polyhedral complex. At each vertex, the **balancing condition** must hold:

```
Σ (weight_e × direction_e) = 0
```

summed over all edges incident to the vertex, where `direction_e` is the primitive integer direction vector. This is the tropical analogue of holomorphicity.

The **genus** is computed topologically: `g = E - V + 1` (for connected curves).

### Tropical Convexity

Tropical convexity generalizes classical convexity:
- The **tropical convex hull** of points `v₁, ..., vₖ` is the set of all tropical linear combinations `maxᵢ(λᵢ + vᵢ)`
- **Tropical halfspaces** are defined by tropical linear inequalities
- The **Hilbert projective metric** `d(a, b) = max(a - b) - min(a - b)` measures tropical distance (invariant under tropical scaling)

### Tropicalization: Classical → Tropical

Tropicalization is the process of converting classical algebraic varieties into tropical ones via the **valuation map**:

```
x ↦ -log_t(|x|)    as t → 0
```

For a polynomial `f = Σ cᵢx^{αᵢ}`, tropicalization gives:

```
Trop(f) = maxᵢ(val(cᵢ) + αᵢ · point)
```

where `val(c) = -log_t(|c|)` is the valuation of the coefficient.

Key constructions:
- **Initial form** `in_w(f)`: the terms of `f` with maximum weighted degree under weight vector `w`
- **Gröbner fan**: the fan of weight vectors giving the same initial form
- **Detropicalization**: reverse the map, `coefficient = base^{-val}`

### Tropical Bézout and Intersection Theory

**Tropical Bézout theorem**: Two tropical curves of degrees `d₁` and `d₂` in the plane intersect in `d₁ · d₂` points, counted with multiplicity.

The **stable intersection** is the limit of intersections under infinitesimal perturbation — this is the correct notion of intersection in tropical geometry (generic intersections are transverse).

**Bernstein's theorem** connects the count to the **mixed volume** of Newton polytopes: the number of common solutions of a polynomial system equals the mixed volume of their Newton polytopes.

### Agent Scheduling as Tropical Optimization

Agent scheduling maps naturally to tropical algebra:

- **Dependency graph** → tropical matrix `D` where `D[i][j] = duration(task_i)` if task `j` depends on task `i`
- **Earliest start times** → tropical matrix powers: `ES = D* ⊗ d`
- **Critical path length** → tropical eigenvalue of `D`
- **Assignment** → tropical bottleneck problem: minimize the maximum cost
- **Makespan** → maximum end time across all assignments

---

## License

MIT
