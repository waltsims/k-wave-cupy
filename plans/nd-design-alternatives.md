# N-Dimensional Design Alternatives Analysis

## The Question

How should we generalize 1D operators to N-dimensions while maintaining clean, readable code?

**Context**: We need to handle gradient/divergence operators, velocity components, and split density for 1D, 2D, and 3D cases.

---

## Alternative 1: Lists of Operators (Proposed in Plan)

### Implementation

```python
# Build operators per dimension
op_grad = [op_grad_x, op_grad_y, op_grad_z][:ndim]
op_div = [op_div_x, op_div_y, op_div_z][:ndim]
u = [ux, uy, uz][:ndim]
rho_split = [rhox, rhoy, rhoz][:ndim]

# Time loop
for i in range(ndim):
    u[i] -= (dt / rho0_staggered[i]) * diff(p, op_grad[i])

for i in range(ndim):
    div_u_i = diff(u[i], op_div[i])
    rho_split[i] -= dt * rho0 * div_u_i * nonlinear_factor(sum(rho_split))

rho_total = sum(rho_split)
p = c0**2 * rho_total
```

### Pros
1. **Explicit dimensionality**: Clear that each dimension has its own operator
2. **Loop-based logic**: Easy to read, step through mentally
3. **Debuggable**: Can inspect `u[0]`, `u[1]` individually
4. **Flexible**: Easy to add per-dimension logic (e.g., different PML per axis)
5. **Memory explicit**: Clear when arrays are allocated

### Cons
1. **Python lists not NumPy**: Loses some array operations (can't use `np.stack` easily)
2. **More verbose**: Explicit loops vs vectorized operations
3. **Line count**: ~40 lines for time loop vs potential ~20 with vectorization

### Code Metrics
- **Cognitive Load**: Low (5/10) - straightforward loops
- **Readability**: High (9/10) - intent is clear
- **Maintainability**: High (8/10) - easy to modify per-dimension
- **Performance**: Medium (6/10) - Python loops slower than vectorized ops

---

## Alternative 2: NumPy Arrays with Axis Parameter

### Implementation

```python
# Stack operators into arrays along a new "dimension axis"
# op_grad.shape = (ndim, Nx, Ny, Nz) where trailing dims are 1 if not used
op_grad = np.stack([op_grad_x, op_grad_y, op_grad_z][:ndim], axis=0)
op_div = np.stack([op_div_x, op_div_y, op_div_z][:ndim], axis=0)

# u.shape = (ndim, Nx, Ny, Nz)
u = np.zeros((ndim,) + grid_shape, dtype=float)
rho_split = np.zeros((ndim,) + grid_shape, dtype=float)

# Time loop - vectorized over dimensions
grad_p = np.array([diff(p, op_grad[i], axis=i) for i in range(ndim)])
u -= (dt / rho0_staggered) * grad_p  # Broadcasting

div_u = np.array([diff(u[i], op_div[i], axis=i) for i in range(ndim)])
rho_split -= dt * rho0 * div_u * nonlinear_factor(rho_split.sum(axis=0))

p = c0**2 * rho_split.sum(axis=0)
```

### Pros
1. **NumPy native**: Uses array operations throughout
2. **Potential performance**: Broadcasting can be faster than loops
3. **Compact**: Fewer lines of code
4. **Vectorized operations**: Can use `np.sum`, `np.stack` naturally

### Cons
1. **Cognitive overhead**: What does `axis=0` mean here? (dimension index, not spatial)
2. **Hidden dimensions**: Shape `(3, 64, 64, 1)` vs `(3, 64, 64)` - confusing
3. **Stacking overhead**: Need to stack/unstack for FFT operations
4. **Mixed indexing**: `u[i]` vs `rho_split.sum(axis=0)` - which axis is which?
5. **Harder to debug**: Can't easily print "just the x-component"
6. **Array creation overhead**: List comprehension still needed for `diff` calls

### Code Metrics
- **Cognitive Load**: High (8/10) - need to track multiple axis meanings
- **Readability**: Medium (5/10) - requires mental translation
- **Maintainability**: Low (4/10) - easy to confuse axes
- **Performance**: High (8/10) - vectorized operations

### Example Confusion

```python
# What does this mean?
grad_p = np.array([diff(p, op_grad[i], axis=i) for i in range(ndim)])
#                                          ^^^^^^ - spatial axis
#                             ^^^^^^^^^ - operator dimension axis
#                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ - list comp!
```

The axis parameter has two different meanings:
- `i` in `op_grad[i]`: which dimension's operator (conceptual)
- `axis=i` in `diff()`: which spatial axis to apply FFT (NumPy)

---

## Alternative 3: Dictionary-Based (Named Dimensions)

### Implementation

```python
# Operators as dictionaries
dim_names = ['x', 'y', 'z'][:ndim]
op_grad = {name: build_op(axis, ...) for axis, name in enumerate(dim_names)}
op_div = {name: build_op(axis, ...) for axis, name in enumerate(dim_names)}

u = {name: np.zeros(grid_shape) for name in dim_names}
rho_split = {name: np.zeros(grid_shape) for name in dim_names}

# Time loop
for dim in dim_names:
    u[dim] -= (dt / rho0_staggered[dim]) * diff(p, op_grad[dim])

for dim in dim_names:
    div_u_dim = diff(u[dim], op_div[dim])
    rho_split[dim] -= dt * rho0 * div_u_dim * nonlinear_factor(sum(rho_split.values()))

p = c0**2 * sum(rho_split.values())
```

### Pros
1. **Self-documenting**: `u['x']` is clearer than `u[0]`
2. **Semantic clarity**: "x dimension" vs "dimension 0"
3. **Robust to reordering**: If we change dimension order, names stay correct
4. **Pythonic**: Dictionaries are idiomatic in Python

### Cons
1. **Overkill for simple case**: More bookkeeping for unclear benefit
2. **Performance hit**: Dictionary lookups slower than list indexing
3. **Harder to vectorize**: Can't easily convert to NumPy operations
4. **More typing**: `rho_split['x']` vs `rho_split[0]`
5. **Ordering unclear**: `sum(rho_split.values())` - what order?

### Code Metrics
- **Cognitive Load**: Medium (6/10) - names help but structure is heavier
- **Readability**: High (8/10) - names are clear
- **Maintainability**: Medium (6/10) - more code to maintain
- **Performance**: Low (4/10) - dictionary overhead

---

## Alternative 4: Object-Oriented (Dimension Objects)

### Implementation

```python
class Dimension:
    def __init__(self, axis, N, dx, dt, c_ref):
        self.axis = axis
        self.N = N
        self.dx = dx
        self.op_grad = self._build_grad_operator()
        self.op_div = self._build_div_operator()
        self.u = np.zeros(grid_shape)
        self.rho = np.zeros(grid_shape)
        self.rho0_staggered = None  # Set later

    def update_velocity(self, p, dt, diff):
        grad_p = diff(p, self.op_grad)
        self.u -= (dt / self.rho0_staggered) * grad_p

    def update_density(self, dt, rho0, nonlinear_factor, diff):
        div_u = diff(self.u, self.op_div)
        self.rho -= dt * rho0 * div_u * nonlinear_factor

# Setup
dimensions = [Dimension(i, dims[i], spacing[i], dt, c_ref) for i in range(ndim)]

# Time loop
for dim in dimensions:
    dim.update_velocity(p, dt, diff)

for dim in dimensions:
    dim.update_density(dt, rho0, nonlinear_factor, diff)

rho_total = sum(dim.rho for dim in dimensions)
p = c0**2 * rho_total
```

### Pros
1. **Encapsulation**: Each dimension manages its own state
2. **Testable**: Can unit test `Dimension` class
3. **Extensible**: Easy to add dimension-specific behavior
4. **Clear ownership**: Each dimension "owns" its velocity and density

### Cons
1. **Over-engineering**: Too much structure for simple numerical code
2. **Performance overhead**: Object attribute access slower than array indexing
3. **Indirection**: Have to look in class definition to understand loop
4. **Not idiomatic for numerical code**: NumPy/SciPy codebases don't use OOP for numerical loops
5. **Hides the math**: Physics equations obscured by method calls
6. **Line count explosion**: ~100 lines for class definition

### Code Metrics
- **Cognitive Load**: High (9/10) - need to understand class structure
- **Readability**: Low (3/10) - math is hidden in methods
- **Maintainability**: Medium (7/10) - well-structured but over-complex
- **Performance**: Low (3/10) - object overhead in tight loop

### Philosophy Violation

This violates the "academic code golf" principle: **the code should be the documentation**. With OOP, the algorithm is scattered across class methods.

---

## Alternative 5: Hybrid - NumPy Arrays with Explicit Indexing

### Implementation

```python
# Store as NumPy arrays but index explicitly
u = np.zeros((ndim,) + grid_shape, dtype=float)
rho_split = np.zeros((ndim,) + grid_shape, dtype=float)
op_grad = [build_op(i, ...) for i in range(ndim)]
op_div = [build_op(i, ...) for i in range(ndim)]

# Time loop - explicit indexing, clear what's happening
for i in range(ndim):
    grad_p = diff(p, op_grad[i])
    u[i] -= (dt / rho0_staggered[i]) * grad_p

for i in range(ndim):
    div_u_i = diff(u[i], op_div[i])
    rho_split[i] -= dt * rho0 * div_u_i * nonlinear_factor(rho_split.sum(axis=0))

rho_total = rho_split.sum(axis=0)
p = c0**2 * rho_total
```

### Pros
1. **Best of both**: NumPy storage + explicit logic
2. **Performance**: Fast NumPy operations where it matters (`sum`, broadcasting)
3. **Readable loops**: Clear iteration over dimensions
4. **Easy aggregation**: `rho_split.sum(axis=0)` is clean
5. **Debuggable**: Can inspect `u[0]` like a list, but it's NumPy

### Cons
1. **Mixed paradigm**: Lists and arrays together might confuse
2. **Axis=0 overload**: Still have "dimension axis" concept
3. **Slightly less explicit**: `u[i]` vs `ux` - requires mental index-to-name mapping

### Code Metrics
- **Cognitive Load**: Medium (5/10) - one axis meaning to track
- **Readability**: High (8/10) - loops are clear, aggregation is clean
- **Maintainability**: High (8/10) - balance of structure and flexibility
- **Performance**: High (8/10) - NumPy operations where beneficial

---

## Alternative 6: Functional Style (Generator-Based)

### Implementation

```python
# Operators and initial state
def velocity_components(ndim, grid_shape):
    return (np.zeros(grid_shape) for _ in range(ndim))

def update_velocity(u_components, p, op_grad, rho0_staggered, dt, diff):
    """Generator that yields updated velocity components."""
    for u_i, op_i, rho0_i in zip(u_components, op_grad, rho0_staggered):
        grad_p = diff(p, op_i)
        yield u_i - (dt / rho0_i) * grad_p

# Time loop
u = list(update_velocity(u, p, op_grad, rho0_staggered, dt, diff))
rho_split = list(update_density(u, op_div, rho0, dt, diff))
p = c0**2 * sum(rho_split)
```

### Pros
1. **Functional purity**: No mutation (kind of - still mutating arrays)
2. **Composable**: Functions can be chained
3. **Memory efficient**: Generators avoid intermediate storage

### Cons
1. **Unnatural for numerical code**: Physics simulations are inherently stateful
2. **Harder to debug**: Can't easily inspect intermediate state
3. **Performance penalty**: Generator overhead + list() calls
4. **Readability loss**: What does `update_velocity()` return? Need to check signature
5. **Unnecessary abstraction**: Doesn't simplify the underlying math

### Code Metrics
- **Cognitive Load**: Very High (9/10) - need to understand generators
- **Readability**: Low (4/10) - obscures the physics
- **Maintainability**: Low (3/10) - hard to modify without breaking chain
- **Performance**: Medium (5/10) - generator overhead

---

## Comparative Analysis

### Summary Table

| Alternative | Lines | Readability | Performance | Maintainability | Cognitive Load | Philosophy Fit |
|-------------|-------|-------------|-------------|-----------------|----------------|----------------|
| 1. Lists    | 40    | 9/10        | 6/10        | 8/10            | 5/10           | ★★★★★          |
| 2. Arrays   | 25    | 5/10        | 8/10        | 4/10            | 8/10           | ★★☆☆☆          |
| 3. Dicts    | 45    | 8/10        | 4/10        | 6/10            | 6/10           | ★★★☆☆          |
| 4. OOP      | 100   | 3/10        | 3/10        | 7/10            | 9/10           | ★☆☆☆☆          |
| 5. Hybrid   | 35    | 8/10        | 8/10        | 8/10            | 5/10           | ★★★★☆          |
| 6. Functional| 50   | 4/10        | 5/10        | 3/10            | 9/10           | ★★☆☆☆          |

### "Academic Code Golf" Philosophy Alignment

**Core Principle**: The code is the documentation. A physicist should be able to read the implementation and immediately recognize the mathematical algorithm.

**Ranked by Philosophy Fit**:

1. **Lists (Alt 1)**: ★★★★★
   - Looks like pseudocode from a textbook
   - Each line maps to a physics equation
   - No hidden behavior

2. **Hybrid (Alt 5)**: ★★★★☆
   - Almost as clear as lists
   - NumPy operations are standard in scientific computing
   - Minor: `axis=0` requires knowing convention

3. **Dicts (Alt 3)**: ★★★☆☆
   - Extra naming helps documentation
   - But adds bookkeeping overhead
   - Not standard in numerical algorithms

4. **Arrays (Alt 2)**: ★★☆☆☆
   - Requires understanding NumPy broadcasting deeply
   - Axis confusion is a real problem
   - Less obvious what's happening

5. **Functional (Alt 6)**: ★★☆☆☆
   - Generators obscure the state evolution
   - Not how physicists think about time-stepping

6. **OOP (Alt 4)**: ★☆☆☆☆
   - Algorithm is fragmented across methods
   - Need to read entire class to understand loop
   - Over-abstracted for a simple numerical kernel

---

## Detailed Comparison: Lists vs Hybrid (Top 2)

### Example: 2D Momentum Equation

**Lists (Alt 1)**:
```python
# Momentum equation: du_i/dt = -(1/rho0) * dp/dx_i
for i in range(ndim):
    grad_p = diff(p, op_grad[i])
    u[i] -= (dt / rho0_staggered[i]) * grad_p
    u[i] = source_u_op[i](t, u[i])
```

**Hybrid (Alt 5)**:
```python
# Momentum equation: du_i/dt = -(1/rho0) * dp/dx_i
for i in range(ndim):
    grad_p = diff(p, op_grad[i])
    u[i] -= (dt / rho0_staggered[i]) * grad_p
    u[i] = source_u_op[i](t, u[i])
```

**Analysis**: They're identical in the time loop! The difference is in setup:

**Lists Setup**:
```python
u = [np.zeros(grid_shape) for _ in range(ndim)]
rho_split = [np.zeros(grid_shape) for _ in range(ndim)]
```

**Hybrid Setup**:
```python
u = np.zeros((ndim,) + grid_shape)
rho_split = np.zeros((ndim,) + grid_shape)
```

### When Hybrid Wins: Aggregation

**Lists**:
```python
rho_total = sum(rho_split)  # Works, but what's the behavior?
# Actually: sum([array1, array2]) -> array1 + array2 (correct!)
```

**Hybrid**:
```python
rho_total = rho_split.sum(axis=0)  # Explicit about summing over dimension axis
```

**Winner**: Hybrid is clearer for aggregation

### When Lists Win: Indexing Clarity

**Lists**:
```python
print(f"X velocity shape: {u[0].shape}")  # Obvious
u_components = u  # Can rename for clarity
```

**Hybrid**:
```python
print(f"X velocity shape: {u[0].shape}")  # Same
# But: u.shape is (3, 64, 64) - need to remember axis=0 is dimension
```

**Winner**: Lists slightly clearer (no axis=0 mental model)

### When Hybrid Wins: Initialization

**Lists**:
```python
# Have to initialize each component separately
for i in range(ndim):
    u[i] += (dt / (2 * rho0_staggered[i])) * diff(p, op_grad[i])
```

**Hybrid**:
```python
# Can broadcast if shapes align
u += (dt / (2 * rho0_staggered)) * grad_p_all  # If we precompute all
# But: still need loop for diff() calls, so no real advantage
```

**Winner**: Tie (both need loops for dimension-specific operations)

---

## Recommendation

### Primary Recommendation: **Alternative 5 (Hybrid)**

**Rationale**:

1. **Performance**: Gets NumPy benefits (fast sum, broadcasting) where it matters
2. **Readability**: Time loop is identical to Lists (9/10 readable)
3. **Pythonic**: NumPy arrays are the standard for multi-dimensional data
4. **Clean aggregation**: `rho_split.sum(axis=0)` is clearer than `sum(rho_split)`
5. **Maintains philosophy**: Code still maps directly to equations

**Key Implementation Points**:
```python
# Setup - NumPy arrays
u = np.zeros((ndim,) + grid_shape, dtype=float)
rho_split = np.zeros((ndim,) + grid_shape, dtype=float)

# Operators - lists (not arrays) because they have different shapes per dim
op_grad = [build_op(i, ...) for i in range(ndim)]

# Time loop - explicit loops for clarity
for i in range(ndim):
    u[i] -= ...  # Clear: updating i-th component

# Aggregation - NumPy methods
rho_total = rho_split.sum(axis=0)  # Clear: summing over dimension axis
```

**Convention to adopt**:
- **axis=0** always means "dimension index" (x, y, z)
- **axes 1, 2, 3** are spatial (Nx, Ny, Nz)
- Document this clearly at top of file

### Secondary Recommendation: **Alternative 1 (Lists)**

Use if:
- Performance of NumPy operations is not critical
- Absolute maximum readability is priority
- Team is less familiar with NumPy axis conventions

**Trade-off**: Slightly more verbose aggregation (`sum(rho_split)` vs `rho_split.sum(axis=0)`)

---

## Implementation Decision Tree

```
Does your team understand NumPy axis conventions?
├─ Yes → Use Hybrid (Alt 5)
│   └─ Benefit: Performance + readability
│
└─ No → Use Lists (Alt 1)
    └─ Benefit: Maximum clarity, easier onboarding
```

**For this project**: Use **Hybrid (Alternative 5)**

**Justification**:
1. Target audience: Researchers familiar with NumPy
2. Performance matters (large 3D grids)
3. axis=0 convention is well-documented
4. Clean aggregation is worth the minor axis complexity

---

## Code Example: Final Recommendation

```python
def simulate(kgrid, medium, source, sensor, backend="auto"):
    """N-D k-Space Pseudospectral Wave Propagator."""
    xp = cp if cp and backend in ("auto", "gpu") else np

    # Parse inputs
    grid_shape, dims, spacing, Nt, dt, c0, rho0, mask = _parse_inputs(kgrid, medium, sensor, xp)
    ndim = len(dims)

    # Build operators (lists - different shapes per dimension)
    op_grad, op_div, source_kappa = _build_kspace_operators(dims, spacing, dt, c_ref, xp)
    diff = lambda f, op: _spectral_diff(f, op, ndim, xp)

    # Initialize fields (NumPy arrays - homogeneous shape)
    p = xp.zeros(grid_shape, dtype=float)
    u = xp.zeros((ndim,) + grid_shape, dtype=float)
    rho_split = xp.zeros((ndim,) + grid_shape, dtype=float)

    # Staggered density (list - computed per dimension)
    rho0_staggered = [_compute_staggered(rho0, axis=i, xp=xp) for i in range(ndim)]

    # Initialize velocity (leapfrog staggering)
    for i in range(ndim):
        u[i] += (dt / (2 * rho0_staggered[i])) * diff(p, op_grad[i])

    for t in range(Nt):
        # Momentum: du_i/dt = -(1/rho0) * dp/dx_i
        for i in range(ndim):
            u[i] -= (dt / rho0_staggered[i]) * diff(p, op_grad[i])
            u[i] = source_u_op[i](t, u[i])

        # Mass: drho_i/dt = -rho0 * du_i/dx_i
        div_u_total = xp.zeros(grid_shape, dtype=float)
        for i in range(ndim):
            div_u_i = diff(u[i], op_div[i])
            div_u_total += div_u_i
            rho_split[i] -= dt * rho0 * div_u_i * nonlinear_factor(rho_split.sum(axis=0))

        # Equation of state
        rho_total = rho_split.sum(axis=0)  # Clear: sum over dimensions
        p = c0**2 * (rho_total + absorption(div_u_total) - dispersion(rho_total) + nonlinearity(rho_total))

        sensor_data[:, t] = p[mask]

    return {"sensor_data": _to_cpu(sensor_data)}
```

**Line count**: ~35 lines for time loop (vs 40 for pure lists)
**Readability**: 8/10 (clear physics, one axis convention to learn)
**Performance**: 8/10 (NumPy operations where beneficial)

---

## Final Answer

**Use Alternative 5 (Hybrid)** with these guidelines:

1. **Store multi-dimensional data as NumPy arrays**: `u`, `rho_split`
2. **Use lists for dimension-specific objects**: `op_grad`, `rho0_staggered`
3. **Explicit loops in time stepping**: Don't try to vectorize over dimensions
4. **NumPy methods for aggregation**: `.sum(axis=0)` over dimensions
5. **Document axis convention**: axis=0 is dimension index, axes 1+ are spatial

This balances readability (the code looks like the equations) with performance (NumPy operations are fast) while maintaining the "academic code golf" philosophy.
