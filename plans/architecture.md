# k-Wave Python Backend Architecture

## Philosophy: Physics as Data, Not Control Flow

The k-Wave MATLAB reference implementation (`kspaceFirstOrder1D.m`) contains ~1000 lines with extensive conditional logic, creating a combinatorial explosion of code paths:
- 2 (linear/nonlinear) × 3 (equation states) × 3 (source modes) = **18 code paths**

Our Python/CuPy implementation avoids this by treating different physics as **composable operators that sum to zero when inactive**, maintaining a single clean code path regardless of enabled features.

## Core Design Principle: Operator Composition

Instead of conditional branches:
```python
# ❌ Bad (MATLAB approach)
if equation_of_state == 'lossless':
    p = c0**2 * rho
elif equation_of_state == 'absorbing':
    p = c0**2 * (rho + tau*term1 - eta*term2)
elif equation_of_state == 'stokes':
    p = c0**2 * (rho + tau*term3)
```

Use additive composition:
```python
# ✅ Good (functional composition)
p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))
# where each operator returns 0 if that physics isn't present
```

## Target Architecture

### Time Loop (8 lines, zero conditionals)
```python
for t in range(Nt):
    sensor_data[:, t] = p[mask]

    # Momentum equation with velocity sources
    grad_p = diff(p, op_grad)
    u -= (dt / rho0_sgx) * (grad_p + source_u(t, u))

    # Mass conservation with nonlinearity
    duxdx = diff(u, op_div)
    rho += dt * (-rho0 * duxdx * nonlinear_factor(rho))

    # Equation of state (all physics terms compose additively)
    p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))

    # Pressure sources (additive or dirichlet)
    p = source_p(t, p)
```

### Operator Factory Pattern

All complexity is isolated in `_build_physics_ops()` during initialization:

```python
def _build_physics_ops(medium, source, k, kappa, xp, Nx, dt):
    """
    Build composable physics operators.
    Each returns 0 (or identity) if feature not present.
    """
    ops = {}

    # --- Absorption (fractional Laplacian) ---
    alpha_coeff = medium.get("alpha_coeff", 0)
    if alpha_coeff > 0:
        # Build k-space operators for absorption/dispersion
        absorb_tau, absorb_eta = _absorption_coeffs(alpha_coeff, alpha_power, c_ref)
        absorb_nabla1 = (xp.abs(k) ** alpha_power) * xp.exp(1j * np.pi * alpha_power / 2)
        absorb_nabla2 = (xp.abs(k) ** alpha_power)

        ops['absorption'] = lambda duxdx: absorb_tau * diff_k(duxdx, absorb_nabla1)
        ops['dispersion'] = lambda rho: absorb_eta * diff_k(rho, absorb_nabla2)
    else:
        ops['absorption'] = lambda duxdx: 0
        ops['dispersion'] = lambda rho: 0

    # --- Nonlinearity (BonA parameter) ---
    BonA = medium.get("BonA", 0)
    if BonA != 0:
        ops['nonlinearity'] = lambda rho: BonA * rho**2 / (2 * rho0)
        ops['nonlinear_factor'] = lambda rho: (2*rho + rho0) / rho0
    else:
        ops['nonlinearity'] = lambda rho: 0
        ops['nonlinear_factor'] = lambda rho: 1.0

    # --- Source Modes (p, u with additive/dirichlet/additive-no-correction) ---
    ops['source_p'] = _build_source_operator(source, 'p', k, dt, xp)
    ops['source_u'] = _build_source_operator(source, 'u', k, dt, xp)

    return ops
```

## Key Benefits

### Eliminates Combinatorial Explosion
- **MATLAB**: 18 distinct code paths in time loop
- **Python**: 1 code path, operators compose additively

### Preserves Readability
- Time loop reads like textbook physics equations
- All complexity isolated in setup phase
- Each operator has single responsibility

### Enables Incremental Development
- Add new physics: Create operator in `_build_physics_ops()`
- Time loop remains unchanged
- Each operator independently testable against MATLAB reference

### Performance
- **Zero-cost abstraction**: Returning 0 is cheaper than branch misprediction
- **GPU-friendly**: Kernels fuse operator chains automatically
- **Memory efficient**: Operators share k-space arrays via closures

## Code Size Target

With **all features** implemented (absorption, nonlinearity, sources, PML):
- **Time loop**: 10 lines (unchanged from current)
- **Setup**: 40 lines (grid, parameters, k-space operators)
- **Operator factory**: 100 lines (~20 lines per operator × 5 operators)
- **Utilities**: 30 lines (absorption coefficients, helpers)
- **Total**: ~180 lines (vs 1000+ lines in MATLAB reference)

## Migration Path

### Week 1: Refactor to Operator Pattern
- Extract current physics into operators
- **Goal**: Same features, cleaner structure
- **Validation**: Existing 2 passing tests remain passing

### Week 2: Add Absorption
- Implement `absorption()` and `dispersion()` operators
- **Tests**: `linear + lossy` and `linear + stokes` cases
- **Expected**: ~30 lines added to operator factory

### Week 3: Add Nonlinearity
- Implement `nonlinearity()` and `nonlinear_factor()` operators
- **Tests**: `nonlinear + lossless/lossy/stokes` cases
- **Expected**: ~20 lines added to operator factory

### Week 4: Add Source Modes
- Implement `source_p()` and `source_u()` with all 3 modes
- **Tests**: All source.p and source.u test cases (42 tests)
- **Expected**: ~40 lines added to operator factory

### Week 5: Add PML Boundaries
- Implement PML as multiplicative operator
- **Tests**: Validate against MATLAB PML behavior
- **Expected**: ~20 lines added to operator factory

## Critical Success Criteria

This architecture **succeeds** if:
- ✅ Time loop stays under 15 lines
- ✅ No `if` statements in time loop (except initialization)
- ✅ Each new feature adds <40 lines to operator factory
- ✅ MATLAB test parity maintained at <1e-14 error

This architecture **fails** if:
- ❌ You add `if feature_enabled:` checks in loop
- ❌ Operators become stateful or have side effects
- ❌ Can't explain each operator in one sentence
- ❌ Test parity degrades

## Operator Contract

Each physics operator must satisfy:

```python
def operator(state_variables) -> term:
    """
    Pure function: Same inputs always give same outputs.
    Returns: Contribution to equation (0 if feature inactive).
    No side effects: Doesn't modify state_variables.
    """
```

**Examples**:
- `absorption(duxdx) -> τ ∇^y (ρ₀ ∂u/∂x)` or 0
- `nonlinearity(rho) -> (B/A) ρ² / (2ρ₀)` or 0
- `source_p(t, p) -> p + source_term` or `p` (identity)

## Testing Strategy

Each operator is independently testable:

```python
def test_absorption_operator():
    """Verify absorption operator matches MATLAB 'absorbing' equation state."""
    # Setup: medium with only absorption enabled
    medium = {"alpha_coeff": 0.5, "alpha_power": 1.5, "c0": 1500}
    ops = _build_physics_ops(medium, {}, k, kappa, np, Nx, dt)

    # Verify: absorption/dispersion return non-zero, nonlinearity returns 0
    assert ops['absorption'](test_duxdx) != 0
    assert ops['dispersion'](test_rho) != 0
    assert ops['nonlinearity'](test_rho) == 0
```

## Comparison: MATLAB vs Python

### MATLAB Time Loop (Lines 700-791, 92 lines)
```matlab
% Nested conditionals for mass conservation
if ~flags.nonlinear
    rhox = pml_x .* ( pml_x .* rhox - dt .* rho0 .* duxdx );
else
    rhox = pml_x .* ( pml_x .* rhox - dt .* (2 .* rhox + rho0) .* duxdx );
end

% 3-way switch for pressure sources
if flags.source_p >= t_index
    switch source.p_mode
        case 'dirichlet'
            rhox(p_source_pos_index) = source.p(...);
        case 'additive'
            source_mat = castZeros([kgrid.Nx, 1]);
            source_mat(p_source_pos_index) = source.p(...);
            source_mat = real(ifft(source_kappa .* fft(source_mat)));
            rhox = rhox + source_mat;
        case 'additive-no-correction'
            rhox(p_source_pos_index) = rhox(...) + source.p(...);
    end
end

% 2×3 nested switches for equation of state
if ~flags.nonlinear
    switch equation_of_state
        case 'lossless'
            p = c0.^2 .* rhox;
        case 'absorbing'
            p = c0.^2 .* (rhox + absorb_tau .* ... - absorb_eta .* ...);
        case 'stokes'
            p = c0.^2 .* (rhox + absorb_tau .* rho0 .* duxdx);
    end
else
    % Duplicate 3-way switch with BonA term added
    switch equation_of_state
        case 'lossless' ...
        case 'absorbing' ...
        case 'stokes' ...
    end
end
```
**Cognitive Load**: High - must trace 18 possible execution paths

### Python Time Loop (8 lines)
```python
for t in range(Nt):
    sensor_data[:, t] = p[mask]
    grad_p = diff(p, op_grad)
    u -= (dt / rho0_sgx) * (grad_p + source_u(t, u))
    duxdx = diff(u, op_div)
    rho += dt * (-rho0 * duxdx * nonlinear_factor(rho))
    p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))
    p = source_p(t, p)
```
**Cognitive Load**: Low - equations match textbook, 1 execution path

## Why This Architecture Wins

### Technical Advantages
1. **Correctness**: Impossible to have "forgot absorption in nonlinear+Stokes case" bugs
2. **Performance**: No branch misprediction, GPU kernel fusion
3. **Maintainability**: New features don't touch existing code
4. **Testability**: Each operator validates independently

### Philosophical Alignment
- **"Code is documentation"**: Time loop IS the physics equations
- **Academic code golf**: Maximum clarity, minimum lines
- **Incremental validation**: Add one operator, test immediately

## Final Recommendation

**Adopt the Operator Factory Pattern immediately.**

The current code (80 lines, 2 features) is already 90% compatible with this architecture. The refactor requires:
- Extract physics into operators (5 lines changed in time loop)
- Add operator factory function (~40 lines)
- Maintain all existing tests passing

The payoff: Seamlessly add 10+ physics features while keeping the time loop pristine and readable forever.
