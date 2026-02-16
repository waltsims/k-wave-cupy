# N-Dimensional (2D/3D) Implementation Plan for k-Wave Python Backend

## Executive Summary

**Goal**: Extend the minimalistic Python/CuPy backend from 1D to N-dimensional (2D/3D) support while maintaining the "academic code golf" philosophy and achieving <1e-14 parity with MATLAB reference implementations.

**Current State**:
- 1D implementation complete: ~65 lines Python, ~40 lines MATLAB
- 84/84 comprehensive 1D tests passing with <1e-15 parity
- Physics features implemented: absorption, dispersion, nonlinearity, heterogeneous media
- Missing: PML boundaries, 2D/3D support

**Key Insight from Analysis**: The 2D/3D MATLAB implementations show that the core algorithm is identical to 1D, with the key differences being:
1. Multiple gradient operators (one per dimension: `ddx_k`, `ddy_k`, `ddz_k`)
2. Multiple velocity components (ux, uy, uz)
3. Split density components (rhox, rhoy, rhoz) for each dimension
4. Sum operations for divergence and equation of state

---

## Part 1: Code Analysis and Architecture Review

### Current 1D Implementation Structure

**Python (`kWavePy.py`) - Core Components**:
```python
# k-space operators (1D)
k = 2 * pi * fftfreq(Nx, d=dx)
kappa = sinc((c_ref * k * dt / 2) / pi)
op_grad = 1j * k * kappa * exp( 1j * k * dx/2)   # Forward shift
op_div  = 1j * k * kappa * exp(-1j * k * dx/2)   # Backward shift

# Time loop - momentum equation
u -= (dt / rho0_staggered) * diff(p, op_grad)

# Time loop - mass conservation
duxdx = diff(u, op_div)
rho -= (dt * rho0) * duxdx * nonlinear_factor(rho)

# Time loop - equation of state
p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))
```

**MATLAB Reference 2D Implementation - Key Patterns**:
```matlab
% k-space operators (2D)
ddx_k_shift_pos = ifftshift( 1i * kx_vec .* exp( 1i * kx_vec * dx/2) );
ddx_k_shift_neg = ifftshift( 1i * kx_vec .* exp(-1i * kx_vec * dx/2) );
ddy_k_shift_pos = ifftshift( 1i * ky_vec .* exp( 1i * ky_vec * dy/2) );
ddy_k_shift_neg = ifftshift( 1i * ky_vec .* exp(-1i * ky_vec * dy/2) );

% Momentum equation (2D)
ux_sgx -= dt * rho0_sgx_inv .* real(ifft2( ddx_k_shift_pos .* kappa .* fft2(p) ))
uy_sgy -= dt * rho0_sgy_inv .* real(ifft2( ddy_k_shift_pos .* kappa .* fft2(p) ))

% Mass conservation (2D)
duxdx = real(ifft2( ddx_k_shift_neg .* kappa .* fft2(ux_sgx) ))
duydy = real(ifft2( ddy_k_shift_neg .* kappa .* fft2(uy_sgy) ))
rhox -= dt * rho0 * duxdx
rhoy -= dt * rho0 * duydy

% Equation of state (2D)
p = c0^2 * (rhox + rhoy + absorption_term - dispersion_term + nonlinear_term)
```

### Mathematical Foundation

The k-space pseudospectral method for coupled first-order acoustic equations:

**1D**:
- Momentum: ∂u/∂t = -(1/ρ₀) ∂p/∂x
- Mass: ∂ρ/∂t = -ρ₀ ∂u/∂x
- State: p = c₀² ρ

**2D**:
- Momentum: ∂ux/∂t = -(1/ρ₀) ∂p/∂x, ∂uy/∂t = -(1/ρ₀) ∂p/∂y
- Mass: ∂ρx/∂t = -ρ₀ ∂ux/∂x, ∂ρy/∂t = -ρ₀ ∂uy/∂y
- State: p = c₀² (ρx + ρy)

**3D**: Same pattern with ux, uy, uz and ρx, ρy, ρz

**Key Observation**: The density is split into dimensional components (rhox, rhoy, rhoz) which are then summed in the equation of state. This is a numerical technique to maintain the staggered grid structure.

---

## Part 2: N-Dimensional Design Strategy

### Design Principle: Dimension-Agnostic Operators

The cleanest generalization uses **lists of operators** indexed by dimension:

```python
# Current 1D (scalar operators)
op_grad = 1j * k * kappa * exp(1j * k * dx/2)
op_div = 1j * k * kappa * exp(-1j * k * dx/2)

# N-D generalization (list of operators)
op_grad = [op_grad_x, op_grad_y, op_grad_z][:ndim]
op_div = [op_div_x, op_div_y, op_div_z][:ndim]
u = [ux, uy, uz][:ndim]
rho_split = [rhox, rhoy, rhoz][:ndim]
```

### Minimal Code Changes Required

**1. Parse Grid Dimensions** (add to `_parse_inputs`):
```python
def _parse_inputs(kgrid, medium, sensor, xp):
    # Unpack grid dimensions (gracefully handle 1D/2D/3D)
    dims = []
    spacing = []
    for dim_name in ['Nx', 'Ny', 'Nz']:
        if hasattr(kgrid, dim_name):
            dims.append(int(getattr(kgrid, dim_name)))
            spacing.append(float(getattr(kgrid, dim_name.replace('N', 'd'))))

    Nt, dt = int(kgrid.Nt), float(kgrid.dt)
    ndim = len(dims)

    # Expand medium parameters to grid shape
    grid_shape = tuple(dims)
    c0 = _expand_to_grid(medium.sound_speed, grid_shape, xp)
    rho0 = _expand_to_grid(_attr(medium, 'density', 1000.0), grid_shape, xp)
    ...
```

**2. Build N-D Operators** (new function):
```python
def _build_kspace_operators(dims, spacing, dt, c_ref, xp):
    """Build k-space operators for each dimension."""
    ndim = len(dims)
    kappa_list = []
    op_grad_list = []
    op_div_list = []

    for i, (N, dx) in enumerate(zip(dims, spacing)):
        # k-space wavenumber for this dimension
        k = 2 * np.pi * xp.fft.fftfreq(N, d=dx)

        # Reshape to broadcast along dimension i
        # (1, 1, N, 1) for dimension 2 in 4D array, etc.
        shape = [1] * ndim
        shape[i] = N
        k = k.reshape(shape)

        # Time-staggering correction (k-space correction)
        kappa = xp.sinc((c_ref * k * dt / 2) / np.pi)

        # Gradient (forward shift) and divergence (backward shift) operators
        op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)
        op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)

        kappa_list.append(kappa)
        op_grad_list.append(op_grad)
        op_div_list.append(op_div)

    # Combine kappa from all dimensions (multiply)
    kappa_nd = kappa_list[0]
    for kappa_d in kappa_list[1:]:
        kappa_nd = kappa_nd * kappa_d

    return kappa_nd, op_grad_list, op_div_list
```

**3. N-D FFT and Spectral Diff** (generalize):
```python
def _spectral_diff(f, op, axis, xp, ndim):
    """Apply spectral operator along specific axis."""
    if ndim == 1:
        return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))
    elif ndim == 2:
        return xp.real(xp.fft.ifft2(op * xp.fft.fft2(f)))
    elif ndim == 3:
        return xp.real(xp.fft.ifftn(op * xp.fft.fftn(f)))
    else:
        raise ValueError(f"Unsupported dimension: {ndim}")
```

**4. Time Loop N-D** (replace scalar with loops):
```python
# Initialize velocity and split density for each dimension
u = [xp.zeros(grid_shape, dtype=float) for _ in range(ndim)]
rho_split = [xp.zeros(grid_shape, dtype=float) for _ in range(ndim)]

for t in range(Nt):
    # Momentum equation: du_i/dt = -(1/rho0) * dp/dx_i
    for i in range(ndim):
        grad_p = diff(p, op_grad_list[i], axis=i)
        u[i] -= (dt / rho0_staggered[i]) * grad_p
        u[i] = source_u_op[i](t, u[i])

    # Mass conservation: drho_i/dt = -rho0 * du_i/dx_i
    div_u_components = []
    for i in range(ndim):
        div_u_i = diff(u[i], op_div_list[i], axis=i)
        div_u_components.append(div_u_i)
        rho_split[i] -= (dt * rho0) * div_u_i * nonlinear_factor(sum(rho_split))
        rho_split[i] = source_p_op(t, rho_split[i])

    # Equation of state
    rho_total = sum(rho_split)
    div_u_total = sum(div_u_components)
    p = c0**2 * (rho_total + absorption(div_u_total) - dispersion(rho_total) + nonlinearity(rho_total))

    sensor_data[:, t] = p[mask]
```

---

## Part 3: Step-by-Step Implementation Plan

### Phase 2.1: Refactor for N-D (No PML)

**Goal**: Make the code dimension-agnostic while maintaining 1D parity

#### Step 1: Refactor Grid Parsing (1-2 hours)

**File**: `k-Wave/python/kWavePy.py`

**Changes**:
```python
def _parse_inputs(kgrid, medium, sensor, xp):
    """Parse inputs for N-dimensional grids."""
    # Extract dimensions and spacing
    dims, spacing = [], []
    for axis in ['x', 'y', 'z']:
        N_attr = f'N{axis}'
        d_attr = f'd{axis}'
        if hasattr(kgrid, N_attr) and getattr(kgrid, N_attr) is not None:
            dims.append(int(getattr(kgrid, N_attr)))
            spacing.append(float(getattr(kgrid, d_attr)))
        else:
            break  # Stop at first missing dimension

    grid_shape = tuple(dims)
    ndim = len(dims)
    Nt, dt = int(kgrid.Nt), float(kgrid.dt)

    # Expand medium parameters to N-D grid
    c0 = _expand_to_grid(medium.sound_speed, grid_shape, xp)
    rho0 = _expand_to_grid(_attr(medium, 'density', 1000.0), grid_shape, xp)

    # Sensor mask handling
    mask_raw = _attr(sensor, 'mask', None)
    if mask_raw is None:
        mask = xp.ones(grid_shape, dtype=bool)
    else:
        mask = xp.array(mask_raw, dtype=bool, copy=True)
        # MATLAB uses Fortran order (column-major)
        if mask.ndim == 1:
            mask = mask.reshape(grid_shape, order='F')
        elif mask.shape != grid_shape:
            raise ValueError(f"Sensor mask shape {mask.shape} doesn't match grid {grid_shape}")

    return grid_shape, dims, spacing, Nt, dt, c0, rho0, mask
```

**Validation**:
- Run existing 1D tests: `runtests('tests/test_1d_parity.m')`
- Expected: All tests still pass (backward compatible)

#### Step 2: Build N-D k-space Operators (2-3 hours)

**File**: `k-Wave/python/kWavePy.py`

**New Function**:
```python
def _build_kspace_operators(dims, spacing, dt, c_ref, xp):
    """Build dimension-specific k-space operators."""
    ndim = len(dims)
    op_grad_list = []
    op_div_list = []

    for axis_idx, (N, dx) in enumerate(zip(dims, spacing)):
        # Wavenumber for this dimension
        k = 2 * np.pi * xp.fft.fftfreq(N, d=dx)

        # Broadcast shape: expand along axis_idx
        shape = [1] * ndim
        shape[axis_idx] = N
        k = k.reshape(shape)

        # k-space correction
        kappa = xp.sinc((c_ref * k * dt / 2) / np.pi)

        # Operators
        op_grad_list.append(1j * k * kappa * xp.exp( 1j * k * dx/2))
        op_div_list.append(1j * k * kappa * xp.exp(-1j * k * dx/2))

    # Compute combined kappa for source corrections
    kappa_nd = op_grad_list[0] / (1j * 2 * np.pi * xp.fft.fftfreq(dims[0], d=spacing[0]).reshape([dims[0]] + [1]*(ndim-1)))
    # ... (need to extract kappa properly)

    return op_grad_list, op_div_list, kappa_nd
```

**Validation**:
- Unit test: verify operators match 1D case when ndim=1
- Check operator shapes are correct for 2D/3D

#### Step 3: Generalize Spectral Diff (1 hour)

**File**: `k-Wave/python/kWavePy.py`

**Replace**:
```python
# Old
def _spectral_diff(f, op, xp):
    return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

# New
def _spectral_diff(f, op, ndim, xp):
    """Apply spectral operator in N-D."""
    if ndim == 1:
        return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))
    elif ndim == 2:
        return xp.real(xp.fft.ifft2(op * xp.fft.fft2(f)))
    elif ndim == 3:
        return xp.real(xp.fft.ifftn(op * xp.fft.fftn(f)))
    else:
        raise ValueError(f"Unsupported ndim: {ndim}")
```

#### Step 4: Generalize Time Loop (3-4 hours)

**File**: `k-Wave/python/kWavePy.py`

**Key Changes**:
```python
def simulate(...):
    # ... parsing ...
    grid_shape, dims, spacing, Nt, dt, c0, rho0, mask = _parse_inputs(...)
    ndim = len(dims)

    # Build operators
    op_grad_list, op_div_list, source_kappa = _build_kspace_operators(dims, spacing, dt, c_ref, xp)
    diff = lambda f, op: _spectral_diff(f, op, ndim, xp)

    # Initialize fields
    p = xp.zeros(grid_shape, dtype=float)
    u = [xp.zeros(grid_shape, dtype=float) for _ in range(ndim)]
    rho_split = [xp.zeros(grid_shape, dtype=float) for _ in range(ndim)]

    # Staggered density (one per dimension)
    rho0_staggered = [_compute_staggered_density(rho0, axis, xp) for axis in range(ndim)]

    # Initialize velocity (staggered leapfrog)
    for axis in range(ndim):
        u[axis] += (dt / (2 * rho0_staggered[axis])) * diff(p, op_grad_list[axis])

    for t in range(Nt):
        # Momentum: du_i/dt = -(1/rho0) dp/dx_i
        for i in range(ndim):
            u[i] -= (dt / rho0_staggered[i]) * diff(p, op_grad_list[i])
            u[i] = source_u_op[i](t, u[i])

        # Mass: drho_i/dt = -rho0 du_i/dx_i
        div_u_components = []
        for i in range(ndim):
            div_u_i = diff(u[i], op_div_list[i])
            div_u_components.append(div_u_i)
            rho_split[i] -= dt * rho0 * div_u_i * nonlinear_factor(sum(rho_split))
            rho_split[i] = source_p_op(t, rho_split[i])

        # Equation of state
        rho_total = sum(rho_split)
        div_u_total = sum(div_u_components)
        p = c0**2 * (rho_total + absorption(div_u_total) - dispersion(rho_total) + nonlinearity(rho_total))

        # p0 override
        if t == 0 and p0_initial is not None:
            p = p0_initial.copy()
            for i in range(ndim):
                rho_split[i] = p0_initial / (c0**2 * ndim)
                u[i] = (dt / (2 * rho0_staggered[i])) * diff(p, op_grad_list[i])

        sensor_data[:, t] = p[mask]

    return {"sensor_data": _to_cpu(sensor_data), "pressure": _to_cpu(p)}
```

**Validation**:
- Run 1D tests: should still pass
- Test basic 2D simulation with homogeneous medium

#### Step 5: Update Absorption/Dispersion for N-D (2 hours)

**File**: `k-Wave/python/kWavePy.py`

**Changes to `_fractional_laplacian`**:
```python
def _fractional_laplacian(k_list, power, xp):
    """N-D fractional Laplacian: |k|^power where k is the magnitude."""
    # Compute k_mag = sqrt(kx^2 + ky^2 + kz^2)
    k_mag_sq = sum(xp.fft.fftshift(k)**2 for k in k_list)
    k_mag = xp.sqrt(k_mag_sq)
    return xp.fft.ifftshift(xp.where(k_mag == 0, 0, k_mag**power))
```

**Changes to absorption operators**: Accept `div_u_total` instead of single derivative

#### Step 6: Update Source Operators for N-D (2 hours)

**File**: `k-Wave/python/kWavePy.py`

**Changes**:
- `source_p_op`: Works on all dimensions (unchanged)
- `source_u_op`: Need separate operators for ux, uy, uz
- Return list of velocity source operators: `[source_ux_op, source_uy_op, source_uz_op]`

#### Step 7: Update MATLAB Wrapper (1 hour)

**File**: `k-Wave/kspaceFirstOrderPy.m`

**Changes**:
```matlab
% Validate inputs
if kgrid.dim > 3
    error('kspaceFirstOrderPy:UnsupportedDimension', 'Only 1D, 2D, 3D supported');
end

% Build kgrid dict with conditional fields
kgrid_args = {'Nx', int64(kgrid.Nx), 'dx', kgrid.dx, 'Nt', int64(kgrid.Nt), 'dt', kgrid.dt};
if kgrid.dim >= 2
    kgrid_args = [kgrid_args, {'Ny', int64(kgrid.Ny), 'dy', kgrid.dy}];
end
if kgrid.dim == 3
    kgrid_args = [kgrid_args, {'Nz', int64(kgrid.Nz), 'dz', kgrid.dz}];
end
k_py = py.dict(pyargs(kgrid_args{:}));

% Add velocity source components for 2D/3D
s_args = {...};
if kgrid.dim >= 2
    s_args = [s_args, {'uy', toNumpy(getField(source, {'uy'}, 0))}];
end
if kgrid.dim == 3
    s_args = [s_args, {'uz', toNumpy(getField(source, {'uz'}, 0))}];
end
```

#### Step 8: Create 2D/3D Shims (30 minutes)

**Files**:
- `tests/shims/kspaceFirstOrder2D.m`
- `tests/shims/kspaceFirstOrder3D.m`

**Content** (identical to 1D shim):
```matlab
function sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, varargin)
    sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin{:});
end
```

**Validation**: Unit tests can now use shims

#### Step 9: Comprehensive Testing (2-3 hours)

**Test Strategy**:

1. **1D Regression**: Ensure all existing 1D tests still pass
   ```bash
   arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "..."
   ```

2. **2D Basic**: Run simple 2D test (homogeneous medium, no absorption)
   ```matlab
   runtests('k-Wave/testing/unit/kspaceFirstOrderPy_binary_sensor_mask_2D.m')
   ```

3. **3D Basic**: Run simple 3D test
   ```matlab
   runtests('k-Wave/testing/unit/kspaceFirstOrderPy_binary_sensor_mask_3D.m')
   ```

4. **CI Integration**: Update `.github/workflows/test.yml` to include 2D/3D patterns

### Phase 2.2: PML Implementation (No N-D yet)

**Goal**: Add Perfectly Matched Layer (PML) absorption for boundaries in 1D first

#### Background: PML Theory

The PML is an absorbing boundary condition that prevents reflections. In k-Wave, it's implemented as a time-domain absorption applied in each dimension:

```
u_new = pml_coeff * (pml_coeff * u_old - dt * grad_p)
```

where `pml_coeff = exp(-alpha * distance_from_edge)` tapers from 1 (interior) to ~0 (boundary).

#### Step 1: Understand MATLAB PML Implementation (1 hour)

**File**: `k-Wave/private/getPML.m`

**Read and extract**:
- How PML coefficients are computed
- How they're applied in the time loop
- Staggered vs non-staggered grids

#### Step 2: Implement PML Coefficient Generation (2 hours)

**File**: `k-Wave/python/kWavePy.py`

**New Function**:
```python
def _get_pml_coeffs(N, pml_size, pml_alpha, dt, c_ref, dx, staggered, xp):
    """Generate PML absorption coefficients (matches MATLAB getPML.m)."""
    if pml_size == 0:
        return xp.ones(N, dtype=float)

    # Create distance profile from edges
    pml_left = xp.arange(pml_size, 0, -1, dtype=float) - (0.5 if staggered else 0)
    pml_right = xp.arange(1, pml_size + 1, dtype=float) - (0.5 if staggered else 0)

    # Apply absorption profile
    pml_left = pml_alpha * (pml_left / pml_size)**4
    pml_right = pml_alpha * (pml_right / pml_size)**4

    # Assemble full profile
    pml = xp.concatenate([pml_left, xp.zeros(N - 2*pml_size), pml_right])

    # Convert to time-domain coefficient
    return xp.exp(-pml * c_ref * dt / dx)
```

#### Step 3: Apply PML in 1D Time Loop (1-2 hours)

**File**: `k-Wave/python/kWavePy.py`

**Changes**:
```python
# Build PML operators
pml_size = int(_attr(medium, 'pml_size', 10))  # Default PML
pml_alpha = float(_attr(medium, 'pml_alpha', 2.0))

pml_u = _get_pml_coeffs(Nx, pml_size, pml_alpha, dt, c_ref, dx, staggered=True, xp=xp)
pml_rho = _get_pml_coeffs(Nx, pml_size, pml_alpha, dt, c_ref, dx, staggered=False, xp=xp)

# In time loop - apply PML to velocity
for i in range(ndim):
    grad_p = diff(p, op_grad_list[i])
    u[i] = pml_u * (pml_u * u[i] - (dt / rho0_staggered[i]) * grad_p)

# In time loop - apply PML to density
for i in range(ndim):
    div_u_i = diff(u[i], op_div_list[i])
    rho_split[i] = pml_rho * (pml_rho * rho_split[i] - dt * rho0 * div_u_i * nonlinear_factor(...))
```

**Validation**:
- Compare 1D with PML enabled vs MATLAB
- Check boundary reflections are suppressed

#### Step 4: Extend PML to N-D (2 hours)

**File**: `k-Wave/python/kWavePy.py`

**Changes**:
```python
# Build PML per dimension
pml_u_list = []
pml_rho_list = []
for axis_idx, (N, dx) in enumerate(zip(dims, spacing)):
    pml_u_axis = _get_pml_coeffs(N, pml_size, pml_alpha, dt, c_ref, dx, staggered=True, xp=xp)
    pml_rho_axis = _get_pml_coeffs(N, pml_size, pml_alpha, dt, c_ref, dx, staggered=False, xp=xp)

    # Broadcast to N-D
    shape = [1] * ndim
    shape[axis_idx] = N
    pml_u_list.append(pml_u_axis.reshape(shape))
    pml_rho_list.append(pml_rho_axis.reshape(shape))

# Apply in time loop
for i in range(ndim):
    u[i] = pml_u_list[i] * (pml_u_list[i] * u[i] - (dt / rho0_staggered[i]) * grad_p)
    rho_split[i] = pml_rho_list[i] * (pml_rho_list[i] * rho_split[i] - dt * rho0 * div_u_i * ...)
```

### Phase 2.3: Test-Driven Feature Completion

**Goal**: Achieve 100% pass rate on 2D/3D unit tests

#### Step 1: Run CI Test Suite (1 hour)

**Command**:
```bash
# Enable shims and run 2D tests
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "..."
```

**Analyze failures**:
- Categorize by feature (sensor types, source modes, physics)
- Prioritize by frequency and importance

#### Step 2: Implement Missing Features Iteratively (10-20 hours)

**Expected Features** (from 1D experience):
- Cartesian sensor masks (interpolation)
- Rectangle sensor regions
- Time reversal
- Different source modes
- Heterogeneous absorption

**Approach**:
1. Pick one failing test
2. Identify missing feature
3. Implement feature minimally
4. Verify test passes
5. Repeat

#### Step 3: Optimize for Readability (2-3 hours)

**After all tests pass**:
- Refactor for clarity
- Remove code duplication
- Add targeted comments for non-obvious sections
- Verify line count stays minimal (~100-150 lines target)

---

## Part 4: Code Quality and Clean Code Principles

### Readability Guidelines

**Current Strengths** (preserve these):
1. Descriptive function names (`_spectral_diff`, `_fractional_laplacian`)
2. Clear parameter names (`op_grad`, `op_div` not `k1`, `k2`)
3. Inline documentation via variable names (code is self-documenting)
4. Single responsibility per function
5. Minimal abstraction - math is visible

**Improvements for N-D**:
1. **Avoid "clever" indexing**: Use explicit loops over dimensions
   ```python
   # Good (readable)
   for i in range(ndim):
       u[i] -= (dt / rho0_staggered[i]) * diff(p, op_grad_list[i])

   # Bad (clever but obscure)
   u = [u[i] - (dt / rho0_staggered[i]) * diff(p, op_grad_list[i]) for i in range(ndim)]
   ```

2. **Broadcast operations explicitly**: Make array shapes clear
   ```python
   # Good
   shape = [1] * ndim
   shape[axis] = N
   k = k.reshape(shape)  # Now broadcasts correctly

   # Bad
   k = k[..., None, None]  # What are these dimensions?
   ```

3. **Name magic numbers**: Avoid unexplained constants
   ```python
   # Good
   DEFAULT_PML_SIZE = 10
   DEFAULT_PML_ALPHA = 2.0

   # Bad
   pml_size = 10  # Where did this come from?
   ```

4. **Physics variable names match literature**:
   - `rho_split` (not `density_components`)
   - `op_grad`, `op_div` (not `k_forward`, `k_backward`)
   - `kappa` for k-space correction (standard k-Wave term)

### Cognitive Load Management

**Keep complexity per function low**:
- `simulate()`: Main orchestration (~30 lines)
- `_build_kspace_operators()`: Operator construction (~20 lines)
- `_parse_inputs()`: Input validation (~30 lines)
- Time loop: ~40 lines (unavoidable physics complexity)

**Total target: 120-150 lines** (up from current 65, but handles 3 dimensions)

### Testing Strategy

**Regression Protection**:
1. All 84 existing 1D tests must pass
2. Add 2D/3D equivalents of key 1D tests
3. CI runs both Python and MATLAB baselines

**Parity Threshold**:
- 1D: <1e-15 (already achieved)
- 2D: <1e-14 (target, may need investigation if not met)
- 3D: <1e-14 (target)

**Test Organization**:
```
k-Wave/testing/unit/
  kspaceFirstOrderPy_binary_sensor_mask_2D.m
  kspaceFirstOrderPy_binary_sensor_mask_3D.m
  kspaceFirstOrderPy_absorption_2D.m
  kspaceFirstOrderPy_pml_2D.m
  ...
```

---

## Part 5: Implementation Checklist and Timeline

### Week 1: N-D Refactoring (No PML)

- [ ] Day 1-2: Refactor `_parse_inputs` for N-D (Step 2.1.1)
- [ ] Day 2-3: Build N-D k-space operators (Step 2.1.2-3)
- [ ] Day 3-4: Generalize time loop (Step 2.1.4)
- [ ] Day 4-5: Update absorption/sources for N-D (Step 2.1.5-6)
- [ ] Day 5: Update MATLAB wrapper, create shims (Step 2.1.7-8)
- [ ] Day 6-7: Test 1D regression, basic 2D/3D (Step 2.1.9)

**Milestone**: 1D tests still pass, basic 2D/3D simulations run

### Week 2: PML and Feature Completion

- [ ] Day 8: Study MATLAB PML (Step 2.2.1)
- [ ] Day 9-10: Implement 1D PML (Step 2.2.2-3)
- [ ] Day 11: Extend PML to N-D (Step 2.2.4)
- [ ] Day 12-14: Run full test suite, implement missing features (Step 2.3.1-2)

**Milestone**: 2D/3D tests passing with PML enabled

### Week 3: Refinement and Documentation

- [ ] Day 15-16: Code review and refactoring (Step 2.3.3)
- [ ] Day 17: Write unit tests for N-D specific features
- [ ] Day 18: Update CI to include 2D/3D tests
- [ ] Day 19: Performance benchmarking (Python vs MATLAB)
- [ ] Day 20: Documentation and examples

**Milestone**: Production-ready N-D implementation

---

## Part 6: Risk Assessment and Mitigation

### Technical Risks

**Risk 1: Array indexing/broadcasting bugs in N-D**
- **Likelihood**: High (common with numpy broadcasting)
- **Impact**: High (silent wrong answers)
- **Mitigation**:
  - Extensive shape assertion tests
  - Visual inspection of intermediate results
  - Compare 2D case with known working MATLAB code step-by-step

**Risk 2: Performance degradation in N-D**
- **Likelihood**: Medium
- **Impact**: Medium (slower than expected)
- **Mitigation**:
  - Profile early (use `cProfile`)
  - Ensure FFT operations use optimal sizes
  - Pre-allocate arrays, avoid copies

**Risk 3: PML implementation differs from MATLAB**
- **Likelihood**: Medium (subtle numerical differences)
- **Impact**: High (boundary reflections, wrong physics)
- **Mitigation**:
  - Unit test PML coefficients against MATLAB output
  - Visual inspection of boundary behavior
  - Compare PML absorption profiles

**Risk 4: Sensor mask handling for 2D/3D**
- **Likelihood**: Medium
- **Impact**: Medium (sensor data extraction broken)
- **Mitigation**:
  - Test all sensor mask types (binary, Cartesian, rectangle)
  - Start with binary only (current support)
  - Add Cartesian later if needed

### Process Risks

**Risk 5: Scope creep - too many features**
- **Likelihood**: High
- **Impact**: High (delayed delivery)
- **Mitigation**:
  - Strict scope: N-D + PML only
  - Defer advanced features (time reversal, etc.) to Phase 3
  - Use test-driven approach - only implement what tests require

**Risk 6: Breaking 1D backward compatibility**
- **Likelihood**: Medium
- **Impact**: Critical (regression)
- **Mitigation**:
  - Run 1D tests after every change
  - CI enforces 1D test passing
  - Keep 1D-specific code path if needed

---

## Part 7: Success Criteria

### Functional Requirements

1. **1D Backward Compatibility**: All 84 existing 1D tests pass with <1e-15 error
2. **2D Support**: Basic 2D simulations run with binary sensor masks
3. **3D Support**: Basic 3D simulations run with binary sensor masks
4. **PML Boundaries**: PML implemented and working in all dimensions
5. **Physics Parity**: Absorption, dispersion, nonlinearity work in N-D

### Non-Functional Requirements

1. **Code Quality**:
   - Python: 120-150 lines (target)
   - MATLAB wrapper: 50-60 lines (target)
   - All functions <30 lines
   - Clear variable names (no abbreviations except standard physics terms)

2. **Performance**:
   - 2D: Within 2x of MATLAB reference
   - 3D: Within 2x of MATLAB reference
   - CuPy acceleration works (optional, already implemented)

3. **Documentation**:
   - All public functions have docstrings
   - Implementation plan (this document) complete
   - Example scripts for 2D and 3D

4. **Testing**:
   - 2D unit tests: At least 5 tests covering core features
   - 3D unit tests: At least 3 tests covering core features
   - CI runs 1D/2D/3D test matrix

### Definition of Done

**Phase 2 Complete** when:
- [ ] All 1D tests pass (<1e-15 error)
- [ ] 2D homogeneous test passes (<1e-14 error)
- [ ] 3D homogeneous test passes (<1e-14 error)
- [ ] PML implemented and tested
- [ ] Code reviewed and refactored for clarity
- [ ] CI updated to run 2D/3D tests
- [ ] This implementation plan archived as reference

---

## Part 8: Open Questions and Design Decisions

### Question 1: How to handle staggered density in N-D?

**Background**: In 1D, density is staggered by averaging neighboring grid points. In 2D/3D, should we stagger per-dimension?

**MATLAB Approach**: Uses per-dimension staggering (rhox, rhoy, rhoz)

**Decision**: Follow MATLAB - use split density components, each staggered in its own dimension

**Rationale**: Maintains numerical stability and grid staggering properties

### Question 2: Should we support Cartesian sensor masks in Phase 2?

**Background**: Current 1D implementation only supports binary masks. 2D/3D often use Cartesian coordinates.

**Trade-off**:
- Pro: More user-friendly, matches MATLAB API
- Con: Adds interpolation complexity (~20 lines)

**Decision**: Defer to Phase 2.3 (test-driven). Implement only if tests require it.

**Rationale**: Keep initial scope minimal, add features as needed

### Question 3: How to handle reference sound speed (`c_ref`)?

**Background**: MATLAB allows `medium.sound_speed_ref` to set reference for k-space operator. Default is `max(c0)`.

**Decision**: Auto-compute as `max(c0)` initially, add parameter support if tests fail

**Rationale**: Simplicity first, configurability later

### Question 4: Should PML size be per-dimension or global?

**Background**: MATLAB allows different PML sizes per dimension (PMLSize = [10 20 15])

**Decision**: Start with global PML size, add per-dimension support if needed

**Rationale**: Minimize parameter complexity for initial implementation

---

## Part 9: Future Work (Phase 3+)

### Phase 3: Advanced Features (Post N-D)

**Not included in this plan**:
1. Time reversal reconstruction
2. Advanced sensor types (directional, frequency response)
3. Elastic wave simulation
4. Nonlinear absorption models
5. GPU optimization (CuPy already works, but profiling/tuning)

### Phase 4: Performance Optimization

**Potential optimizations**:
1. JIT compilation (Numba)
2. Fused operations (reduce FFT calls)
3. Multi-GPU support
4. Sparse sensor masks (don't record full grid)

### Phase 5: API Enhancement

**Usability improvements**:
1. Progress bars for long simulations
2. Checkpoint/resume capability
3. Automatic grid sizing
4. Parameter validation and helpful error messages

---

## Conclusion

This implementation plan provides a clear, incremental path to N-dimensional support while maintaining the minimalistic "academic code golf" philosophy. The key insights are:

1. **Dimension-agnostic operators**: Lists of per-dimension operators generalize cleanly
2. **Backward compatibility first**: Ensure 1D tests always pass
3. **Test-driven development**: Let failing tests guide feature implementation
4. **Clean code over cleverness**: Explicit loops beat list comprehensions for readability
5. **Incremental delivery**: Week 1 gets basic N-D working, Week 2 adds PML and features

**Estimated Total Effort**: 15-20 person-days (3-4 weeks part-time)

**Key Deliverables**:
- N-D Python backend (~150 lines)
- Updated MATLAB wrapper (~60 lines)
- 2D/3D test shims
- Unit tests for 2D/3D
- CI integration
- This implementation plan

**Success Metric**: 100% of critical tests passing with <1e-14 parity between Python and MATLAB in all dimensions.
