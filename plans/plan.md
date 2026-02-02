# Plan for CuPy-Accelerated k-Wave Backend (Minimalistic)

## Goal
Create a **minimalistic, "academic code golf" style** Python backend for k-Wave. 
- **Compute Engine**: Pure Python/NumPy/CuPy PSTD solver.
- **Integration**: Direct in-memory calls from MATLAB (`kspaceFirstOrderPy.m` -> `py.kWavePy.simulate`).
- **Philosophy**: Simple, interpretable code. No bloat. The code is the documentation.
- **Testing**: Verify using the existing MATLAB test suite.

## Background: Numerical Model (from k-Wave Manual)
The compute engine implements the **k-space pseudospectral method** for the coupled first-order acoustic equations.

### Governing Equations
$$ \frac{\partial p}{\partial t} = -\rho_0 c_0^2 \nabla \cdot \mathbf{u} $$
$$ \frac{\partial \mathbf{u}}{\partial t} = -\frac{1}{\rho_0} \nabla p $$

### Discretization
- **Spatial**: Fourier collocation spectral method.
  $$ \frac{\partial f}{\partial x} = \mathcal{F}^{-1} \{ i k_x \mathcal{F} \{ f \} \} $$
  *(Refinement: Add k-space correction $\text{sinc}(k \Delta x/2)$ later)*
- **Temporal**: Staggered Leapfrog.
  1. $\mathbf{u}^{n+1/2} = \mathbf{u}^{n-1/2} - \Delta t \rho_0^{-1} \nabla p^n$
  2. $p^{n+1} = p^n - \Delta t \rho_0 c_0^2 \nabla \cdot \mathbf{u}^{n+1/2}$

## Phases


### Phase 1: The "Steel Thread" (Interop & 1D)
*Goal: Prove the data travels correctly and the simplest physics works.*

- [x] **Risk Spike: Data Layout**: 
    - Create `tests/test_interop_sanity.m`: Send a non-symmetric 2D array to Python, modify it (e.g., `A[0,1] = 99`), return it, and verify indices match MATLAB expectations (Row vs Column major check).
- [x] **1D Interface Test (Failing)**: 
    - Create `testing/unit/test_interface_1D.m` invoking `kspaceFirstOrderPy` with a 1D grid. Assert it fails cleanly (TDD).
- [x] **Minimal 1D Engine (Python)**: 
    - Implement `kWavePy.simulate` to handle 1D arrays `(Nx, 1)`.
    - Implement minimal 1D spectral derivatives `d/dx = ifft(ik * fft)`.
- [x] **MATLAB Wrapper (1D Support)**: 
    - Update `kspaceFirstOrderPy.m` to detect 1D inputs and marshal them correctly (handling MATLAB's trailing singleton dimensions).
- [x] **Verify 1D Parity**: 
    - Run a simple pulse propagation test. Compare `norm(p_matlab - p_python) < tol`.
    - Achieved parity < 1e-15 by aligning sinc normalization, k-space operator ordering, and time loop recording.

### Phase 1.5: The "Zen Garden" (Refactoring)
*Goal: Ensure code is minimalistic, expressive, and easily interpretable.*
- [x] **Academic Code Golf**:
    - Reduced `kWavePy.py` to ~65 lines (was 135).
    - Reduced `kspaceFirstOrderPy.m` to ~40 lines (was 108).
    - Unified validation and data loading logic.
    - Verified 1D parity logic preservation.
- [ ] **Initial Perf Eval**: Compare runtimes between python and matlab impl.

### Phase 2: Feature Parity via Test-Based Development
*Goal: Generalize the engine to support standard k-Wave features using the existing test suite.*

#### Strategy: The "Shim" Architecture
Instead of modifying existing tests, we will inject a "Shim" path that redirects standard k-Wave function calls (e.g., `kspaceFirstOrder1D`) to our Python wrapper (`kspaceFirstOrderPy`).

- [x] **1D Shim Validation**:
    - Create `tests/shims/kspaceFirstOrder1D.m` which simply calls `kspaceFirstOrderPy`.
    - Relax `kspaceFirstOrderPy.m` argument parsing to ignore unsupported flags (like `PMLSize`, `PlotSim`).
    - **Verify**: Run `kspaceFirstOrder1D_check_source_scaling_p.m` with shims enabled. It should pass or fail only on numerics.
    - ✅ Test passed! Shim architecture validated successfully.
- [ ] **N-Dimensional Python Upgrade**:
    - Refactor `kWavePy.py` to handle N-dimensions dynamically (generic `op_grad` and `op_div` lists).
    - Ensure `Nx, Ny, Nz` unpacking handles missing dimensions gracefully.
    - **Verify**: Ensure 1D tests still pass with the N-D code.
- [ ] **2D/3D Shims**:
    - Add `tests/shims/kspaceFirstOrder2D.m` and `kspaceFirstOrder3D.m`.
    - **Verify**: Run `kspaceFirstOrder2D_check_source_scaling_p.m` and 3D equivalents.
- [ ] **Missing Physics Implementation**:
    - Implement **PML** (Perfectly Matched Layer) to fix boundary reflection failures.
    - Implement **Loss** (Power Law Absorption) to fix regression tests.
    - Implement **Nonlinearity** (BonA) to fix high-intensity tests.
- [ ] **CI Integration**:
    - Update GitHub Actions to run the full test suite twice: once normally, and once with the Shim path injected.

### Phase 3: Acceleration
- [x] **CuPy Backend**:
    - Add logic to switch `xp = cupy` if available.
    - Verify on GPU node.

## Next Steps
1. Execute the **1D Shim Validation** (create shim, update parser).
2. Refactor `kWavePy.py` for N-Dimensions.
3. Enable 2D/3D shims and begin implementing missing physics (PML first).
