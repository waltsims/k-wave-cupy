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
- [ ] **Verify 1D Parity**: 
    - Run a simple pulse propagation test. Compare `norm(p_matlab - p_python) < tol`.

### Phase 2: Feature Parity (2D/3D & PML)
*Goal: Generalize the engine to support standard k-Wave features.*

- [ ] **Test Adapter Implementation**:
    - Create `kWaveTesterPy.m` that runs existing MATLAB unit tests but injects `kspaceFirstOrderPy` as the solver.
- [ ] **2D/3D Generalization**:
    - Expand `kWavePy.py` to handle N-dimensions.
    - Validate row/column major indexing for 2D/3D FFTs.
- [ ] **PML Implementation**:
    - Implement the Perfectly Matched Layer (essential for non-periodic simulation).
    - Verify against `acousticFieldPropagator` tests.

### Phase 3: Acceleration
- [x] **CuPy Backend**:
    - Add logic to switch `xp = cupy` if available.
    - Verify on GPU node.

## Next Steps
1. Run the MATLAB tests (`tests/test_interop_sanity.m`, `tests/test_1d_parity.m`, and `testing/unit/test_interface_1D.m`) to confirm interop, parity, and the 1D wrapper.
2. Begin Phase 2: add a MATLAB test adapter `kWaveTesterPy.m` that swaps in the Python solver for existing unit tests.
