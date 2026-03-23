# NumPy/CuPy Parity Test Suite via Reference Data

## Context

The Python backend (`kWavePy.py`) uses `self.xp` to dispatch operations to NumPy or CuPy. We need a **transitive parity chain**: MATLAB ≈ NumPy ≈ CuPy.

- **CI (MATLAB available)**: Generate MATLAB reference data, then compare NumPy against it
- **GPU device (no MATLAB)**: Compare CuPy output against NumPy output directly
- If reference data files exist, NumPy is validated against MATLAB; if CuPy is available, CuPy is validated against NumPy

## Overview of Changes

| File | Action | Purpose |
|------|--------|---------|
| `tests/test_parity.py` | NEW | Main pytest file — 252 parametrized scenarios |
| `tests/conftest.py` | NEW | sys.path setup for kWavePy import |
| `k-Wave/testing/unit/generate_reference_data_2D.m` | NEW | Generate 2D MATLAB reference data (84 scenarios) |
| `k-Wave/testing/unit/generate_reference_data_3D.m` | NEW | Generate 3D MATLAB reference data (84 scenarios) |
| `pyproject.toml` | MODIFY | Add `dev` dependencies (pytest, scipy) |
| `.github/workflows/run_tests.yml` | MODIFY | Add 2D/3D reference data generation + Python parity job |

## Detailed Design

### 1. `tests/test_parity.py` — Main Test File

**Scenario parameterization** — Combinatorial product of 4 axes:

| Axis | Values | Count |
|------|--------|-------|
| Linearity | `linear`, `nonlinear` | 2 |
| Absorption | `lossless`, `lossy` (alpha_power=1.5), `stokes` (alpha_power=2) | 3 |
| Source | `p0`, `p_add`, `p_add_nocorr`, `p_dirichlet`, `ux_add`, `ux_add_nocorr`, `ux_dirichlet` | 7 |
| Medium | `homogeneous`, `heterogeneous` | 2 |

**Per dimension**: 2 × 3 × 7 × 2 = **84 scenarios**. Across 1D/2D/3D = **252 total**.

Additionally, targeted tests for features not covered by the combinatorial matrix:
- `uy` velocity source (2D/3D only, one scenario each) — validates multi-axis velocity dispatch
- `uz` velocity source (3D only, one scenario) — validates z-axis velocity
- `sensor.record_start_index` (one scenario per dim) — validates delayed recording
- Velocity recording with `sensor.record = ('p', 'ux')` (1D), `('p', 'ux', 'uy')` (2D), `('p', 'u')` (3D) — validates velocity output parity across backends

**Extra scenarios: 9** → **261 total**.

Scenario parameters match the existing MATLAB test (`kspaceFirstOrder1D_compare_plane_waves.m`) exactly:
- c0=1500, rho0=1000, BonA0=10, alpha0=5, alpha_power=1.5
- c1=2000, rho1=1200, BonA1=5, alpha1=2 (heterogeneous)
- source_strength=5e6, source_freq=200, CFL=0.1

**Grid sizes** (small for speed):
- 1D: Nx=64, dx=1, Nt=600, PML disabled (matches MATLAB test exactly)
- 2D: Nx=Ny=32, dx=dy=1, Nt=50, PML enabled (size=10, alpha=2)
- 3D: Nx=Ny=Nz=16, dx=dy=dz=1, Nt=30, PML enabled (size=10, alpha=2)

**Source/sensor in N-D**: Point source and point sensor, generalized:
- 1D source at `x=4`, sensor at `x=Nx-10` (matching MATLAB's `PML_size - 6` and `Nx - PML_size`)
- 2D source at `(4, Ny//2)`, sensor at `(Nx-10, Ny//2)`
- 3D source at `(4, Ny//2, Nz//2)`, sensor at `(Nx-10, Ny//2, Nz//2)`

**Helper function** `build_scenario(ndim, linearity, absorption, source_type, medium_type)`:
- Returns `(kgrid_dict, medium_dict, source_dict, sensor_dict)` ready for `kWavePy.simulate_from_dicts()`
- Constructs arrays as NumPy with Fortran order where needed for multi-D
- Heterogeneous media: step function at `x = N//2` (extended along all other dims)

**Three test functions**:

```python
def available_backends():
    """Returns ["cpu"] or ["cpu", "gpu"] depending on CuPy availability."""
    return ["cpu", "gpu"] if cp is not None else ["cpu"]

@pytest.mark.parametrize("scenario", ALL_SCENARIOS, ids=scenario_ids)
def test_numpy_valid(scenario):
    """NumPy backend produces valid (finite, non-trivial) output."""
    result = run_scenario(scenario, backend="cpu")
    assert np.all(np.isfinite(result["sensor_data"]))
    assert np.max(np.abs(result["sensor_data"])) > 0

@pytest.mark.parametrize("scenario", ALL_SCENARIOS, ids=scenario_ids)
@pytest.mark.parametrize("backend", available_backends())
def test_matlab_parity(scenario, backend):
    """Each available backend (NumPy and CuPy) matches MATLAB reference."""
    ref_data = load_reference_data(scenario)
    if ref_data is None:
        pytest.skip("No MATLAB reference data available")
    result = run_scenario(scenario, backend=backend)
    rel_error = compute_rel_error(ref_data, result["sensor_data"])
    assert rel_error < MATLAB_PARITY_THRESHOLD

@pytest.mark.parametrize("scenario", ALL_SCENARIOS, ids=scenario_ids)
@pytest.mark.skipif(cp is None, reason="CuPy not available")
def test_numpy_cupy_parity(scenario):
    """CuPy output matches NumPy (fallback when no MATLAB reference data)."""
    np_result = run_scenario(scenario, backend="cpu")
    cp_result = run_scenario(scenario, backend="gpu")
    rel_error = compute_rel_error(np_result["sensor_data"], cp_result["sensor_data"])
    assert rel_error < CUPY_PARITY_THRESHOLD
```

**Key insight**: `test_matlab_parity` is parametrized over `backend`, so when CuPy is available AND reference data exists, CuPy is validated directly against MATLAB — not just transitively via NumPy. The `test_numpy_cupy_parity` serves as a fallback for GPU machines without reference data.

**Reference data loading** — `load_reference_data(scenario)`:
- Looks for `.mat` files in `k-Wave/testing/unit/`:
  - `kspaceFirstOrder1D_reference_data.mat` (existing, 84 test results)
  - `kspaceFirstOrder2D_reference_data.mat` (new)
  - `kspaceFirstOrder3D_reference_data.mat` (new)
- Uses `scipy.io.loadmat()` to read `.mat` files
- Returns `None` if file doesn't exist (test skipped)
- Maps scenario index to the correct entry in the reference cell array

**Parity thresholds**:
- `MATLAB_PARITY_THRESHOLD = 1e-14` (matches existing 1D test; applies to both NumPy and CuPy vs MATLAB)
- `CUPY_PARITY_THRESHOLD = 1e-12` (conservative fallback; cuFFT vs pocketfft may differ at ULP level)

### 2. `tests/conftest.py` — Path Setup

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "k-Wave" / "python"))
```

### 3. MATLAB Reference Data Generators (2D/3D)

**`k-Wave/testing/unit/generate_reference_data_2D.m`** and **`generate_reference_data_3D.m`**:

Same 84-scenario structure as `kspaceFirstOrder1D_compare_plane_waves.m`, adapted for 2D/3D:
- 2D: `kWaveGrid(32, 1, 32, 1)`, point source at `(4, 16)`, sensor at `(22, 16)`
- 3D: `kWaveGrid(16, 1, 16, 1, 16, 1)`, point source at `(4, 8, 8)`, sensor at `(6, 8, 8)`
- Same physics parameter values (c0, rho0, BonA0, alpha0, etc.)
- Same source/medium configuration logic (switch on test_num)
- PML enabled (size=10 for 2D, size=4 for 3D to fit in 16-point grid)
- Saves `reference_results` cell array to `kspaceFirstOrder{2D,3D}_reference_data.mat`
- Uses default `save()` format (not `-v7.3`) for scipy compatibility

**Note**: 3D grid is 16×16×16 with PML size=4 (since PML_size must be < N/2). This matches the Python test grid.

### 4. `pyproject.toml` — Dependencies

```toml
[optional-dependencies]
cuda = ["cupy>=13.6.0"]
dev = ["pytest>=8.0", "scipy>=1.10"]
```

### 5. `.github/workflows/run_tests.yml` — CI Changes

Add two new components:

**a) Extend reference data generation** to include 2D/3D:
```yaml
generate-reference-data:
  # Existing 1D generation PLUS new 2D/3D scripts
  steps:
    - ... (existing steps)
    - name: Generate 2D/3D reference data
      run: |
        addpath('k-Wave');
        cd('k-Wave/testing/unit');
        generate_reference_data_2D;
        generate_reference_data_3D;
    - name: Upload reference data
      uses: actions/upload-artifact@v4
      with:
        name: reference_data
        path: |
          k-Wave/testing/unit/kspaceFirstOrder1D_reference_data.mat
          k-Wave/testing/unit/kspaceFirstOrder2D_reference_data.mat
          k-Wave/testing/unit/kspaceFirstOrder3D_reference_data.mat
```

**b) New Python parity test job**:
```yaml
python-parity-tests:
  name: Python parity tests
  runs-on: ubuntu-latest
  needs: generate-reference-data
  steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v5
      with:
        python-version: '3.10'
    - run: pip install numpy scipy pytest
    - name: Download reference data
      uses: actions/download-artifact@v4
      with:
        name: reference_data
        path: k-Wave/testing/unit
    - run: pytest tests/test_parity.py -v
```

This job runs without MATLAB — it only needs the reference data artifacts. The `test_matlab_numpy_parity` tests compare NumPy against downloaded reference data. The `test_numpy_cupy_parity` tests auto-skip (no CuPy in CI).

## Feature Coverage vs k-Wave Examples

Our scenarios cover **every feature that kWavePy.py implements** (including velocity recording from Plan 2):
- All source types: p0, time-varying pressure (3 modes), velocity ux/uy/uz (3 modes each)
- All medium types: homogeneous, heterogeneous (c, rho, alpha, BonA)
- All absorption models: lossless, power-law, Stokes
- Nonlinearity (B/A)
- PML boundaries
- Binary sensor masks
- sensor.record_start_index
- Velocity recording: `sensor.record = ('p', 'ux', ...)` — compares velocity outputs across backends

Features from k-Wave examples that are **out of scope** (not yet in kWavePy.py):
Cartesian sensor masks, kWaveArray/kWaveTransducer, axisymmetric mode, elastic waves, thermal diffusion, single precision, directional sensors.

## Dependency on velocity-recording.md

This plan depends on velocity recording being implemented first. The velocity recording scenarios in the parity test compare `result['ux']`, `result['uy']`, `result['uz']` outputs across backends.

## Verification

| Environment | What runs | Expected |
|---|---|---|
| **Local (no ref data, no GPU)** | `test_numpy_valid` only | 252 pass, MATLAB/CuPy tests skip |
| **CI (ref data, no GPU)** | `test_numpy_valid` + `test_matlab_parity[cpu]` | 252 + 252 = 504 pass |
| **GPU (no ref data)** | `test_numpy_valid` + `test_numpy_cupy_parity` | 252 + 252 = 504 pass |
| **GPU (with ref data)** | All three tests, both backends | 252 + 504 + 252 = 1008 pass |

Commands:
```bash
# Local: numpy-only smoke test
pytest tests/test_parity.py::test_numpy_valid -v

# CI: after downloading reference data artifacts
pytest tests/test_parity.py -v

# GPU machine (with ref data copied from CI):
pip install cupy-cuda12x scipy pytest
pytest tests/test_parity.py -v
```
