# Testing Strategy for Python/CuPy Backend

## Overview

The k-Wave Python/CuPy backend is validated using a **dual-path testing strategy** that runs tests against both the standard MATLAB implementation and the Python backend simultaneously.

## CI Testing Architecture

Uses a **matrix strategy** to run unit tests with both backends in parallel:

```yaml
matrix:
  include:
    - backend: matlab
      shim_path: ""
      test_pattern: "kspaceFirstOrder1D"  # Temporarily limited to 1D for fast iteration
      artifact_suffix: ""
    - backend: python-1d
      shim_path: "tests/shims"
      test_pattern: "kspaceFirstOrder1D"
      artifact_suffix: "_python_backend"
```

### Matrix Configuration 1: Standard MATLAB (Baseline)
- **Purpose**: Ensure standard MATLAB functionality remains intact
- **Scope**: **1D tests only during development** (5 tests, ~minutes vs 2.5 hours for full suite)
- **Path**: Standard k-Wave functions only
- **Artifact**: `unit_test_results`

### Matrix Configuration 2: Python Backend (Development)
- **Purpose**: Track Python backend implementation progress
- **Scope**: 1D tests only (`kspaceFirstOrder1D*`)
- **Path**: Shims intercept calls and route to Python backend
- **Artifact**: `unit_test_results_python_backend`

### Optimization for Fast Iteration

**Current (Development Phase):**
- Both backends run only 1D tests (5 tests)
- Regression tests run only on `main` branch
- ~5-10 minute CI runs for fast feedback

**Future (After 1D Completion):**
- Restore `test_pattern: ""` for full test suite
- Remove regression test branch condition
- Enable 2D/3D testing

### Benefits of Matrix Approach
- **DRY**: Single job definition runs with multiple configurations
- **Parallel Execution**: Both backends tested simultaneously
- **Easy Expansion**: Adding 2D/3D is just another matrix entry
- **Clear Separation**: `fail-fast: false` ensures both runs complete independently

## Shim Architecture

The shim mechanism allows existing k-Wave tests to transparently use the Python backend:

```matlab
% tests/shims/kspaceFirstOrder1D.m
function sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, varargin)
    % Redirect to Python backend
    sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin{:});
end
```

When `tests/shims` is added to the MATLAB path **before** `k-Wave`, function calls are intercepted and routed to the Python implementation.

## Development Workflow

### Phase 1: 1D Feature Completion
1. Run CI to identify failing 1D tests with Python backend
2. Implement missing features (PML, absorption, etc.)
3. Iterate until all 1D tests pass
4. **Goal**: 100% parity between MATLAB and Python for 1D

### Phase 2: 2D/3D Expansion
1. Add `kspaceFirstOrder2D.m` and `kspaceFirstOrder3D.m` shims
2. Update Python engine for N-dimensional support
3. Update CI to run 2D/3D tests
4. Iterate until full dimensional parity achieved

### Phase 3: Performance Optimization
1. Enable CuPy GPU acceleration
2. Benchmark against C++/CUDA backends
3. Optimize bottlenecks

## Running Tests Locally

### Standard MATLAB tests
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('k-Wave'); cd('k-Wave/testing/unit'); runUnitTests();"
```

### Python backend tests (with shim)
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('k-Wave'); addpath('tests/shims'); cd('k-Wave/testing/unit'); runUnitTests('kspaceFirstOrder1D');"
```

**Important**: Add k-Wave BEFORE shims. In MATLAB, the last `addpath()` goes to the front of the search path, so this ensures shims are found first.

### Single shim validation test
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "run('tests/run_shim_validation.m')"
```

## Test Status Tracking

- **Baseline**: Standard MATLAB tests must continue passing
- **Python Backend**: Track implementation progress via test pass rate
- **No Regression**: Python backend failures don't block CI (separate job)
- **Visibility**: Both test results uploaded as artifacts for comparison

## Benefits

1. **Safe Development**: Python backend work can't break MATLAB functionality
2. **Clear Progress**: Test pass rate shows implementation completeness
3. **Early Detection**: CI catches regressions in either implementation
4. **Minimal Overhead**: Shim mechanism requires no test modifications
5. **Fast Iteration**: Limited test scope during development enables rapid feedback

## Restoring Full Testing

When 1D implementation is complete, restore full testing by:

1. **Remove test pattern limit** in `.github/workflows/run_tests.yml`:
   ```yaml
   - backend: matlab
     shim_path: ""
     test_pattern: ""  # Change from "kspaceFirstOrder1D" to "" for all tests
   ```

2. **Remove regression test condition**:
   ```yaml
   regression-tests:
     name: Regression tests
     runs-on: ubuntu-latest
     # Remove the "if: github.ref == 'refs/heads/main'" line
   ```

3. **Add 2D/3D matrix entries** when ready:
   ```yaml
   - backend: python-2d
     shim_path: "tests/shims"
     test_pattern: "kspaceFirstOrder2D"
     artifact_suffix: "_python_backend_2d"
   ```
