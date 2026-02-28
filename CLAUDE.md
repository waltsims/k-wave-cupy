# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the k-Wave MATLAB Toolbox with an added CuPy-accelerated Python backend. The project implements acoustic wave simulation using k-space pseudospectral methods for medical ultrasound and photoacoustic applications.

**Project Structure:**
- `k-Wave/` - Main MATLAB toolbox (legacy acoustic simulation toolkit)
  - `python/kWavePy.py` - Python/CuPy compute engine implementation
  - `kspaceFirstOrderPy.m` - MATLAB wrapper for Python backend
- `plans/` - Development plan for Python/CuPy integration
- `tests/` - Python backend integration and parity tests
- `pyproject.toml` - Python project configuration (using uv)

## Key Architecture Components

### Core Simulation Functions
The main simulation functions follow the pattern `kspaceFirstOrder[1D|2D|3D][C|G].m`:
- **kspaceFirstOrder1D/2D/3D.m** - Pure MATLAB implementations
- **kspaceFirstOrder2DC/3DC.m** - C++ backend wrappers  
- **kspaceFirstOrder2DG/3DG.m** - CUDA/GPU backend wrappers
- **kspaceFirstOrderAS.m** - Angular spectrum method for time-varying sources

### Key Classes
- **kWaveGrid** - Grid definition and coordinate system management
- **kWaveArray** - Array transducer modeling (sources/sensors)
- **kWaveTransducer** - Transducer simulation and beamforming
- **kWaveDiffusion** - Thermal diffusion modeling

### Numerical Method
The toolbox implements a **k-space pseudospectral method** solving coupled first-order acoustic equations:
- Spatial derivatives: Fourier collocation spectral method
- Temporal scheme: Staggered leapfrog with k-space correction
- Absorption: Fractional Laplacian for power-law absorption
- Boundaries: Split-field Perfectly Matched Layer (PML)

## Development Commands

### Testing Commands
```matlab
% Run all unit tests
cd('k-Wave/testing/unit')
runUnitTests()

% Run specific test pattern
runUnitTests('kspace*')

% Run regression tests 
cd('k-Wave/testing/regression')
runRegressionTests()

% Run comprehensive test suite
cd('k-Wave/testing/kWaveTester')
kWaveTester(struct('force_plot_off', true))
```

### Writing Unit Tests

Unit tests in `k-Wave/testing/unit/` follow a standard format for consistency:

**Function Signature:**
```matlab
function test_pass = functionName_test_description(plot_comparisons, plot_simulations)
```

**Standard Structure:**
```matlab
function test_pass = kspaceFirstOrder2D_my_test(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Brief description of what this test validates.
%
% ABOUT:
%     author      - Your Name
%     date        - 1st January 2026
%     last update - 1st January 2026
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

% [Standard k-Wave license block]

% check for plot inputs, and set to false if nargin is zero
if nargin == 0
    plot_comparisons = false;
    plot_simulations = false;
end

% set pass variable
test_pass = true;

% =========================================================================
% TEST LOGIC
% =========================================================================

% ... test code here ...

% check conditions and set test_pass = false if any fail
if (some_error_condition)
    test_pass = false;
end

% plot comparison results if requested
if plot_comparisons
    figure;
    % ... plotting code ...
end
```

**Key Conventions:**
- Function name format: `functionUnderTest_description_of_test.m`
- Return boolean `test_pass` (true = pass, false = fail)
- Two optional plotting arguments: `plot_comparisons` and `plot_simulations`
- If `nargin == 0`, set both plotting flags to `false` (or `true` for development)
- Use `%#ok<*NOPRT>` to suppress output warnings if needed
- Set `test_pass = true` at start, set to `false` only when checks fail
- For comparison tests, define a `comparison_thresh` variable (e.g., `1e-14`)

**Python Backend Tests:**
For tests involving `kspaceFirstOrderPy`, add Python availability check:
```matlab
% check Python availability
try
    env = pyenv;
    if strcmp(env.Status, "NotLoaded")
        py.list();
    end
catch
    disp('Python not available. Skipping test.');
    return;
end
```

**Example Test Names:**
- `kspaceFirstOrder2D_compare_plane_waves.m`
- `kspaceFirstOrderPy_sensor_masks.m`
- `kspaceFirstOrderPy_record_options.m`

### Running Examples

Example scripts are located in `k-Wave/examples/` and demonstrate various features of the toolbox. Examples can be run in two ways:

#### Method 1: Direct execution (recommended for interactive use)
```matlab
% Run example directly
addpath('k-Wave');
run('k-Wave/examples/example_ivp_homogeneous_medium_py.m');
```

#### Method 2: Via runtests (for automated testing)
```bash
# Run example as a test (requires path setup in the example)
matlab -batch "runtests('k-Wave/examples/example_ivp_homogeneous_medium_py.m')"
```

**Important:** Examples that include path setup (`addpath(fullfile(fileparts(mfilename('fullpath')), '..'));`) can be run via `runtests()` without errors. Examples without path setup will fail with `Undefined function` errors when run via `runtests()`.

**Fix for path-related errors:** Add the following line after `clearvars;` in the example:
```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
```

This ensures the k-Wave toolbox is accessible regardless of how the example is invoked.

### Current Development Status

The project is implementing a **minimalistic Python/CuPy backend** following the plan in `plans/plan.md`. The goal is to create:

1. **Python Compute Engine**: Pure NumPy/CuPy implementation of the k-space PSTD solver
2. **MATLAB Integration**: Direct in-memory calls via `kspaceFirstOrderPy.m`
3. **Testing**: Validation against existing MATLAB test suite

### Development Progress
- **Phase 1** ✅ **COMPLETE**: Data interop verification and 1D implementation achieved <1e-15 parity
- **Phase 1.5** ✅ **COMPLETE**: Code refactoring to ~65 lines (Python) and ~40 lines (MATLAB)
- **Phase 2** 🚧 **IN PROGRESS**: 2D/3D generalization with PML boundaries
- **Phase 3** ✅ **READY**: CuPy GPU acceleration infrastructure implemented

### Python Backend Current Limitations

The Python/CuPy backend (`kspaceFirstOrderPy`) is under active development and has the following limitations:

**Sensor Masks:**
- ✅ **Supported**: Binary grid masks (e.g., `sensor.mask = ones(Nx, Ny)` or `makeCircle()`)
- ✅ **Supported**: Cartesian coordinate masks (e.g., `makeCartCircle()` output, Delaunay interpolation)
- ✅ **Supported**: Empty sensor (defaults to full-grid recording)

**Other Limitations:**
- Advanced sensor types (directional, frequency response) not implemented

When creating examples or tests for the Python backend, ensure they only use supported features.

### File Patterns and Conventions

- **Examples**: `k-Wave/examples/example_*.m` - Documented usage examples
- **Unit Tests**: `k-Wave/testing/unit/kspaceFirstOrder*_*.m` - Individual feature tests
- **Helper Functions**: `k-Wave/private/` - Internal implementation details
- **Utilities**: Grid creation (`make*.m`), post-processing (`*Plot.m`, `*Filter.m`)

### Important Notes

- This is an **acoustic simulation** toolkit - only work on defensive security tools
- The codebase follows MATLAB conventions (camelCase, function documentation)
- Test functions expect signature: `pass = testFunction(plot_simulations, plot_comparisons)`
- Grid indexing: (x,y,z) with MATLAB's column-major order
- Time stepping uses staggered grids (pressure and velocity at different temporal positions)

### Integration Points

When working on the Python backend:
- Data exchange happens through MATLAB's Python interface
- Key concern: Row-major (Python) vs column-major (MATLAB) array indexing
- Validation against `acousticFieldPropagator` tests for correctness
- Performance benchmarking against existing C++/CUDA codes

### MATLAB + Python Setup (arm64)

- Supported Python for R2024b: 3.9, 3.10, 3.11, 3.12 (arm64). Do **not** use Intel builds.
- Project venv: `.venv310` (Python 3.10.16, NumPy 2.2.6) created via `uv venv --python 3.10 .venv310` and `uv pip install --python .venv310/bin/python numpy`.
- **Architecture Requirement**: On Apple Silicon systems, MATLAB must run with ARM64 architecture. If your terminal is emulating x86_64 (Rosetta), prefix commands with `arch -arm64`:
  ```bash
  arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "run('script.m')"
  ```
- If MATLAB reports the environment is already loaded, restart MATLAB and rerun the command.

#### Run all Python-specific tests (no shims needed)
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "
  pyenv('Version', fullfile(pwd,'.venv310','bin','python'));
  addpath('k-Wave');
  cd('k-Wave/testing/unit');
  runUnitTests('kspaceFirstOrderPy')"
```

#### Generate 2D/3D reference data (MATLAB backend, no shims)
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "
  addpath('k-Wave');
  cd('k-Wave/testing/unit');
  kspaceFirstOrder2D_compare_plane_waves(false, false)"
```

#### Run 1D/2D/3D parity tests (Python vs MATLAB reference, via shims)
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "
  pyenv('Version', fullfile(pwd,'.venv310','bin','python'));
  addpath('k-Wave'); addpath('tests/shims');
  cd('k-Wave/testing/unit');
  kspaceFirstOrder1D_compare_plane_waves(false, false)"
```
Replace `1D` with `2D` or `3D` for other dimensions.

#### Run orchestrator (all Py tests or per-dimension parity)
```bash
arch -arm64 /Applications/MATLAB_R2024b.app/bin/matlab -batch "
  pyenv('Version', fullfile(pwd,'.venv310','bin','python'));
  addpath('k-Wave'); addpath('tests');
  run_py_tests()"
```
Or `run_py_tests(1)` / `run_py_tests(2)` / `run_py_tests(3)` for per-dimension parity.

- To visualise differences, run `tests/plot_1d_parity.m` (saves `tests/plots/1d_parity.png`).

### Project Management

- **Dependencies**: Managed via `uv` package manager (see `pyproject.toml`)
- **Current Focus**: Implementing 2D/3D support and PML boundaries (see `plans/plan.md`)
- **Code Philosophy**: Minimalistic "academic code golf" approach - the code is the documentation
- **Performance Target**: Match or exceed C++/CUDA implementations while maintaining interpretability

## Bug workflow (mandatory)

When the user reports a bug, do **not** start by proposing fixes.

### 1) Reproduce-first gate
Before any fix attempt:
- Create or identify a **minimal reproducible test** that fails and demonstrates the bug.
- Prefer automated tests (unit/integration/regression). If automation is not feasible, write a deterministic repro script and expected output.
- Record: environment assumptions, exact command(s) to run, and expected vs actual behavior.

**Hard rule:** No code changes are allowed until a failing repro exists.

### 2) Subagent fix attempts (after repro exists)
Once a failing test exists, spawn subagents to attempt fixes in parallel:
- Each subagent must propose a fix **and** explain how it will make the repro test pass.
- Each subagent must provide either:
  - a patch + evidence the test would pass, or
  - a reasoning-based proof tied directly to the failing test and change.

### 3) Merge criteria
Accept a fix only if:
- The original repro test now passes.
- No related tests regress (or at least, the agent explains any new failures).
- The repro test is kept as a **regression test** (do not delete it).

### 4) Output format when a bug is reported
Always respond in this order:
1. “Repro/Test Plan” (what test you’ll add/run)
2. “Current hypothesis” (optional, but no fixes yet)
3. “Subagent tasks” (what each subagent will try after repro is in place)
4. “Fix + Proof” (only after the failing test exists)
