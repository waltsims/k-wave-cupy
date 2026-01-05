# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the k-Wave MATLAB Toolbox with an added CuPy-accelerated Python backend. The project implements acoustic wave simulation using k-space pseudospectral methods for medical ultrasound and photoacoustic applications.

**Project Structure:**
- `k-Wave/` - Main MATLAB toolbox (legacy acoustic simulation toolkit)
- `plans/` - Development plan for Python/CuPy integration
- `tests/` - Unit and regression tests (currently empty, uses k-Wave's internal test suite)

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

### Current Development Focus

The project is implementing a **minimalistic Python/CuPy backend** following the plan in `plans/plan.md`. The goal is to create:

1. **Python Compute Engine**: Pure NumPy/CuPy implementation of the k-space PSTD solver
2. **MATLAB Integration**: Direct in-memory calls via `kspaceFirstOrderPy.m`
3. **Testing**: Validation against existing MATLAB test suite

### Development Phases
1. **Phase 1**: Data interop verification and 1D implementation
2. **Phase 2**: 2D/3D generalization with PML boundaries  
3. **Phase 3**: CuPy GPU acceleration

### File Patterns and Conventions

- **Examples**: `k-Wave/examples/example_*.m` - Documented usage examples
- **Unit Tests**: `k-Wave/testing/unit/*_test*.m` - Individual feature tests
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
- Run tests with MATLAB using this venv:
  ```
  "/Applications/MATLAB_R2024b.app/bin/matlab" -batch "pyenv('Version', fullfile(pwd,'.venv310','bin','python')); addpath('k-Wave'); addpath('tests'); addpath('k-Wave/testing/unit'); runtests({'tests/test_interop_sanity','k-Wave/testing/unit/test_interface_1D'});"
  ```
- If MATLAB reports the environment is already loaded, restart MATLAB and rerun the command.
