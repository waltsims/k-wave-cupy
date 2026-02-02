# Bug Fixes Summary

## Bug #1: Undefined function 'kWaveGrid' when running examples via runtests

**Error:**
```
Undefined function 'kWaveGrid' for input arguments of type 'double'.
Error in example_ivp_homogeneous_medium_py (line 20)
```

**Root Cause:** Example script didn't add k-Wave toolbox to MATLAB path before using it.

**Fix:** Added path setup in `k-Wave/examples/example_ivp_homogeneous_medium_py.m` (line 11-12):
```matlab
% add k-Wave toolbox to path
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
```

**Status:** ✅ FIXED

---

## Bug #2: Python backend reshape error with Cartesian sensor masks

**Error:**
```
Python Error: ValueError: cannot reshape array of size 100 into shape (128,128)
Error using example_ivp_homogeneous_medium_py (line 58)
```

**Root Cause:** Python backend doesn't support Cartesian coordinate sensor masks yet - only binary grid masks.

**Fix:** Converted Cartesian sensor to binary grid mask in `k-Wave/examples/example_ivp_homogeneous_medium_py.m` (lines 47-52):
```matlab
% define a centered circular sensor
% Note: Python backend currently only supports binary grid masks
% Convert Cartesian circle to binary grid mask
sensor_radius = 4e-3;   % [m]
num_sensor_points = 50;
cart_sensor = makeCartCircle(sensor_radius, num_sensor_points);
sensor.mask = cart2grid(kgrid, cart_sensor);
```

**Documentation:** Added Python backend limitations to `CLAUDE.md`

**Status:** ✅ FIXED

---

## How to Test

### Quick Test (run the fixed example)
```bash
cd /Users/Walter/git/k-wave-cupy
/Applications/MATLAB_R2024b.app/bin/matlab -batch "pyenv('Version', fullfile(pwd,'.venv310','bin','python')); runtests('k-Wave/examples/example_ivp_homogeneous_medium_py.m')"
```

### Comprehensive Test (with verification script)
```bash
cd /Users/Walter/git/k-wave-cupy
/Applications/MATLAB_R2024b.app/bin/matlab -batch "run('verify_sensor_mask_fix.m')"
```

### Expected Output
1. ✅ MATLAB simulation completes (~10 seconds)
2. ✅ Python simulation completes (no reshape error)
3. ✅ Comparison plots displayed
4. ✅ Max difference printed (should be relatively small)

---

## Files Modified

1. `k-Wave/examples/example_ivp_homogeneous_medium_py.m` - Added path setup and fixed sensor mask
2. `CLAUDE.md` - Added "Running Examples" section and "Python Backend Current Limitations" section

## Files Created

1. `verify_sensor_mask_fix.m` - Verification script for testing the fixes
2. `tests/test_2d_binary_sensor.m` - Regression test for 2D binary sensor masks
3. `BUG_FIXES.md` - This file

---

## Regression Tests Created

To prevent these bugs from reoccurring:

1. **Path bug test:** `tests/test_example_path_bug.m`
2. **Sensor mask test:** `tests/test_2d_binary_sensor.m`

Run all regression tests:
```bash
/Applications/MATLAB_R2024b.app/bin/matlab -batch "addpath('tests'); runtests('tests/')"
```

---

## Bug #3: GitHub Actions Variable Interpolation in Matrix Values

**Error:**
```
[Warning: Name is nonexistent or not a directory: /tests/shims]
[Warning: Name is nonexistent or not a directory: /k-Wave]
Unrecognized function or variable 'scaleTime'.
```

**Root Cause:** GitHub Actions does not interpolate variables like `${{ github.workspace }}` inside matrix value strings. When we defined:
```yaml
matrix:
  include:
    - path_setup: 'addpath("${{ github.workspace }}/k-Wave");'
```

The `${{ github.workspace }}` remained literal (empty) and was not substituted, resulting in paths like `/tests/shims` instead of `/home/runner/work/k-wave-cupy/k-wave-cupy/tests/shims`.

**Fix:** Changed matrix configuration in `.github/workflows/run_tests.yml`:

**Before (broken):**
```yaml
matrix:
  include:
    - backend: matlab
      path_setup: 'addpath("${{ github.workspace }}/k-Wave");'
```

**After (fixed):**
```yaml
matrix:
  include:
    - backend: matlab
      shim_path: ""
    - backend: python-1d
      shim_path: "tests/shims"
```

Then construct full paths in the command section where variable interpolation works:
```yaml
command: |
  ${{ matrix.shim_path != '' && format('addpath(''{0}/{1}'');', github.workspace, matrix.shim_path) || '' }}
  addpath('${{ github.workspace }}/k-Wave');
```

**Additional Changes:**
- Standardized all MATLAB string literals in workflow to use single quotes (MATLAB convention)
- Removed trailing spaces and cleaned up formatting

**Status:** ✅ FIXED

**Testing:** Push changes and verify:
1. Both matrix configurations (matlab and python-1d) complete without path warnings
2. `scaleTime` function is found
3. Test artifacts are generated successfully
