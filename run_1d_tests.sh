#!/bin/bash
# Fast local test runner for 1D k-Wave Python backend

set -e

MATLAB="/Applications/MATLAB_R2024b.app/bin/matlab"

echo "Running 1D unit tests with Python backend..."
echo ""

arch -arm64 "$MATLAB" -nojvm -batch "pyenv('Version', fullfile(pwd,'.venv310','bin','python')); addpath('k-Wave'); addpath('k-Wave/testing'); addpath('tests/shims'); cd('k-Wave/testing/unit'); test_struct = runUnitTests('kspaceFirstOrder1D', false); show_test_results(test_struct);"
