% Script to validate 1D shim architecture
%
% This script tests whether the shim mechanism correctly redirects
% kspaceFirstOrder1D calls to kspaceFirstOrderPy.
%
% USAGE:
%   Run from repository root:
%   matlab -batch "run('tests/run_shim_validation.m')"

clearvars;

% Determine repository root (script location is tests/)
repo_root = fileparts(fileparts(mfilename('fullpath')));
cd(repo_root);

% Setup Python environment
pyenv('Version', fullfile(repo_root, '.venv310', 'bin', 'python'));

% Add shim path FIRST (to override kspaceFirstOrder1D)
addpath(fullfile(repo_root, 'tests', 'shims'));

% Add k-Wave paths
addpath(fullfile(repo_root, 'k-Wave'));
addpath(fullfile(repo_root, 'k-Wave', 'testing', 'unit'));

fprintf('\n========================================\n');
fprintf('1D Shim Validation Test\n');
fprintf('========================================\n\n');

% Run the test
try
    fprintf('Running kspaceFirstOrder1D_check_source_scaling_p with shim...\n');
    test_result = kspaceFirstOrder1D_check_source_scaling_p(false, false);

    if test_result
        fprintf('\n✓ TEST PASSED\n');
    else
        fprintf('\n✗ TEST FAILED (numerics)\n');
        fprintf('This is expected if Python backend lacks required features.\n');
    end

catch e
    fprintf('\n✗ ERROR: %s\n', e.message);
    fprintf('\nFull error report:\n');
    fprintf('%s\n', getReport(e));
end

fprintf('\n========================================\n');
fprintf('Shim validation complete\n');
fprintf('========================================\n\n');
