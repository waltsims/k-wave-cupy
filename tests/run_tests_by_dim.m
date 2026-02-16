function results = run_tests_by_dim(dim, verbose)
%RUN_TESTS_BY_DIM Run k-Wave tests for a specific dimension using Python backend.
%
% USAGE:
%   results = run_tests_by_dim(1)       % Run all 1D tests
%   results = run_tests_by_dim(2)       % Run all 2D tests
%   results = run_tests_by_dim(3)       % Run all 3D tests
%   results = run_tests_by_dim(1, true) % Verbose output

if nargin < 2, verbose = false; end

% Setup paths
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'k-Wave'));
addpath(fullfile(project_root, 'tests', 'shims'));

% Setup Python (use PYTHON_PATH env var if set, otherwise local venv)
python_path = getenv('PYTHON_PATH');
if isempty(python_path)
    python_path = fullfile(project_root, '.venv310', 'bin', 'python');
end
pyenv('Version', python_path);

% Get test directory
test_dir = fullfile(project_root, 'k-Wave', 'testing', 'unit');

% Define tests by dimension
tests_1d = {
    'kspaceFirstOrder1D_compare_plane_waves'
    'kspaceFirstOrder1D_check_source_scaling_p'
    'kspaceFirstOrder1D_check_source_scaling_ux'
    'kspaceFirstOrder1D_p0_vs_two_clicks'
};

tests_2d = {
    'kspaceFirstOrderPy_binary_sensor_mask_2D'
    'kspaceFirstOrder2D_check_source_scaling_p'
    'kspaceFirstOrder2D_check_source_scaling_ux'
    'kspaceFirstOrder2D_check_source_scaling_uy'
    'kspaceFirstOrder2D_p0_vs_two_clicks'
};

tests_3d = {
    'kspaceFirstOrderPy_binary_sensor_mask_3D'
    'kspaceFirstOrder3D_check_source_scaling_p'
    'kspaceFirstOrder3D_check_source_scaling_ux'
    'kspaceFirstOrder3D_check_source_scaling_uy'
    'kspaceFirstOrder3D_check_source_scaling_uz'
    'kspaceFirstOrder3D_p0_vs_two_clicks'
};

% Select tests based on dimension
switch dim
    case 1
        tests = tests_1d;
        dim_name = '1D';
    case 2
        tests = tests_2d;
        dim_name = '2D';
    case 3
        tests = tests_3d;
        dim_name = '3D';
    otherwise
        error('Dimension must be 1, 2, or 3');
end

% Run tests
fprintf('\n=== Running %s Tests (Python Backend) ===\n\n', dim_name);
n_tests = numel(tests);
results = struct('name', {}, 'passed', {}, 'error', {});
for j = 1:n_tests
    results(j).name = tests{j};
    results(j).passed = false;
    results(j).error = '';
end

cd(test_dir);
for i = 1:n_tests
    test_name = tests{i};
    fprintf('[%d/%d] %s... ', i, n_tests, test_name);

    try
        pass = feval(test_name, false, false);
        results(i).passed = pass;
        if pass
            fprintf('PASSED\n');
        else
            fprintf('FAILED\n');
        end
    catch ME
        results(i).passed = false;
        results(i).error = ME.message;
        fprintf('ERROR: %s\n', ME.message);
        if verbose
            disp(getReport(ME));
        end
    end
end

% Summary
n_passed = sum([results.passed]);
fprintf('\n=== %s Summary: %d/%d tests passed ===\n', dim_name, n_passed, n_tests);

if n_passed < n_tests
    fprintf('\nFailed tests:\n');
    for i = 1:n_tests
        if ~results(i).passed
            fprintf('  - %s', results(i).name);
            if ~isempty(results(i).error)
                fprintf(' (Error: %s)', results(i).error);
            end
            fprintf('\n');
        end
    end
end

end
