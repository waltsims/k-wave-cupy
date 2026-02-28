function results = run_py_tests(dim)
%RUN_PY_TESTS Run Python backend tests.
%
% USAGE:
%   run_py_tests()      Run all kspaceFirstOrderPy_* unit tests
%   run_py_tests(dim)   Run kspaceFirstOrder{dim}D_compare_plane_waves + shimmed upstream tests
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

% Setup paths
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'k-Wave'));
test_dir = fullfile(project_root, 'k-Wave', 'testing', 'unit');

% Setup Python
python_path = getenv('PYTHON_PATH');
if isempty(python_path)
    python_path = fullfile(project_root, '.venv310', 'bin', 'python');
end
pyenv('Version', python_path);

if nargin == 0
    % Mode 1: Run all kspaceFirstOrderPy_* tests (no shims needed)
    files = dir(fullfile(test_dir, 'kspaceFirstOrderPy_*.m'));
    tests = cellfun(@(f) f(1:end-2), {files.name}, 'UniformOutput', false);
    dim_name = 'Py';
else
    % Mode 2: Run dimension-specific parity tests (needs shims)
    addpath(fullfile(project_root, 'tests', 'shims'));

    % Auto-discover shimmed tests for this dimension
    pattern = sprintf('kspaceFirstOrder%dD_', dim);
    files = dir(fullfile(test_dir, [pattern '*.m']));
    tests = cellfun(@(f) f(1:end-2), {files.name}, 'UniformOutput', false);

    % Exclude reference data .mat files (just in case)
    tests = tests(~contains(tests, 'reference_data'));

    dim_name = sprintf('%dD', dim);
end

% Run tests
fprintf('\n=== Running %s Tests (Python Backend) ===\n\n', dim_name);
n_tests = numel(tests);
results = struct('name', {}, 'passed', {}, 'error', {});

cd(test_dir);
for i = 1:n_tests
    test_name = tests{i};
    results(i).name = test_name;
    results(i).passed = false;
    results(i).error = '';

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
