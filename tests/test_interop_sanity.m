function tests = test_interop_sanity
%TEST_INTEROP_SANITY Validate MATLAB <-> Python array layout.
%
% Sends a non-symmetric array to Python, mutates it there, and checks the
% modified value lands in the expected MATLAB position (column-major
% layout).
%
% This is a guard against silent transposes when moving data between
% MATLAB and NumPy/CuPy.
tests = functiontests(localfunctions);
end

function testPythonSeesColumnMajor(testCase)
% ensure we can talk to Python
try
    module = importPythonModule(testCase);
catch ME
    testCase.assumeFail("Python unavailable: " + ME.message);
end

% asymmetric shape to catch row/column swaps
matlab_array = reshape(1:6, [2, 3]);

result = module.interop_sanity(py.numpy.array(matlab_array, pyargs('order', 'F')));
returned = double(result);

expected = matlab_array;
expected(1, 2) = 99; % Python modifies A[0, 1]

testCase.verifyEqual(returned, expected);
end

function module = importPythonModule(testCase)
% Add the repo python folder to sys.path and load kWavePy.
repo_root = fileparts(fileparts(mfilename('fullpath')));
python_path = fullfile(repo_root, 'k-Wave', 'python');

try
    % n.b. py.sys.path returns a list object
    paths = cell(py.sys.path);
    if ~any(strcmp(paths, python_path))
        insert(py.sys.path, int32(0), python_path);
    end
catch ME
    testCase.assumeFail("Failed to configure Python path: " + ME.message);
end

module = py.importlib.import_module('kWavePy');
py.importlib.reload(module);
end
