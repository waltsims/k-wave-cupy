function test_pass = kspaceFirstOrderPy_interop_sanity(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Validate MATLAB <-> Python array layout (column-major interop) and
%     optionally test CuPy GPU backend parity against CPU.
%
%     Sends a non-symmetric array to Python, mutates it there, and checks
%     the modified value lands in the expected MATLAB position. This is a
%     guard against silent transposes when moving data between MATLAB and
%     NumPy/CuPy.
%
% ABOUT:
%     author      - k-Wave-CuPy Development
%     date        - 27th February 2026
%     last update - 27th February 2026
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

% check for plot inputs, and set to false if nargin is zero
if nargin == 0
    plot_comparisons = false;
    plot_simulations = false;
end

% set pass variable
test_pass = true;

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

% =========================================================================
% IMPORT PYTHON MODULE
% =========================================================================

% add k-Wave/python to sys.path and load kWavePy
repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
python_path = fullfile(repo_root, 'k-Wave', 'python');

paths = cell(py.sys.path);
if ~any(strcmp(paths, python_path))
    insert(py.sys.path, int32(0), python_path);
end

module = py.importlib.import_module('kWavePy');
py.importlib.reload(module);

% =========================================================================
% TEST 1: COLUMN-MAJOR INTEROP
% =========================================================================

disp('--- Sub-test 1: Column-major interop ---');

% asymmetric shape to catch row/column swaps
matlab_array = reshape(1:6, [2, 3]);

% send to Python and call interop_sanity (modifies A[0,1])
result = module.interop_sanity(py.numpy.array(matlab_array, pyargs('order', 'F')));
returned = double(result);

expected = matlab_array;
expected(1, 2) = 99;  % Python modifies A[0, 1]

if ~isequal(returned, expected)
    test_pass = false;
    disp('FAILED: Returned array does not match expected (column-major mismatch)');
    disp('Expected:');
    disp(expected);
    disp('Got:');
    disp(returned);
else
    disp('PASSED: Column-major interop verified');
end

% =========================================================================
% TEST 2: CUPY GPU BACKEND (OPTIONAL)
% =========================================================================

disp('--- Sub-test 2: CuPy GPU backend ---');

% CuPy GPU backend check (optional)
try
    py.importlib.import_module('cupy');
    has_cupy = true;
catch
    has_cupy = false;
end

if has_cupy
    disp('CuPy available - testing GPU backend...');

    % run a trivial 1D simulation with both backends
    kgrid = kWaveGrid(32, 1);
    kgrid.setTime(20, 0.1 / 1500);
    medium.sound_speed = 1500;
    medium.density = 1000;
    source.p0 = zeros(32, 1);
    source.p0(10) = 1;
    sensor.mask = zeros(32, 1);
    sensor.mask(20) = 1;

    data_cpu = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLSize', 0, 'Smooth', false);
    data_gpu = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLSize', 0, 'Smooth', false, 'DataCast', 'cupy');

    gpu_err = max(abs(data_cpu(:) - data_gpu(:))) / max(abs(data_cpu(:)));
    if gpu_err > 1e-12
        test_pass = false;
        disp(['FAILED: CuPy vs CPU mismatch: ' num2str(gpu_err, '%.2e')]);
    else
        disp(['PASSED: CuPy matches CPU (err=' num2str(gpu_err, '%.2e') ')']);
    end

    % plot comparison if requested
    if plot_comparisons
        figure;
        plot(data_cpu(:), 'b-', 'LineWidth', 1.5);
        hold on;
        plot(data_gpu(:), 'r--', 'LineWidth', 1.5);
        legend('CPU', 'CuPy GPU');
        title('1D Simulation: CPU vs CuPy');
        xlabel('Time Step');
        ylabel('Pressure');
    end
else
    disp('CuPy not available - skipping GPU backend check.');
end
