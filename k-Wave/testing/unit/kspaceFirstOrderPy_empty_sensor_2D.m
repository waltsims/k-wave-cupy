function test_pass = kspaceFirstOrderPy_empty_sensor_2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Test that empty sensor (sensor = []) defaults to recording pressure
%     at all grid points, matching the result from an explicit full-grid
%     binary sensor mask.
%
% ABOUT:
%     author      - k-Wave-CuPy Development
%     date        - 22nd February 2026

% check for plot inputs
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

comparison_thresh = 1e-10;

% =========================================================================
% SIMULATION SETUP
% =========================================================================

% create small computational grid
Nx = 64;
Ny = 64;
dx = 0.1e-3;
dy = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% homogeneous medium
medium.sound_speed = 1500;
medium.density = 1000;

% disc initial pressure source
source.p0 = 5 * makeDisc(Nx, Ny, Nx/2, Ny/2, 5);

% short simulation
kgrid.makeTime(medium.sound_speed);
kgrid.setTime(min(50, kgrid.Nt), kgrid.dt);

% =========================================================================
% RUN WITH EMPTY SENSOR
% =========================================================================

disp('Running Python backend with empty sensor...');
sensor_data_empty = kspaceFirstOrderPy(kgrid, medium, source, [], 'PMLSize', 0, 'Smooth', false);

% check that output is not empty
if isempty(sensor_data_empty)
    test_pass = false;
    disp('FAILED: Empty sensor returned empty output');
    return;
end

% check expected size: (Nx*Ny) x Nt
expected_size = [Nx * Ny, kgrid.Nt];
if ~isequal(size(sensor_data_empty), expected_size)
    test_pass = false;
    fprintf('FAILED: Expected size [%d %d], got [%d %d]\n', ...
        expected_size(1), expected_size(2), size(sensor_data_empty, 1), size(sensor_data_empty, 2));
    return;
end

% =========================================================================
% RUN WITH EXPLICIT FULL-GRID SENSOR
% =========================================================================

disp('Running Python backend with explicit full-grid sensor...');
sensor.mask = ones(Nx, Ny);
sensor_data_full = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLSize', 0, 'Smooth', false);

% =========================================================================
% COMPARISON
% =========================================================================

diff = abs(sensor_data_empty(:) - sensor_data_full(:));
max_diff = max(diff);
rel_err = max_diff / max(abs(sensor_data_full(:)));

fprintf('Max absolute difference: %.2e\n', max_diff);
fprintf('Max relative error: %.2e\n', rel_err);

if rel_err > comparison_thresh
    test_pass = false;
    disp('FAILED: Empty sensor result differs from explicit full-grid sensor');
else
    disp('PASSED: Empty sensor defaults to full-grid recording');
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    subplot(1, 3, 1);
    imagesc(reshape(sensor_data_empty(:, end), [Nx, Ny]));
    title('Empty Sensor');
    colorbar; axis image;

    subplot(1, 3, 2);
    imagesc(reshape(sensor_data_full(:, end), [Nx, Ny]));
    title('Full-Grid Sensor');
    colorbar; axis image;

    subplot(1, 3, 3);
    imagesc(reshape(sensor_data_empty(:, end) - sensor_data_full(:, end), [Nx, Ny]));
    title('Difference');
    colorbar; axis image;
end
