function test_pass = kspaceFirstOrderPy_sensor_masks(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Consolidated test for binary, Cartesian, and empty sensor masks.
%     Compares Python backend against MATLAB reference for binary and
%     Cartesian masks, and verifies empty sensor self-consistency.
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
% SHARED SIMULATION SETUP
% =========================================================================

Nx = 64;
Ny = 64;
dx = 0.1e-3;
dy = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = 5 * makeDisc(Nx, Ny, Nx/2, Ny/2, 5);

kgrid.makeTime(medium.sound_speed);
kgrid.setTime(min(50, kgrid.Nt), kgrid.dt);

comparison_thresh = 1e-10;

% =========================================================================
% SUB-TEST 1: BINARY SENSOR MASK
% =========================================================================

disp('--- Sub-test 1: Binary sensor mask ---');

sensor_binary.mask = ones(Nx, Ny);

% MATLAB reference
disp('Running MATLAB reference (binary mask)...');
sensor_data_matlab_binary = kspaceFirstOrder2D(kgrid, medium, source, sensor_binary, ...
    'PMLInside', false, 'PMLSize', 0, 'PlotSim', plot_simulations, 'Smooth', false, ...
    'DataCast', 'off');

% Python backend
disp('Running Python backend (binary mask)...');
sensor_data_python_binary = kspaceFirstOrderPy(kgrid, medium, source, sensor_binary, ...
    'PMLSize', 0, 'Smooth', false);

% compare
diff_binary = abs(sensor_data_matlab_binary(:) - sensor_data_python_binary(:));
max_diff_binary = max(diff_binary);
rel_err_binary = max_diff_binary / max(abs(sensor_data_matlab_binary(:)));

fprintf('  Max absolute difference: %.2e\n', max_diff_binary);
fprintf('  Max relative error:      %.2e\n', rel_err_binary);

if rel_err_binary > comparison_thresh
    test_pass = false;
    disp('FAILED: Binary mask - relative error exceeds threshold');
else
    disp('PASSED: Binary mask - Python matches MATLAB reference');
end

% =========================================================================
% SUB-TEST 2: CARTESIAN SENSOR MASK (self-consistency)
% =========================================================================

disp('--- Sub-test 2: Cartesian sensor mask ---');

sensor_radius = 20 * dx;
num_sensor_points = 50;
angles = linspace(0, 2*pi, num_sensor_points + 1);
angles = angles(1:end-1);
sensor_cart.mask = [sensor_radius * cos(angles); sensor_radius * sin(angles)];

% Python backend with Cartesian mask
disp('Running Python backend (Cartesian mask)...');
sensor_data_python_cart = kspaceFirstOrderPy(kgrid, medium, source, sensor_cart, ...
    'PMLSize', 0, 'Smooth', false);

% Manual bilinear interpolation of Python binary output at Cartesian points
x_vec = ((0:Nx-1) - Nx/2) * dx;
y_vec = ((0:Ny-1) - Ny/2) * dy;
manual_cart = zeros(num_sensor_points, kgrid.Nt);
for t = 1:kgrid.Nt
    p_field = reshape(sensor_data_python_binary(:, t), [Nx, Ny]);
    F = griddedInterpolant({x_vec, y_vec}, p_field, 'linear');
    manual_cart(:, t) = F(sensor_cart.mask(1,:)', sensor_cart.mask(2,:)');
end

% compare Python Cartesian output against bilinear reference
diff_cart = abs(sensor_data_python_cart(:) - manual_cart(:));
max_diff_cart = max(diff_cart);
rel_err_cart = max_diff_cart / max(abs(manual_cart(:)));

fprintf('  Max absolute difference: %.2e\n', max_diff_cart);
fprintf('  Max relative error:      %.2e\n', rel_err_cart);
fprintf('  Python sensor_data size: %s\n', mat2str(size(sensor_data_python_cart)));

if rel_err_cart > comparison_thresh
    test_pass = false;
    disp('FAILED: Cartesian mask - Python bilinear interp does not match reference');
else
    disp('PASSED: Cartesian mask - Python bilinear interp matches reference');
end

% =========================================================================
% SUB-TEST 3: EMPTY SENSOR (SELF-CONSISTENCY)
% =========================================================================

disp('--- Sub-test 3: Empty sensor ---');

% run with empty sensor
disp('Running Python backend with empty sensor...');
sensor_data_empty = kspaceFirstOrderPy(kgrid, medium, source, [], 'PMLSize', 0, 'Smooth', false);

% check output is not empty
if isempty(sensor_data_empty)
    test_pass = false;
    disp('FAILED: Empty sensor returned empty output');
else
    % check expected size: (Nx*Ny) x Nt
    expected_size = [Nx * Ny, kgrid.Nt];
    if ~isequal(size(sensor_data_empty), expected_size)
        test_pass = false;
        fprintf('FAILED: Expected size [%d %d], got [%d %d]\n', ...
            expected_size(1), expected_size(2), ...
            size(sensor_data_empty, 1), size(sensor_data_empty, 2));
    else
        % run with explicit full-grid sensor
        disp('Running Python backend with explicit full-grid sensor...');
        sensor_full.mask = ones(Nx, Ny);
        sensor_data_full = kspaceFirstOrderPy(kgrid, medium, source, sensor_full, ...
            'PMLSize', 0, 'Smooth', false);

        % compare
        diff_empty = abs(sensor_data_empty(:) - sensor_data_full(:));
        max_diff_empty = max(diff_empty);
        rel_err_empty = max_diff_empty / max(abs(sensor_data_full(:)));

        fprintf('  Max absolute difference: %.2e\n', max_diff_empty);
        fprintf('  Max relative error:      %.2e\n', rel_err_empty);

        if rel_err_empty > comparison_thresh
            test_pass = false;
            disp('FAILED: Empty sensor result differs from explicit full-grid sensor');
        else
            disp('PASSED: Empty sensor defaults to full-grid recording');
        end
    end
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    % binary mask: MATLAB / Python / difference
    subplot(2, 3, 1);
    imagesc(reshape(sensor_data_matlab_binary(:, end), [Nx, Ny]));
    title('Binary: MATLAB');
    colorbar; axis image;

    subplot(2, 3, 2);
    imagesc(reshape(sensor_data_python_binary(:, end), [Nx, Ny]));
    title('Binary: Python');
    colorbar; axis image;

    subplot(2, 3, 3);
    imagesc(reshape(sensor_data_matlab_binary(:, end) - sensor_data_python_binary(:, end), [Nx, Ny]));
    title('Binary: Difference');
    colorbar; axis image;

    % Cartesian mask: MATLAB / Python / difference
    subplot(2, 3, 4);
    plot(sensor_data_matlab_cart(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data_python_cart(1, :), 'r--', 'LineWidth', 1.5);
    legend('MATLAB', 'Python');
    title('Cartesian: Sensor 1');
    xlabel('Time Step');
    ylabel('Pressure');

    subplot(2, 3, 5);
    plot(sensor_data_matlab_cart(25, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data_python_cart(25, :), 'r--', 'LineWidth', 1.5);
    legend('MATLAB', 'Python');
    title('Cartesian: Sensor 25');
    xlabel('Time Step');
    ylabel('Pressure');

    subplot(2, 3, 6);
    imagesc(abs(sensor_data_matlab_cart - sensor_data_python_cart));
    title('Cartesian: Abs Difference');
    colorbar;
    xlabel('Time Step');
    ylabel('Sensor Point');
end
