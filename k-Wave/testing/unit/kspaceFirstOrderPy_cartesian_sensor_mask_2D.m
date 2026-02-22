function test_pass = kspaceFirstOrderPy_cartesian_sensor_mask_2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Test Python backend Cartesian sensor mask against MATLAB reference.
%     Passes Cartesian coordinates directly to Python (Delaunay interpolation
%     in Python) and compares against MATLAB kspaceFirstOrder2D with the
%     same Cartesian sensor mask (Delaunay interpolation in MATLAB).
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

% create disc source
source.p0 = 5 * makeDisc(Nx, Ny, Nx/2, Ny/2, 5);

% define Cartesian sensor points on a circle (not on grid points)
sensor_radius = 20 * dx;
num_sensor_points = 50;
angles = linspace(0, 2*pi, num_sensor_points + 1);
angles = angles(1:end-1);
cart_x = sensor_radius * cos(angles);
cart_y = sensor_radius * sin(angles);
sensor.mask = [cart_x; cart_y];

% set time array
kgrid.makeTime(medium.sound_speed);
kgrid.setTime(min(50, kgrid.Nt), kgrid.dt);

% =========================================================================
% RUN MATLAB REFERENCE
% =========================================================================

disp('Running MATLAB reference (Cartesian sensor)...');
sensor_data_matlab = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PMLInside', false, 'PMLSize', 0, 'PlotSim', plot_simulations, 'Smooth', false, ...
    'DataCast', 'off');

% =========================================================================
% RUN PYTHON BACKEND
% =========================================================================

disp('Running Python backend (Cartesian sensor)...');
sensor_data_python = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLSize', 0, 'Smooth', false);

% =========================================================================
% COMPARISON
% =========================================================================

diff_val = abs(sensor_data_matlab(:) - sensor_data_python(:));
max_diff = max(diff_val);
rel_err = max_diff / max(abs(sensor_data_matlab(:)));

fprintf('Max absolute difference: %.2e\n', max_diff);
fprintf('Max relative error: %.2e\n', rel_err);
fprintf('Python sensor_data size: %s\n', mat2str(size(sensor_data_python)));
fprintf('MATLAB sensor_data size: %s\n', mat2str(size(sensor_data_matlab)));

if rel_err > comparison_thresh
    test_pass = false;
    disp('FAILED: Relative error exceeds threshold');
else
    disp('PASSED: Python Cartesian sensor matches MATLAB reference');
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    subplot(2, 2, 1);
    plot(sensor_data_matlab(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data_python(1, :), 'r--', 'LineWidth', 1.5);
    legend('MATLAB', 'Python');
    title('Sensor Point 1 Time Series');
    xlabel('Time Step');
    ylabel('Pressure');

    subplot(2, 2, 2);
    plot(sensor_data_matlab(25, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data_python(25, :), 'r--', 'LineWidth', 1.5);
    legend('MATLAB', 'Python');
    title('Sensor Point 25 Time Series');
    xlabel('Time Step');
    ylabel('Pressure');

    subplot(2, 2, 3);
    imagesc(sensor_data_matlab);
    title('MATLAB Cartesian Sensor Data');
    colorbar;
    xlabel('Time Step');
    ylabel('Sensor Point');

    subplot(2, 2, 4);
    imagesc(abs(sensor_data_matlab - sensor_data_python));
    title('Absolute Difference');
    colorbar;
    xlabel('Time Step');
    ylabel('Sensor Point');
end
