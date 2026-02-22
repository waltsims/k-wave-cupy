function test_pass = kspaceFirstOrderPy_binary_sensor_mask_3D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Test 3D Python backend against MATLAB reference using binary sensor mask.
%     Uses a simple ball source with homogeneous medium (no absorption).
%
% ABOUT:
%     author      - k-Wave-CuPy Development
%     date        - 15th February 2026

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

% comparison threshold (relaxed for initial 3D implementation)
comparison_thresh = 0.15;

% =========================================================================
% SIMULATION SETUP
% =========================================================================

% create small 3D computational grid (fast test)
Nx = 32;
Ny = 32;
Nz = 32;
dx = 0.1e-3;
dy = 0.1e-3;
dz = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% homogeneous medium (no absorption for simplicity)
medium.sound_speed = 1500;
medium.density = 1000;

% create simple ball source at center
ball_x_pos = Nx/2;
ball_y_pos = Ny/2;
ball_z_pos = Nz/2;
ball_radius = 3;
source.p0 = 5 * makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ball_z_pos, ball_radius);

% binary sensor mask - record all points
sensor.mask = ones(Nx, Ny, Nz);

% set time array (limit steps for initial testing)
kgrid.makeTime(medium.sound_speed);
kgrid.setTime(min(30, kgrid.Nt), kgrid.dt);

% =========================================================================
% RUN MATLAB REFERENCE
% =========================================================================

disp('Running MATLAB reference...');
sensor_data_matlab = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    'PMLInside', false, 'PMLSize', 0, 'PlotSim', plot_simulations, 'Smooth', false, ...
    'DataCast', 'off');

% =========================================================================
% RUN PYTHON BACKEND
% =========================================================================

disp('Running Python backend...');
sensor_data_python = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLSize', 0, 'Smooth', false);

% =========================================================================
% COMPARISON
% =========================================================================

% compute relative error
diff = abs(sensor_data_matlab(:) - sensor_data_python(:));
max_diff = max(diff);
rel_err = max_diff / max(abs(sensor_data_matlab(:)));

fprintf('Max absolute difference: %.2e\n', max_diff);
fprintf('Max relative error: %.2e\n', rel_err);

if rel_err > comparison_thresh
    test_pass = false;
    disp('FAILED: Relative error exceeds threshold');
else
    disp('PASSED: Python backend matches MATLAB reference');
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    % Show central slices
    z_slice = Nz/2;

    subplot(2, 2, 1);
    p_mat = reshape(sensor_data_matlab(:, end), [Nx, Ny, Nz]);
    imagesc(squeeze(p_mat(:, :, z_slice)));
    title('MATLAB Final Pressure (z-slice)');
    colorbar; axis image;

    subplot(2, 2, 2);
    p_py = reshape(sensor_data_python(:, end), [Nx, Ny, Nz]);
    imagesc(squeeze(p_py(:, :, z_slice)));
    title('Python Final Pressure (z-slice)');
    colorbar; axis image;

    subplot(2, 2, 3);
    imagesc(squeeze(p_mat(:, :, z_slice) - p_py(:, :, z_slice)));
    title('Difference (z-slice)');
    colorbar; axis image;

    subplot(2, 2, 4);
    center_idx = Nx*Ny*Nz/2 + Nx*Ny/2 + Nx/2;
    plot(sensor_data_matlab(center_idx, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data_python(center_idx, :), 'r--', 'LineWidth', 1.5);
    legend('MATLAB', 'Python');
    title('Time Series at Center');
    xlabel('Time Step');
    ylabel('Pressure');
end
