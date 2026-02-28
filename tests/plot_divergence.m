% PLOT_DIVERGENCE  Compare MATLAB and Python solver outputs over time.
%
% Runs identical simulations with full-grid sensors, recording pressure at
% every grid point and every timestep.  Plots:
%   1. Total acoustic energy (sum p^2) vs timestep
%   2. Energy difference (MATLAB - Python) vs timestep
%   3. Max |p| vs timestep
%   4. Relative error (max |diff| / max |p_matlab|) vs timestep
%
% For lossless media the total energy should remain approximately constant
% (PML absorbs some), so this reveals whether either solver has a spurious
% energy gain or loss relative to the other.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'k-Wave'));

fprintf('\n=== Solver Divergence Plots ===\n');

%% 1D: p0 source, lossless, 600 steps
fprintf('\nRunning 1D comparison (Nx=64, Nt=600) ...\n');
[pm_1d, pp_1d] = run_1d(64, 600);

%% 2D: p0 source, lossless, 200 steps (longer run to see divergence growth)
fprintf('\nRunning 2D comparison (64x64, Nt=200) ...\n');
[pm_2d, pp_2d] = run_2d(64, 64, 200);

%% 2D: time-varying source, lossless, 200 steps
fprintf('\nRunning 2D comparison with source.p (64x64, Nt=200) ...\n');
[pm_2d_src, pp_2d_src] = run_2d_source(64, 64, 200);

%% Plot
plot_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plot_dir, 'dir'), mkdir(plot_dir); end

plot_comparison(pm_1d, pp_1d, '1D (Nx=64, p0, lossless)', ...
    fullfile(plot_dir, 'divergence_1d.png'));
plot_comparison(pm_2d, pp_2d, '2D (64x64, p0, lossless)', ...
    fullfile(plot_dir, 'divergence_2d.png'));
plot_comparison(pm_2d_src, pp_2d_src, '2D (64x64, source.p, lossless)', ...
    fullfile(plot_dir, 'divergence_2d_source.png'));

fprintf('\nPlots saved to tests/plots/\n');

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [p_matlab, p_python] = run_1d(Nx, Nt)
    dx = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    source.p0 = zeros(Nx, 1);
    source.p0(Nx/4) = 1;

    % Full-grid sensor: record pressure at every grid point
    sensor.mask = ones(Nx, 1);
    sensor.record = {'p'};

    % MATLAB native
    sensor_data_m = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...
        'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;  % (Nx, Nt)

    % Python backend
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor);
    p_python = sensor_data_p.p;  % (Nx, Nt)

    fprintf('  1D: p_matlab size = [%s], p_python size = [%s]\n', ...
        num2str(size(p_matlab)), num2str(size(p_python)));
end

function [p_matlab, p_python] = run_2d(Nx, Ny, Nt)
    dx = 1e-4; dy = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    source.p0 = zeros(Nx, Ny);
    source.p0(Nx/2, Ny/2) = 1;

    % Full-grid sensor: record pressure at every grid point
    sensor.mask = ones(Nx, Ny);
    sensor.record = {'p'};

    % MATLAB native
    sensor_data_m = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PMLInside', true, 'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;  % (Nx*Ny, Nt)

    % Python backend
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLInside', true);
    p_python = sensor_data_p.p;  % (Nx*Ny, Nt)

    fprintf('  2D: p_matlab size = [%s], p_python size = [%s]\n', ...
        num2str(size(p_matlab)), num2str(size(p_python)));
end

function [p_matlab, p_python] = run_2d_source(Nx, Ny, Nt)
    dx = 1e-4; dy = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    % Time-varying pressure source
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(Nx/4, Ny/2) = 1;
    source.p = sin(2 * pi * 1e6 * (0:Nt-1) * dt);
    source.p_mode = 'additive-no-correction';

    % Full-grid sensor
    sensor.mask = ones(Nx, Ny);
    sensor.record = {'p'};

    % MATLAB native
    sensor_data_m = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PMLInside', true, 'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;

    % Python backend
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLInside', true);
    p_python = sensor_data_p.p;

    fprintf('  2D src: p_matlab size = [%s], p_python size = [%s]\n', ...
        num2str(size(p_matlab)), num2str(size(p_python)));
end

function plot_comparison(p_matlab, p_python, title_str, save_path)
    Nt = size(p_matlab, 2);
    t = 1:Nt;

    % Acoustic energy: proportional to sum(p^2) over space at each timestep
    energy_m = sum(p_matlab.^2, 1);
    energy_p = sum(p_python.^2, 1);

    % Max |p| at each timestep
    maxp_m = max(abs(p_matlab), [], 1);
    maxp_p = max(abs(p_python), [], 1);

    % Pointwise difference metrics at each timestep
    abs_diff = abs(p_matlab - p_python);
    max_diff = max(abs_diff, [], 1);
    rel_err = max_diff ./ max(maxp_m, eps);

    % Mean pressure (should be ~0 for symmetric p0)
    meanp_m = mean(p_matlab, 1);
    meanp_p = mean(p_python, 1);

    figure('Position', [50 50 1200 900], 'Visible', 'off');

    % --- Row 1: Energy ---
    subplot(3, 2, 1);
    plot(t, energy_m, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, energy_p, 'r--', 'LineWidth', 1.5);
    xlabel('Timestep'); ylabel('sum(p^2)');
    title('Total Acoustic Energy');
    legend('MATLAB', 'Python', 'Location', 'best');
    grid on;

    subplot(3, 2, 2);
    plot(t, energy_m - energy_p, 'k-', 'LineWidth', 1.2);
    xlabel('Timestep'); ylabel('\Delta sum(p^2)');
    title('Energy Difference (MATLAB - Python)');
    grid on;

    % --- Row 2: Max |p| ---
    subplot(3, 2, 3);
    plot(t, maxp_m, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, maxp_p, 'r--', 'LineWidth', 1.5);
    xlabel('Timestep'); ylabel('max |p|');
    title('Peak Pressure');
    legend('MATLAB', 'Python', 'Location', 'best');
    grid on;

    subplot(3, 2, 4);
    semilogy(t, max_diff, 'k-', 'LineWidth', 1.2);
    xlabel('Timestep'); ylabel('max |p_M - p_P|');
    title('Max Absolute Difference');
    grid on;

    % --- Row 3: Relative error & Mean pressure ---
    subplot(3, 2, 5);
    semilogy(t, rel_err, 'k-', 'LineWidth', 1.2);
    xlabel('Timestep'); ylabel('Relative error');
    title('Relative Error (max |diff| / max |p_M|)');
    grid on;

    subplot(3, 2, 6);
    plot(t, meanp_m, 'b-', 'LineWidth', 1.5); hold on;
    plot(t, meanp_p, 'r--', 'LineWidth', 1.5);
    xlabel('Timestep'); ylabel('mean(p)');
    title('Mean Pressure (Domain Average)');
    legend('MATLAB', 'Python', 'Location', 'best');
    grid on;

    sgtitle(sprintf('MATLAB vs Python Divergence: %s', title_str));

    saveas(gcf, save_path);
    fprintf('  Saved: %s\n', save_path);
    close(gcf);
end
