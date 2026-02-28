% PLOT_DIVERGENCE_SCENARIOS  Full-grid divergence plots for specific 2D scenarios.
%
% Runs selected scenarios from kspaceFirstOrder2D_compare_plane_waves with
% a full-grid sensor to capture pressure at every grid point and timestep.
% Produces energy, peak pressure, and error growth plots to visualize
% exactly how MATLAB and Python diverge.
%
% Selected scenarios:
%   #80: nonlinear + stokes + source.ux (additive) + heterogeneous  (worst: 2.1e-13)
%   #40: linear + stokes + source.ux (add-no-corr) + heterogeneous  (1.8e-13)
%   #24: linear + lossy + source.ux (additive) + heterogeneous      (1.7e-13)
%   #1:  linear + lossless + source.p0 + homogeneous                (best: 1.0e-14)

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'k-Wave'));

fprintf('\n=== 2D Scenario Divergence Plots ===\n');

% =========================================================================
% SIMULATION PARAMETERS (matching compare_plane_waves)
% =========================================================================

Nx = 32; dx = 1; Ny = 32; dy = 1;
PML_size = 10; PML_alpha_def = 2;
c0 = 1500; rho0 = 1000; BonA0 = 10; alpha0 = 5; y = 1.5;
c1 = 2000; rho1 = 1200; BonA1 = 5; alpha1 = 2;
interface_position = Nx / 2;
cfl = 0.1; Nt = 50;
dt = cfl * dx / c0;
t_array = 0:dt:(Nt - 1) * dt;
source_strength = 5e6;
source_position_x = PML_size - 6;
source_position_y = Ny / 2;
source_freq = 200;
source_signal = source_strength * sin(2 * pi * source_freq * t_array);

input_args_base = {'PMLSize', PML_size, 'Smooth', false, ...
    'UsekSpace', true, 'UseSG', true, 'PlotSim', false};

% =========================================================================
% RUN SCENARIOS
% =========================================================================

scenarios = {
    1,  'linear + lossless + p0 + homo (best: ~1e-14)';
    24, 'linear + lossy + ux(add) + hetero (~1.7e-13)';
    40, 'linear + stokes + ux(add-nc) + hetero (~1.8e-13)';
    80, 'nonlinear + stokes + ux(add) + hetero (worst: ~2.1e-13)';
};

plot_dir = fullfile(fileparts(mfilename('fullpath')), 'plots');
if ~exist(plot_dir, 'dir'), mkdir(plot_dir); end

all_pm = cell(size(scenarios, 1), 1);
all_pp = cell(size(scenarios, 1), 1);

for s = 1:size(scenarios, 1)
    test_num = scenarios{s, 1};
    label = scenarios{s, 2};
    fprintf('\n--- Scenario %d: %s ---\n', test_num, label);

    [pm, pp] = run_scenario(test_num);
    all_pm{s} = pm;
    all_pp{s} = pp;
    fname = sprintf('divergence_2d_scenario_%02d.png', test_num);
    plot_comparison(pm, pp, sprintf('Scenario %d: %s', test_num, label), ...
        fullfile(plot_dir, fname));
end

% Combined error-growth comparison across all scenarios
plot_error_comparison(scenarios, all_pm, all_pp, plot_dir);

fprintf('\nAll plots saved to tests/plots/\n');

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [p_matlab, p_python] = run_scenario(test_num)
    % Access outer workspace variables via nested function scoping
    % (MATLAB requires assignin/evalin workaround for scripts, but since
    % this is a script we just use the workspace directly)

    evalin('caller', 'clear source medium');

    % Get shared params from caller workspace
    Nx = evalin('caller', 'Nx'); Ny = evalin('caller', 'Ny');
    dx = evalin('caller', 'dx'); dy = evalin('caller', 'dy');
    dt = evalin('caller', 'dt'); Nt = evalin('caller', 'Nt');
    t_array = evalin('caller', 't_array');
    c0 = evalin('caller', 'c0'); rho0 = evalin('caller', 'rho0');
    c1 = evalin('caller', 'c1'); rho1 = evalin('caller', 'rho1');
    BonA0 = evalin('caller', 'BonA0'); BonA1 = evalin('caller', 'BonA1');
    alpha0 = evalin('caller', 'alpha0'); alpha1 = evalin('caller', 'alpha1');
    y = evalin('caller', 'y');
    interface_position = evalin('caller', 'interface_position');
    source_strength = evalin('caller', 'source_strength');
    source_position_x = evalin('caller', 'source_position_x');
    source_position_y = evalin('caller', 'source_position_y');
    source_signal = evalin('caller', 'source_signal');
    input_args_base = evalin('caller', 'input_args_base');
    PML_alpha_def = evalin('caller', 'PML_alpha_def');

    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.t_array = t_array;

    % Medium
    medium.sound_speed = c0;
    medium.density = rho0;
    PML_alpha = PML_alpha_def;

    if any(test_num == [15:28, 57:70])
        medium.alpha_coeff = alpha0;
        medium.alpha_power = y;
    elseif any(test_num == [29:42, 71:84])
        medium.alpha_coeff = alpha0;
        medium.alpha_power = 2;
        medium.alpha_mode = 'stokes';
    end

    if test_num >= 43 && test_num <= 84
        medium.BonA = BonA0;
    end

    % Heterogeneous (even test numbers)
    if ~rem(test_num, 2)
        medium.sound_speed = c0 * ones(Nx, Ny);
        medium.density = rho0 * ones(Nx, Ny);
        medium.sound_speed(interface_position:end, :) = c1;
        medium.density(interface_position:end, :) = rho1;
        if isfield(medium, 'alpha_coeff')
            medium.alpha_coeff = alpha0 * ones(Nx, Ny);
            medium.alpha_coeff(interface_position:end, :) = alpha1;
        end
        if isfield(medium, 'BonA')
            medium.BonA = BonA0 * ones(Nx, Ny);
            medium.BonA(interface_position:end, :) = BonA1;
        end
    end

    % Source
    p0_tests = [1, 2, 15, 16, 29, 30, 43, 44, 57, 58, 71, 72];
    p_tests = [3:8, 17:22, 31:36, 45:50, 59:64, 73:78];
    u_tests = [9:14, 23:28, 37:42, 51:56, 65:70, 79:84];
    additive_tests = [3,4,9,10,17,18,23,24,31,32,37,38,45,46,51,52,59,60,65,66,73,74,79,80];
    additive_no_correction_tests = additive_tests + 2;
    dirichlet_tests = additive_tests + 4;

    if any(test_num == p0_tests)
        source.p0 = zeros(Nx, Ny);
        source.p0(source_position_x, source_position_y) = source_strength;
    elseif any(test_num == p_tests)
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(source_position_x, source_position_y) = 1;
        source.p = source_signal;
        if any(test_num == dirichlet_tests)
            source.p_mode = 'dirichlet';
        elseif any(test_num == additive_tests)
            source.p_mode = 'additive';
            PML_alpha = 0;
        else
            source.p_mode = 'additive-no-correction';
        end
    elseif any(test_num == u_tests)
        source.u_mask = zeros(Nx, Ny);
        source.u_mask(source_position_x, source_position_y) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(test_num == dirichlet_tests)
            source.u_mode = 'dirichlet';
        elseif any(test_num == additive_tests)
            source.u_mode = 'additive';
        else
            source.u_mode = 'additive-no-correction';
        end
    end

    % Full-grid sensor
    sensor.mask = ones(Nx, Ny);
    sensor.record = {'p'};

    input_args = [input_args_base, {'PMLAlpha', PML_alpha}];

    % MATLAB
    sensor_data_m = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        input_args{:}, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;

    % Python
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLSize', 10, 'PMLAlpha', PML_alpha, 'Smooth', false);
    p_python = sensor_data_p.p;
end

function plot_comparison(p_matlab, p_python, title_str, save_path)
    Nt = size(p_matlab, 2);
    t = 1:Nt;

    energy_m = sum(p_matlab.^2, 1);
    energy_p = sum(p_python.^2, 1);

    maxp_m = max(abs(p_matlab), [], 1);
    maxp_p = max(abs(p_python), [], 1);

    abs_diff = abs(p_matlab - p_python);
    max_diff = max(abs_diff, [], 1);
    rel_err = max_diff ./ max(maxp_m, eps);

    meanp_m = mean(p_matlab, 1);
    meanp_p = mean(p_python, 1);

    figure('Position', [50 50 1200 900], 'Visible', 'off');

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

    sgtitle(title_str);
    saveas(gcf, save_path);
    fprintf('  Saved: %s\n', save_path);
    close(gcf);
end

function plot_error_comparison(scenarios, all_pm, all_pp, plot_dir)
    % Overlay the relative error growth for all scenarios on one plot
    figure('Position', [50 50 800 500], 'Visible', 'off');
    colors = {'b', 'r', [0.9 0.6 0], [0.5 0 0.5]};  % blue, red, orange, purple
    styles = {'-', '--', '-.', ':'};

    for s = 1:size(scenarios, 1)
        test_num = scenarios{s, 1};
        label = scenarios{s, 2};

        pm = all_pm{s};
        pp = all_pp{s};
        Nt = size(pm, 2);
        t = 1:Nt;

        maxp_m = max(abs(pm), [], 1);
        max_diff = max(abs(pm - pp), [], 1);
        rel_err = max_diff ./ max(maxp_m, eps);

        semilogy(t, rel_err, 'Color', colors{s}, 'LineStyle', styles{s}, ...
            'LineWidth', 1.5, 'DisplayName', sprintf('#%d: %s', test_num, label));
        hold on;
    end

    xlabel('Timestep');
    ylabel('Relative Error');
    title('Error Growth: MATLAB vs Python (Selected 2D Scenarios)');
    legend('Location', 'southeast', 'FontSize', 7);
    grid on;
    ylim([1e-17, 1e-11]);

    saveas(gcf, fullfile(plot_dir, 'divergence_2d_comparison.png'));
    fprintf('  Saved: divergence_2d_comparison.png\n');
    close(gcf);
end
