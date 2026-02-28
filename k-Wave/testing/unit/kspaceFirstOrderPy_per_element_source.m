function test_pass = kspaceFirstOrderPy_per_element_source(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Tests that 2D source signal arrays (n_src x Nt) produce different
%     waveforms at different source points. Verifies per-element source
%     injection by checking that three point sources with amplitudes in
%     a 1:2:3 ratio produce distinguishable recorded signals.
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
% SIMULATION SETUP
% =========================================================================

Nx = 32;
Ny = 32;
dx = 1;
dy = 1;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1500;
medium.density = 1000;

Nt = 30;
kgrid.setTime(Nt, 0.1 * dx / medium.sound_speed);

% create 3 point sources at different locations
source.p_mask = zeros(Nx, Ny);
src_positions = [8, 16; 16, 16; 24, 16];  % 3 source points
for k = 1:3
    source.p_mask(src_positions(k, 1), src_positions(k, 2)) = 1;
end

% each source gets a different amplitude sinusoid
t_array = (0:Nt-1) * kgrid.dt;
freq = 200;
source.p = [1e6 * sin(2*pi*freq*t_array); ...
            2e6 * sin(2*pi*freq*t_array); ...
            3e6 * sin(2*pi*freq*t_array)];
source.p_mode = 'additive';

% sensor at the same 3 locations
sensor.mask = source.p_mask;

% =========================================================================
% RUN SIMULATION
% =========================================================================

disp('Running Python backend with per-element source...');
sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
    'PMLSize', 0, 'Smooth', false);

% =========================================================================
% VERIFICATION
% =========================================================================

% --- Check 1: Output has 3 rows (one per sensor point) ---
disp('--- Check 1: Output size ---');
n_sensors = size(sensor_data, 1);
if n_sensors ~= 3
    test_pass = false;
    fprintf('FAILED: Expected 3 sensor rows, got %d\n', n_sensors);
else
    disp('PASSED: sensor_data has 3 rows (one per sensor point)');
end

% --- Check 2: Recorded signals are NOT all identical ---
disp('--- Check 2: Signals are distinguishable ---');
maxvals = max(abs(sensor_data), [], 2);

fprintf('  Max amplitude sensor 1: %.4e\n', maxvals(1));
fprintf('  Max amplitude sensor 2: %.4e\n', maxvals(2));
fprintf('  Max amplitude sensor 3: %.4e\n', maxvals(3));

% check that signals differ (ratio of max to min amplitude > 1.5)
amplitude_ratio = max(maxvals) / min(maxvals);
fprintf('  Amplitude ratio (max/min): %.2f\n', amplitude_ratio);

if amplitude_ratio < 1.5
    test_pass = false;
    disp('FAILED: Per-element sources produced nearly identical signals');
else
    disp('PASSED: Per-element sources produced distinguishable signals');
end

% --- Check 3: Pairwise signals are different ---
disp('--- Check 3: Pairwise signal differences ---');
for i = 1:2
    for j = (i+1):3
        diff_ij = max(abs(sensor_data(i, :) - sensor_data(j, :)));
        if diff_ij < 1e-10
            test_pass = false;
            fprintf('FAILED: Sensors %d and %d are identical (max diff=%.2e)\n', i, j, diff_ij);
        else
            fprintf('PASSED: Sensors %d and %d differ (max diff=%.2e)\n', i, j, diff_ij);
        end
    end
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    subplot(2, 1, 1);
    plot(sensor_data(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(sensor_data(2, :), 'r-', 'LineWidth', 1.5);
    plot(sensor_data(3, :), 'g-', 'LineWidth', 1.5);
    legend('Source 1 (1x)', 'Source 2 (2x)', 'Source 3 (3x)');
    title('Recorded Pressure at Each Sensor');
    xlabel('Time Step');
    ylabel('Pressure');

    subplot(2, 1, 2);
    bar(maxvals);
    xlabel('Sensor Index');
    ylabel('Max Absolute Pressure');
    title('Peak Amplitude per Sensor');
    set(gca, 'XTickLabel', {'1x amp', '2x amp', '3x amp'});
end
