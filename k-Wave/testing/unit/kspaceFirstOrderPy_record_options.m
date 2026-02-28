function test_pass = kspaceFirstOrderPy_record_options(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Tests sensor.record variants using the Python backend. Validates
%     aggregate records (p_max, p_min, p_rms), velocity aggregates,
%     final-state snapshots, record expansion (u -> ux, uy), and
%     intensity records. All checks are Python-only self-consistency.
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

source.p0 = zeros(Nx, Ny);
source.p0(Nx/2, Ny/2) = 1e6;

Nt = 20;
kgrid.setTime(Nt, 0.1 * dx / medium.sound_speed);

comparison_thresh = 1e-10;

sim_args = {'PMLSize', 0, 'Smooth', false};

% =========================================================================
% SUB-TEST 1: AGGREGATE RECORDS (p_max, p_min, p_rms)
% =========================================================================

disp('--- Sub-test 1: Aggregate records (p, p_max, p_min, p_rms) ---');

sensor1.mask = zeros(Nx, Ny);
sensor1.mask(Nx/4, Ny/2) = 1;
sensor1.mask(3*Nx/4, Ny/2) = 1;
sensor1.record = {'p', 'p_max', 'p_min', 'p_rms'};

out1 = kspaceFirstOrderPy(kgrid, medium, source, sensor1, sim_args{:});

% verify p_max == max(p, [], 2)
% Note: (:) forces column vectors to avoid row/column broadcast mismatches
% when Python returns 1D arrays that MATLAB may interpret as row vectors.
expected_p_max = max(out1.p, [], 2);
err_p_max = max(abs(out1.p_max(:) - expected_p_max(:)));
if err_p_max > comparison_thresh
    test_pass = false;
    fprintf('FAILED: p_max mismatch (err=%.2e)\n', err_p_max);
else
    fprintf('PASSED: p_max matches max(p) (err=%.2e)\n', err_p_max);
end

% verify p_min == min(p, [], 2)
expected_p_min = min(out1.p, [], 2);
err_p_min = max(abs(out1.p_min(:) - expected_p_min(:)));
if err_p_min > comparison_thresh
    test_pass = false;
    fprintf('FAILED: p_min mismatch (err=%.2e)\n', err_p_min);
else
    fprintf('PASSED: p_min matches min(p) (err=%.2e)\n', err_p_min);
end

% verify p_rms == sqrt(mean(p.^2, 2))
expected_p_rms = sqrt(mean(out1.p.^2, 2));
err_p_rms = max(abs(out1.p_rms(:) - expected_p_rms(:)));
if err_p_rms > comparison_thresh
    test_pass = false;
    fprintf('FAILED: p_rms mismatch (err=%.2e)\n', err_p_rms);
else
    fprintf('PASSED: p_rms matches sqrt(mean(p^2)) (err=%.2e)\n', err_p_rms);
end

% =========================================================================
% SUB-TEST 2: VELOCITY AGGREGATES (u_max, u_rms)
% =========================================================================

disp('--- Sub-test 2: Velocity aggregates (u_non_staggered, u_max, u_rms) ---');

% NOTE: MATLAB 'u' records staggered velocity, but Python computes u_max
% from colocated (non-staggered) data. Use 'u_non_staggered' so the
% recorded time series and aggregates refer to the same colocated field.
sensor2.mask = zeros(Nx, Ny);
sensor2.mask(Nx/4, Ny/2) = 1;
sensor2.mask(3*Nx/4, Ny/2) = 1;
sensor2.record = {'u_non_staggered', 'u_max', 'u_rms'};

out2 = kspaceFirstOrderPy(kgrid, medium, source, sensor2, sim_args{:});

% verify ux_non_staggered, uy_non_staggered fields exist
if ~isfield(out2, 'ux_non_staggered') || ~isfield(out2, 'uy_non_staggered')
    test_pass = false;
    disp('FAILED: Missing ux_non_staggered or uy_non_staggered fields');
else
    disp('PASSED: ux_non_staggered and uy_non_staggered fields present');
end

% verify ux_max, uy_max fields exist and match colocated time series
if ~isfield(out2, 'ux_max') || ~isfield(out2, 'uy_max')
    test_pass = false;
    disp('FAILED: Missing ux_max or uy_max fields');
else
    expected_ux_max = max(out2.ux_non_staggered, [], 2);
    err_ux_max = max(abs(out2.ux_max(:) - expected_ux_max(:)));
    if err_ux_max > comparison_thresh
        test_pass = false;
        fprintf('FAILED: ux_max mismatch (err=%.2e)\n', err_ux_max);
    else
        fprintf('PASSED: ux_max matches max(ux_non_staggered) (err=%.2e)\n', err_ux_max);
    end

    expected_uy_max = max(out2.uy_non_staggered, [], 2);
    err_uy_max = max(abs(out2.uy_max(:) - expected_uy_max(:)));
    if err_uy_max > comparison_thresh
        test_pass = false;
        fprintf('FAILED: uy_max mismatch (err=%.2e)\n', err_uy_max);
    else
        fprintf('PASSED: uy_max matches max(uy_non_staggered) (err=%.2e)\n', err_uy_max);
    end
end

% verify ux_rms, uy_rms fields exist and match colocated time series
if ~isfield(out2, 'ux_rms') || ~isfield(out2, 'uy_rms')
    test_pass = false;
    disp('FAILED: Missing ux_rms or uy_rms fields');
else
    expected_ux_rms = sqrt(mean(out2.ux_non_staggered.^2, 2));
    err_ux_rms = max(abs(out2.ux_rms(:) - expected_ux_rms(:)));
    if err_ux_rms > comparison_thresh
        test_pass = false;
        fprintf('FAILED: ux_rms mismatch (err=%.2e)\n', err_ux_rms);
    else
        fprintf('PASSED: ux_rms matches sqrt(mean(ux_ns^2)) (err=%.2e)\n', err_ux_rms);
    end
end

% =========================================================================
% SUB-TEST 3: FINAL-STATE SNAPSHOTS (p_final, u_final)
% =========================================================================

disp('--- Sub-test 3: Final-state snapshots (p_final, u_final) ---');

sensor3.mask = zeros(Nx, Ny);
sensor3.mask(Nx/4, Ny/2) = 1;
sensor3.mask(3*Nx/4, Ny/2) = 1;
sensor3.record = {'p_final', 'u_final'};

out3 = kspaceFirstOrderPy(kgrid, medium, source, sensor3, sim_args{:});

% verify p_final exists and has correct size (full grid snapshot: Nx*Ny)
if ~isfield(out3, 'p_final')
    test_pass = false;
    disp('FAILED: p_final field missing');
else
    p_final_numel = numel(out3.p_final);
    if p_final_numel == Nx * Ny
        fprintf('PASSED: p_final has correct size (%d elements)\n', p_final_numel);
    else
        test_pass = false;
        fprintf('FAILED: p_final expected %d elements, got %d\n', Nx * Ny, p_final_numel);
    end
end

% verify u_final expands to ux_final, uy_final
if ~isfield(out3, 'ux_final') || ~isfield(out3, 'uy_final')
    test_pass = false;
    disp('FAILED: u_final did not expand to ux_final and uy_final');
else
    ux_final_numel = numel(out3.ux_final);
    if ux_final_numel == Nx * Ny
        fprintf('PASSED: ux_final has correct size (%d elements)\n', ux_final_numel);
    else
        test_pass = false;
        fprintf('FAILED: ux_final expected %d elements, got %d\n', Nx * Ny, ux_final_numel);
    end
end

% =========================================================================
% SUB-TEST 4: RECORD EXPANSION (u -> ux, uy)
% =========================================================================

disp('--- Sub-test 4: Record expansion (u -> ux, uy) ---');

sensor4.mask = zeros(Nx, Ny);
sensor4.mask(Nx/4, Ny/2) = 1;
sensor4.mask(3*Nx/4, Ny/2) = 1;
sensor4.record = {'u'};

out4 = kspaceFirstOrderPy(kgrid, medium, source, sensor4, sim_args{:});

if ~isfield(out4, 'ux')
    test_pass = false;
    disp('FAILED: ux field missing when recording u');
else
    disp('PASSED: ux field present');
end

if ~isfield(out4, 'uy')
    test_pass = false;
    disp('FAILED: uy field missing when recording u');
else
    disp('PASSED: uy field present');
end

% verify ux has expected shape: n_sensor_points x Nt
expected_rows = 2;  % two sensor points
if isfield(out4, 'ux')
    if size(out4.ux, 1) ~= expected_rows || size(out4.ux, 2) ~= Nt
        test_pass = false;
        fprintf('FAILED: ux size expected [%d %d], got [%d %d]\n', ...
            expected_rows, Nt, size(out4.ux, 1), size(out4.ux, 2));
    else
        fprintf('PASSED: ux has correct size [%d %d]\n', size(out4.ux, 1), size(out4.ux, 2));
    end
end

% =========================================================================
% SUB-TEST 5: INTENSITY RECORDS (I, I_avg)
% =========================================================================

disp('--- Sub-test 5: Intensity records (p, u_non_staggered, I, I_avg) ---');

sensor5.mask = zeros(Nx, Ny);
sensor5.mask(Nx/4, Ny/2) = 1;
sensor5.mask(3*Nx/4, Ny/2) = 1;
sensor5.record = {'p', 'u_non_staggered', 'I', 'I_avg'};

out5 = kspaceFirstOrderPy(kgrid, medium, source, sensor5, sim_args{:});

% verify Ix and Iy fields exist
if ~isfield(out5, 'Ix') || ~isfield(out5, 'Iy')
    test_pass = false;
    disp('FAILED: Ix or Iy fields missing');
else
    disp('PASSED: Ix and Iy fields present');
end

% verify Ix_avg and Iy_avg fields exist
if ~isfield(out5, 'Ix_avg') || ~isfield(out5, 'Iy_avg')
    test_pass = false;
    disp('FAILED: Ix_avg or Iy_avg fields missing');
else
    disp('PASSED: Ix_avg and Iy_avg fields present');
end

% verify u_non_staggered expanded to ux_non_staggered, uy_non_staggered
if ~isfield(out5, 'ux_non_staggered') || ~isfield(out5, 'uy_non_staggered')
    test_pass = false;
    disp('FAILED: ux_non_staggered or uy_non_staggered fields missing');
else
    disp('PASSED: ux_non_staggered and uy_non_staggered fields present');
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    % plot pressure time series from sub-test 1
    subplot(2, 2, 1);
    plot(out1.p(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(out1.p(2, :), 'r-', 'LineWidth', 1.5);
    yline(out1.p_max(1), 'b--');
    yline(out1.p_max(2), 'r--');
    yline(out1.p_min(1), 'b:');
    yline(out1.p_min(2), 'r:');
    legend('Sensor 1 p', 'Sensor 2 p', 'S1 p\_max', 'S2 p\_max', 'S1 p\_min', 'S2 p\_min');
    title('Pressure with Aggregates');
    xlabel('Time Step');
    ylabel('Pressure');

    % plot velocity from sub-test 2
    subplot(2, 2, 2);
    plot(out2.ux_non_staggered(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(out2.uy_non_staggered(1, :), 'r-', 'LineWidth', 1.5);
    legend('ux sensor 1', 'uy sensor 1');
    title('Velocity Components (colocated)');
    xlabel('Time Step');
    ylabel('Velocity');

    % plot final state from sub-test 3
    subplot(2, 2, 3);
    imagesc(reshape(out3.p_final, [Nx, Ny]));
    title('p\_final Snapshot');
    colorbar; axis image;

    % plot intensity from sub-test 5
    subplot(2, 2, 4);
    plot(out5.Ix(1, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(out5.Iy(1, :), 'r-', 'LineWidth', 1.5);
    yline(out5.Ix_avg(1), 'b--');
    yline(out5.Iy_avg(1), 'r--');
    legend('Ix', 'Iy', 'Ix\_avg', 'Iy\_avg');
    title('Intensity');
    xlabel('Time Step');
    ylabel('Intensity');
end
