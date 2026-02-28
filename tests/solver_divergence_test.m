% SOLVER_DIVERGENCE_TEST  Pinpoint the source of MATLAB/Python parity gap.
%
% The FFT accumulation test (fft_parity_test.m) shows FFTW vs pocketfft
% differences are ~4e-15 after 50 round-trips, but the solver parity gap
% is ~1e-13. This test isolates which solver operations amplify the error.
%
% Strategy: Run identical simulations, extract per-timestep sensor data,
% and measure the growth rate of divergence. Then compare against a
% PML-free simulation to test whether PML feedback amplifies the gap.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'k-Wave'));

fprintf('\n=== Solver Divergence Test ===\n');

%% Test 1: With PML (standard configuration matching parity tests)
fprintf('\n--- Test 1: Standard 2D simulation (PML inside, 32x32, 50 steps) ---\n');
[p_m1, p_p1] = run_comparison(64, 64, 50, true, 'additive-no-correction');

%% Test 2: Without PML — isolate whether PML amplifies the gap
fprintf('\n--- Test 2: No PML (PMLInside=false, expanded grid) ---\n');
[p_m2, p_p2] = run_comparison(64, 64, 50, false, 'additive-no-correction');

%% Test 3: Source-free (p0 only) — no source injection path
fprintf('\n--- Test 3: Initial pressure only (no time-varying source) ---\n');
[p_m3, p_p3] = run_comparison_p0(64, 64, 50, true);

%% Test 4: Longer 1D for reference
fprintf('\n--- Test 4: 1D simulation (64 pts, 600 steps) ---\n');
[p_m4, p_p4] = run_comparison_1d(64, 600);

fprintf('\n=== Done ===\n');

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function [p_matlab, p_python] = run_comparison(Nx, Ny, Nt, pml_inside, p_mode)
    dx = 1e-4; dy = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    source.p_mask = zeros(Nx, Ny);
    source.p_mask(Nx/4, Ny/2) = 1;
    source.p = sin(2 * pi * 1e6 * (0:Nt-1) * dt);
    source.p_mode = p_mode;

    sensor.mask = zeros(Nx, Ny);
    sensor.mask(3*Nx/4, Ny/2) = 1;
    sensor.record = {'p'};

    % MATLAB
    sensor_data_m = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PMLInside', pml_inside, 'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;

    % Python
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLInside', pml_inside);
    p_python = sensor_data_p.p;

    print_divergence(p_matlab, p_python, Nt);
end

function [p_matlab, p_python] = run_comparison_p0(Nx, Ny, Nt, pml_inside)
    dx = 1e-4; dy = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    % Initial pressure only — no time-varying source
    source.p0 = zeros(Nx, Ny);
    source.p0(Nx/2, Ny/2) = 1;

    sensor.mask = zeros(Nx, Ny);
    sensor.mask(3*Nx/4, Ny/2) = 1;
    sensor.record = {'p'};

    % MATLAB
    sensor_data_m = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PMLInside', pml_inside, 'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;

    % Python
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLInside', pml_inside);
    p_python = sensor_data_p.p;

    print_divergence(p_matlab, p_python, Nt);
end

function [p_matlab, p_python] = run_comparison_1d(Nx, Nt)
    dx = 1e-4; dt = 1e-8;
    kgrid = kWaveGrid(Nx, dx);
    kgrid.setTime(Nt, dt);

    medium.sound_speed = 1500;
    medium.density = 1000;

    source.p0 = zeros(Nx, 1);
    source.p0(Nx/4) = 1;

    sensor.mask = zeros(Nx, 1);
    sensor.mask(3*Nx/4) = 1;
    sensor.record = {'p'};

    % MATLAB
    sensor_data_m = kspaceFirstOrder1D(kgrid, medium, source, sensor, ...
        'PlotSim', false, 'DataCast', 'double');
    p_matlab = sensor_data_m.p;

    % Python
    sensor_data_p = kspaceFirstOrderPy(kgrid, medium, source, sensor);
    p_python = sensor_data_p.p;

    print_divergence(p_matlab, p_python, Nt);
end

function print_divergence(p_matlab, p_python, Nt)
    abs_diff = abs(p_matlab - p_python);
    p_max = max(abs(p_matlab(:)));
    if p_max < eps
        fprintf('  Signal too small to compute relative error.\n');
        return;
    end
    rel_diff = abs_diff / p_max;

    % Print at key timesteps
    steps = unique([1, 2, 5, 10, round(Nt/4), round(Nt/2), round(3*Nt/4), Nt]);
    steps = steps(steps <= length(p_matlab));
    fprintf('  t     | Rel. error   | Abs. error\n');
    fprintf('  ------|--------------|-------------\n');
    for t = steps
        fprintf('  %-5d | %.3e  | %.3e\n', t, rel_diff(t), abs_diff(t));
    end
    fprintf('  Peak rel. error: %.3e (at t=%d)\n', max(rel_diff), find(rel_diff == max(rel_diff), 1));
end
