function test_pass = kspaceFirstOrderPy_pml_inside_false_no_final(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Test that PMLInside=false works when sensor.record contains no
%     _final fields (e.g., only 'p'). Regression test for a crash where
%     final_fields{1} was indexed on an empty cell array.
%
% ABOUT:
%     author      - k-Wave-CuPy Development
%     date        - 21st March 2026
%     last update - 21st March 2026
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

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

% =========================================================================
% TEST: PMLInside=false with sensor.record = {'p'} (no _final fields)
% =========================================================================

Nx = 64;
Ny = 64;
dx = 0.1e-3;
dy = 0.1e-3;
pml = 10;

kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.makeTime(1500);
kgrid.setTime(min(20, kgrid.Nt), kgrid.dt);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = zeros(Nx, Ny);
source.p0(Nx/2, Ny/2) = 5;

% Record only time-series, no _final fields
sensor.record = {'p'};

disp('Running PMLInside=false with record={''p''} (no _final fields)...');
try
    result = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLSize', pml, 'PMLInside', false, 'Smooth', false);
catch ME
    fprintf('FAILED: Crashed with error: %s\n', ME.message);
    test_pass = false;
    return;
end

% Verify we got time-series data with correct dimensions
fprintf('sensor_data.p size: [%s]\n', num2str(size(result.p)));

if isempty(result.p)
    fprintf('FAILED: result.p is empty\n');
    test_pass = false;
else
    disp('PASSED: PMLInside=false with no _final fields runs without error');
end
