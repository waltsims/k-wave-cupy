function test_pass = kspaceFirstOrderPy_pml_inside_false_2D(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Test that PMLInside=false produces identical results to manually
%     expanding the grid with PMLInside=true (default). Both cases present
%     the same expanded grid to the Python solver, so p_final should match
%     to machine precision.
%
% ABOUT:
%     author      - k-Wave-CuPy Development
%     date        - 22nd February 2026
%     last update - 22nd February 2026
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

comparison_thresh = 1e-14;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

pml = 10;
Nx = 64;
Ny = 64;
dx = 0.1e-3;
dy = 0.1e-3;

% homogeneous medium
medium.sound_speed = 1500;
medium.density = 1000;

% =========================================================================
% BASELINE: Manually expanded grid with PMLInside = true (default)
% =========================================================================

kgrid_big = kWaveGrid(Nx + 2*pml, dx, Ny + 2*pml, dy);
kgrid_big.makeTime(medium.sound_speed);
kgrid_big.setTime(min(50, kgrid_big.Nt), kgrid_big.dt);

% p0: single point source at physical center
source_big.p0 = zeros(Nx + 2*pml, Ny + 2*pml);
source_big.p0(Nx/2 + pml, Ny/2 + pml) = 5;

sensor_big.record = {'p_final'};

disp('Running baseline (expanded grid, PMLInside=true)...');
result_big = kspaceFirstOrderPy(kgrid_big, medium, source_big, sensor_big, ...
    'PMLSize', pml, 'Smooth', false);

% =========================================================================
% TEST: Original grid with PMLInside = false
% =========================================================================

kgrid_small = kWaveGrid(Nx, dx, Ny, dy);
kgrid_small.setTime(kgrid_big.Nt, kgrid_big.dt);

% p0: same physical location (center of user domain)
source_small.p0 = zeros(Nx, Ny);
source_small.p0(Nx/2, Ny/2) = 5;

sensor_small.record = {'p_final'};

disp('Running test (original grid, PMLInside=false)...');
result_small = kspaceFirstOrderPy(kgrid_small, medium, source_small, sensor_small, ...
    'PMLSize', pml, 'PMLInside', false, 'Smooth', false);

% =========================================================================
% COMPARISON
% =========================================================================

% p_final from Python already trims PML, so both should be Nx x Ny
fprintf('Baseline p_final size: [%s]\n', num2str(size(result_big.p_final)));
fprintf('Test p_final size:     [%s]\n', num2str(size(result_small.p_final)));

diff_val = abs(result_big.p_final(:) - result_small.p_final(:));
max_diff = max(diff_val);
fprintf('Max absolute difference: %.2e\n', max_diff);

if max_diff > comparison_thresh
    test_pass = false;
    fprintf('FAILED: Max difference %.2e exceeds threshold %.2e\n', max_diff, comparison_thresh);
else
    disp('PASSED: PMLInside=false matches manually expanded grid');
end

% =========================================================================
% VISUALISATION
% =========================================================================

if plot_comparisons
    figure;

    subplot(1, 3, 1);
    imagesc(result_big.p_final);
    title('Baseline (expanded grid)');
    colorbar; axis image;

    subplot(1, 3, 2);
    imagesc(result_small.p_final);
    title('Test (PMLInside=false)');
    colorbar; axis image;

    subplot(1, 3, 3);
    imagesc(result_big.p_final - result_small.p_final);
    title('Difference');
    colorbar; axis image;
end
