function test_pass = kspaceFirstOrder1D_photoacoustic_parity(plot_comparisons, plot_simulations)
%KSPACEFIRSTORDER1D_PHOTOACOUSTIC_PARITY Compare Python and MATLAB backends
%   with default settings (including smoothing).
%
% DESCRIPTION:
%     Runs the 1D portion of example_ivp_photoacoustic_waveforms through
%     both the Python and MATLAB backends and checks that the results
%     match. This catches mismatches caused by differing default
%     preprocessing (e.g., smoothing of p0).
%
% ABOUT:
%     author      - Walter
%     date        - 22nd February 2026
%     last update - 22nd February 2026
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

if nargin == 0
    plot_comparisons = false;
    plot_simulations = false;
end

test_pass = true;
comparison_thresh = 1e-6;

% Check Python availability
try
    py.numpy.array(1);
catch
    warning('Python/NumPy not available. Skipping test.');
    return;
end

% =========================================================================
% SETUP (from example_ivp_photoacoustic_waveforms, 1D section)
% =========================================================================

Nx = 64;
x = 1e-3;
dx = x / Nx;

medium.sound_speed = 1500;

source_radius = 2;
source_sensor_distance = 10;

source.p0 = zeros(Nx, 1);
source.p0(Nx/2 - source_radius:Nx/2 + source_radius) = 1;

sensor.mask = zeros(Nx, 1);
sensor.mask(Nx/2 + source_sensor_distance) = 1;

dt = 2e-9;
t_end = 300e-9;
kgrid = kWaveGrid(Nx, dx);
kgrid.setTime(round(t_end / dt) + 1, dt);

% =========================================================================
% RUN BOTH BACKENDS WITH DEFAULT SETTINGS
% =========================================================================

% Python backend (via wrapper, should apply smoothing by default)
py_data = kspaceFirstOrderPy(kgrid, medium, source, sensor);

% MATLAB backend with default settings (smooths p0 by default)
% PlotSim must be off for headless testing; all other defaults preserved
mat_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, 'PlotSim', false);

% =========================================================================
% COMPARE
% =========================================================================

diff = max(abs(double(py_data(:)) - double(mat_data(:))));
fprintf('  Max difference (default settings): %e\n', diff);

if diff > comparison_thresh
    fprintf('FAIL: Python/MATLAB mismatch exceeds threshold %e\n', comparison_thresh);
    test_pass = false;
end

if plot_comparisons
    figure;
    subplot(3,1,1); plot(py_data.'); title('Python (default)');
    subplot(3,1,2); plot(mat_data.'); title('MATLAB (default)');
    subplot(3,1,3); plot(py_data.' - mat_data.'); title('Difference');
    sgtitle('Photoacoustic Parity (1D)');
end

end
