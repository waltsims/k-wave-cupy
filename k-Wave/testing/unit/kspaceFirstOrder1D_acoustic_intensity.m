function test_pass = kspaceFirstOrder1D_acoustic_intensity(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Compare acoustic intensity computed by Python's acoustic_intensity()
%     against MATLAB's fourierShift-based computation. Both backends run
%     with {'p', 'u_non_staggered'}, then intensity is computed post-hoc.
%
% ABOUT:
%     author      - Walter
%     date        - 21st February 2026
%     last update - 21st February 2026

if nargin == 0
    plot_comparisons = false;
    plot_simulations = false;
end

test_pass = true;
comparison_thresh = 1e-6;

% =========================================================================
% SETUP
% =========================================================================

Nx = 64;
dx = 1e-4;

kgrid = kWaveGrid(Nx, dx);
kgrid.makeTime(1500, [], 2e-6);

medium.sound_speed = 1500;
medium.density     = 1000;

source.p0 = zeros(Nx, 1);
source.p0(Nx/2) = 1;

sensor.mask   = ones(Nx, 1);
sensor.record = {'p', 'u_non_staggered'};

% =========================================================================
% MATLAB REFERENCE (compute intensity manually via fourierShift)
% =========================================================================

opts = {'PlotSim', false, 'PMLAlpha', 0, 'PMLSize', 0, ...
        'PMLInside', true, 'DataCast', 'off', 'DisplayMask', 'off', ...
        'PlotPML', false, 'Smooth', false};
mat = kspaceFirstOrder1D(kgrid, medium, source, sensor, opts{:});

% Compute intensity: shift velocity forward by dt/2, then I = p * u
ux_shifted = fourierShift(mat.ux_non_staggered, 1/2);
Ix_matlab = mat.p .* ux_shifted;
Ix_avg_matlab = mean(Ix_matlab, 2);

% =========================================================================
% PYTHON BACKEND (compute intensity via acoustic_intensity helper)
% =========================================================================

py_result = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
    'PMLSize', 0, 'PMLAlpha', 0, 'Smooth', false);

% Call Python acoustic_intensity helper
toNumpy = @(x) py.numpy.array(double(x), pyargs('order', 'F'));
py_dict = py.dict(pyargs('p', toNumpy(py_result.p), ...
                          'ux', toNumpy(py_result.ux_non_staggered)));
py_I = py.kWavePy.acoustic_intensity(py_dict);
Ix_py = double(py_I{'Ix'});
Ix_avg_py = double(py_I{'Ix_avg'});

% =========================================================================
% COMPARE
% =========================================================================

% Time-varying intensity
diff_Ix = max(abs(Ix_matlab(:) - Ix_py(:)));
fprintf('  Ix max diff:      %e\n', diff_Ix);
if diff_Ix > comparison_thresh
    test_pass = false;
    fprintf('  FAILED: Ix exceeds threshold\n');
end

% Time-averaged intensity
diff_Ix_avg = max(abs(Ix_avg_matlab(:) - Ix_avg_py(:)));
fprintf('  Ix_avg max diff:  %e\n', diff_Ix_avg);
if diff_Ix_avg > comparison_thresh
    test_pass = false;
    fprintf('  FAILED: Ix_avg exceeds threshold\n');
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
    figure;
    subplot(2,3,1); imagesc(Ix_matlab); title('MATLAB Ix'); colorbar;
    subplot(2,3,2); imagesc(Ix_py);     title('Python Ix'); colorbar;
    subplot(2,3,3); imagesc(Ix_matlab - Ix_py); title('Diff Ix'); colorbar;
    subplot(2,3,4); plot(Ix_avg_matlab); title('MATLAB Ix\_avg');
    subplot(2,3,5); plot(Ix_avg_py);     title('Python Ix\_avg');
    subplot(2,3,6); plot(Ix_avg_matlab - Ix_avg_py); title('Diff Ix\_avg');
end

end
