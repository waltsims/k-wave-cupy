function test_pass = kspaceFirstOrder1D_velocity_recording(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Compare velocity recording between MATLAB and Python backends.
%     Tests both staggered ('u') and colocated ('u_non_staggered') velocity.
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
sensor.record = {'p', 'u', 'u_non_staggered'};

% =========================================================================
% MATLAB REFERENCE
% =========================================================================

opts = {'PlotSim', false, 'PMLAlpha', 0, 'PMLSize', 0, ...
        'PMLInside', true, 'DataCast', 'off', 'DisplayMask', 'off', ...
        'PlotPML', false, 'Smooth', false};
mat = kspaceFirstOrder1D(kgrid, medium, source, sensor, opts{:});

% =========================================================================
% PYTHON BACKEND
% =========================================================================

py_result = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
    'PMLSize', 0, 'PMLAlpha', 0, 'Smooth', false);

% =========================================================================
% COMPARE
% =========================================================================

% Pressure
diff_p = max(abs(mat.p(:) - py_result.p(:)));
fprintf('  Pressure max diff:              %e\n', diff_p);
if diff_p > comparison_thresh
    test_pass = false;
    fprintf('  FAILED: pressure exceeds threshold\n');
end

% Staggered velocity (ux)
diff_ux = max(abs(mat.ux(:) - py_result.ux(:)));
fprintf('  Staggered velocity max diff:    %e\n', diff_ux);
if diff_ux > comparison_thresh
    test_pass = false;
    fprintf('  FAILED: staggered ux exceeds threshold\n');
end

% Colocated velocity (ux_non_staggered)
diff_ux_ns = max(abs(mat.ux_non_staggered(:) - py_result.ux_non_staggered(:)));
fprintf('  Colocated velocity max diff:    %e\n', diff_ux_ns);
if diff_ux_ns > comparison_thresh
    test_pass = false;
    fprintf('  FAILED: colocated ux exceeds threshold\n');
end

% =========================================================================
% PLOT
% =========================================================================

if plot_comparisons
    figure;
    subplot(2,3,1); imagesc(mat.ux);      title('MATLAB ux (staggered)');
    subplot(2,3,2); imagesc(py_result.ux); title('Python ux (staggered)');
    subplot(2,3,3); imagesc(mat.ux - py_result.ux); title('Diff staggered'); colorbar;
    subplot(2,3,4); imagesc(mat.ux_non_staggered);      title('MATLAB ux\_ns');
    subplot(2,3,5); imagesc(py_result.ux_non_staggered); title('Python ux\_ns');
    subplot(2,3,6); imagesc(mat.ux_non_staggered - py_result.ux_non_staggered); title('Diff colocated'); colorbar;
end

end
