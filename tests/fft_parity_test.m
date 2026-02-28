% FFT_PARITY_TEST  Compare MATLAB (FFTW) vs NumPy (pocketfft) FFT outputs.
%
% Tests the hypothesis that the ~1e-12 parity gap between MATLAB and Python
% solvers in 2D comes from differences between FFT library implementations
% (FFTW vs pocketfft), not from a model error.
%
% Methodology:
%   1. Generate deterministic test arrays matching simulation grid sizes
%   2. Pass identical data to both MATLAB fft/fft2 and numpy.fft.fft/fft2
%   3. Measure element-wise relative differences
%   4. Repeat for multiple FFT round-trips (fft -> ifft -> fft -> ...) to
%      simulate error accumulation over time steps

% check for plot inputs
if ~exist('plot_results', 'var')
    plot_results = true;
end

%% Setup
rng(42);  % deterministic

% Match the grid sizes used in the parity tests
Nx_1d = 64;       % 1D test grid
Nx_2d = 32;       % 2D test grid
Ny_2d = 32;
n_roundtrips = 50; % matches 2D test Nt=50

fprintf('\n=== FFT Parity Test: MATLAB (FFTW) vs NumPy (pocketfft) ===\n\n');

%% 1D single FFT comparison
data_1d = randn(Nx_1d, 1);

% MATLAB FFT
fft_matlab_1d = fft(data_1d);

% NumPy FFT (same input)
fft_numpy_1d = double(py.numpy.array(py.numpy.fft.fft(py.numpy.array(data_1d'))));

% Compare
diff_1d = abs(fft_matlab_1d(:) - fft_numpy_1d(:));
rel_diff_1d = diff_1d ./ max(abs(fft_matlab_1d(:)), eps);
max_rel_1d = max(rel_diff_1d);
mean_rel_1d = mean(rel_diff_1d);

fprintf('--- 1D FFT (N=%d) ---\n', Nx_1d);
fprintf('  Max  element-wise relative diff: %.2e\n', max_rel_1d);
fprintf('  Mean element-wise relative diff: %.2e\n', mean_rel_1d);
fprintf('  Expected bound (eps*log2(N)):    %.2e\n', eps * log2(Nx_1d));

%% 2D single FFT comparison
data_2d = randn(Nx_2d, Ny_2d);

% MATLAB FFT2
fft_matlab_2d = fft2(data_2d);

% NumPy FFT2 — pass as column-major (Fortran order)
py_data_2d = py.numpy.array(data_2d(:)', pyargs('dtype', 'float64'));
py_data_2d = py_data_2d.reshape(int32([Nx_2d, Ny_2d]), pyargs('order', 'F'));
fft_numpy_2d_py = py.numpy.fft.fft2(py_data_2d);

% Convert back to MATLAB
fft_numpy_2d_real = double(py.numpy.real(fft_numpy_2d_py).flatten(pyargs('order', 'F')));
fft_numpy_2d_imag = double(py.numpy.imag(fft_numpy_2d_py).flatten(pyargs('order', 'F')));
fft_numpy_2d = reshape(complex(fft_numpy_2d_real, fft_numpy_2d_imag), [Nx_2d, Ny_2d]);

% Compare
diff_2d = abs(fft_matlab_2d(:) - fft_numpy_2d(:));
rel_diff_2d = diff_2d ./ max(abs(fft_matlab_2d(:)), eps);
max_rel_2d = max(rel_diff_2d);
mean_rel_2d = mean(rel_diff_2d);

fprintf('\n--- 2D FFT (%dx%d) ---\n', Nx_2d, Ny_2d);
fprintf('  Max  element-wise relative diff: %.2e\n', max_rel_2d);
fprintf('  Mean element-wise relative diff: %.2e\n', mean_rel_2d);
fprintf('  Expected bound (eps*log2(N^2)):  %.2e\n', eps * log2(Nx_2d * Ny_2d));

%% Accumulation test: repeated FFT->IFFT round-trips
% Simulates how per-FFT errors compound through time-stepping.
% Start with identical data, apply fft->multiply->ifft N times in both
% MATLAB and NumPy, measure how the results diverge.

% Use a spectral operator similar to what the solver applies (k-space shift)
kx = (-Nx_2d/2:Nx_2d/2-1) / Nx_2d;
ky = (-Ny_2d/2:Ny_2d/2-1) / Ny_2d;
[KX, KY] = meshgrid(kx, ky);
KX = ifftshift(KX'); KY = ifftshift(KY');
spectral_op = exp(-0.01 * (KX.^2 + KY.^2));  % mild spectral filter

field_matlab = data_2d;

% Setup Python side
py_field = py.numpy.array(data_2d(:)', pyargs('dtype', 'float64'));
py_field = py_field.reshape(int32([Nx_2d, Ny_2d]), pyargs('order', 'F'));
py_spectral_op = py.numpy.array(spectral_op(:)', pyargs('dtype', 'float64'));
py_spectral_op = py_spectral_op.reshape(int32([Nx_2d, Ny_2d]), pyargs('order', 'F'));

divergence = zeros(n_roundtrips, 1);

for k = 1:n_roundtrips
    % MATLAB: fft2 -> multiply -> ifft2
    field_matlab = real(ifft2(spectral_op .* fft2(field_matlab)));

    % NumPy: same operation
    py_fft = py.numpy.fft.fft2(py_field);
    py_filtered = py.numpy.multiply(py_spectral_op, py_fft);
    py_field = py.numpy.real(py.numpy.fft.ifft2(py_filtered));

    % Extract NumPy result for comparison
    py_vals = double(py_field.flatten(pyargs('order', 'F')));
    field_numpy = reshape(py_vals, [Nx_2d, Ny_2d]);

    % Measure divergence
    rel_err = max(abs(field_matlab(:) - field_numpy(:))) / max(abs(field_matlab(:)));
    divergence(k) = rel_err;
end

fprintf('\n--- Accumulation over %d FFT round-trips (%dx%d) ---\n', n_roundtrips, Nx_2d, Ny_2d);
fprintf('  After  1 round-trip:  %.2e\n', divergence(1));
fprintf('  After 10 round-trips: %.2e\n', divergence(10));
fprintf('  After 25 round-trips: %.2e\n', divergence(25));
fprintf('  After %d round-trips: %.2e\n', n_roundtrips, divergence(end));

%% 1D accumulation for comparison
field_matlab_1d = data_1d;
py_field_1d = py.numpy.array(data_1d', pyargs('dtype', 'float64'));
spectral_op_1d = exp(-0.01 * ((-Nx_1d/2:Nx_1d/2-1) / Nx_1d).^2);
spectral_op_1d = ifftshift(spectral_op_1d(:));
py_spectral_op_1d = py.numpy.array(spectral_op_1d', pyargs('dtype', 'float64'));

n_roundtrips_1d = 600;  % matches 1D test Nt=600
divergence_1d = zeros(n_roundtrips_1d, 1);

for k = 1:n_roundtrips_1d
    field_matlab_1d = real(ifft(spectral_op_1d .* fft(field_matlab_1d)));
    py_fft_1d = py.numpy.fft.fft(py_field_1d);
    py_field_1d = py.numpy.ascontiguousarray(py.numpy.real(py.numpy.fft.ifft(py.numpy.multiply(py_spectral_op_1d, py_fft_1d))));
    py_vals_1d = double(py_field_1d);
    field_numpy_1d = py_vals_1d(:);
    rel_err_1d = max(abs(field_matlab_1d(:) - field_numpy_1d(:))) / max(abs(field_matlab_1d(:)));
    divergence_1d(k) = rel_err_1d;
end

fprintf('\n--- Accumulation over %d FFT round-trips (1D, N=%d) ---\n', n_roundtrips_1d, Nx_1d);
fprintf('  After   1 round-trip:  %.2e\n', divergence_1d(1));
fprintf('  After 100 round-trips: %.2e\n', divergence_1d(100));
fprintf('  After 300 round-trips: %.2e\n', divergence_1d(300));
fprintf('  After %d round-trips: %.2e\n', n_roundtrips_1d, divergence_1d(end));

%% Summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Single FFT relative diff:   1D = %.2e,  2D = %.2e\n', max_rel_1d, max_rel_2d);
fprintf('After full sim round-trips:  1D(%d) = %.2e,  2D(%d) = %.2e\n', ...
    n_roundtrips_1d, divergence_1d(end), n_roundtrips, divergence(end));
fprintf('Solver parity observed:      1D ~ 1e-15,  2D ~ 1e-13\n');
fprintf('\nIf accumulated FFT differences match solver parity, the\n');
fprintf('discrepancy is explained by FFTW vs pocketfft rounding.\n');

%% Plot
if plot_results
    figure('Position', [100 100 900 400]);

    subplot(1, 2, 1);
    semilogy(1:n_roundtrips_1d, divergence_1d, 'b-', 'LineWidth', 1.5);
    xlabel('FFT round-trips');
    ylabel('Max relative difference');
    title(sprintf('1D (N=%d): FFTW vs pocketfft', Nx_1d));
    grid on;
    yline(1e-15, 'r--', '1D solver parity', 'LineWidth', 1);
    ylim([1e-17, 1e-10]);

    subplot(1, 2, 2);
    semilogy(1:n_roundtrips, divergence, 'r-', 'LineWidth', 1.5);
    xlabel('FFT round-trips');
    ylabel('Max relative difference');
    title(sprintf('2D (%dx%d): FFTW vs pocketfft', Nx_2d, Ny_2d));
    grid on;
    yline(1e-13, 'r--', '2D solver parity', 'LineWidth', 1);
    ylim([1e-17, 1e-10]);

    sgtitle('FFT Library Divergence: MATLAB (FFTW) vs NumPy (pocketfft)');

    % Save
    saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'plots', 'fft_parity.png'));
    fprintf('\nPlot saved to tests/plots/fft_parity.png\n');
end
