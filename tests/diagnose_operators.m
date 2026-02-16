%DIAGNOSE_OPERATORS Compare k-space operators between Python and MATLAB.

clearvars;
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'k-Wave'));

% Setup Python
pyenv('Version', fullfile(project_root, '.venv310', 'bin', 'python'));
module_dir = fullfile(project_root, 'k-Wave', 'python');
if ~any(strcmp(cell(py.sys.path), module_dir))
    insert(py.sys.path, int32(0), module_dir);
end

%% Simple 2D grid
Nx = 8; Ny = 8;
dx = 0.1e-3; dy = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
c_ref = 1500;
dt = 2e-8;

fprintf('=== K-VECTOR COMPARISON (Nx=%d) ===\n', Nx);

% MATLAB k-vectors (centered order)
fprintf('\nMATLAB kx_vec (centered):\n');
fprintf('  %.4f', kgrid.kx_vec'); fprintf('\n');

% MATLAB k-vectors after ifftshift (FFT order)
kx_fft_order = ifftshift(kgrid.kx_vec);
fprintf('MATLAB ifftshift(kx_vec) (FFT order):\n');
fprintf('  %.4f', kx_fft_order'); fprintf('\n');

% Python k-vector (FFT order via fftfreq)
k_py = 2 * pi * double(py.numpy.fft.fftfreq(int64(Nx), dx));
fprintf('Python fftfreq k (FFT order):\n');
fprintf('  %.4f', k_py); fprintf('\n');

% Check if they match
diff_k = max(abs(kx_fft_order - k_py'));
fprintf('\nMax diff between MATLAB and Python k-vectors: %.2e\n', diff_k);

%% Compare sinc values
fprintf('\n=== SINC/KAPPA COMPARISON ===\n');

% MATLAB sinc (note: MATLAB sinc(x) = sin(pi*x)/(pi*x), so sinc(x/pi) = sin(x)/x)
kappa_matlab = sinc(c_ref * kgrid.kx_vec * dt / 2 / pi);
kappa_matlab_fft = ifftshift(kappa_matlab);
fprintf('MATLAB kappa (centered): '); fprintf('%.6f ', kappa_matlab'); fprintf('\n');
fprintf('MATLAB kappa (FFT order): '); fprintf('%.6f ', kappa_matlab_fft'); fprintf('\n');

% Python: np.sinc(x) = sin(pi*x)/(pi*x), same as MATLAB
% Our code: xp.sinc((c_ref * k * dt / 2) / np.pi)
% This means sinc(arg/pi) where arg = c_ref * k * dt / 2
arg_py = c_ref * k_py * dt / 2;
kappa_py = double(py.numpy.sinc(arg_py / pi));
fprintf('Python kappa (FFT order): '); fprintf('%.6f ', kappa_py); fprintf('\n');

diff_kappa = max(abs(kappa_matlab_fft - kappa_py'));
fprintf('Max diff in kappa: %.2e\n', diff_kappa);

%% Compare full gradient operator
fprintf('\n=== GRADIENT OPERATOR COMPARISON ===\n');

% MATLAB: 1i * k * exp(1i * k * dx/2) * kappa
op_grad_matlab = ifftshift(1i * kgrid.kx_vec .* exp(1i * kgrid.kx_vec * dx/2) .* kappa_matlab);
fprintf('MATLAB op_grad (real): '); fprintf('%.6f ', real(op_grad_matlab')); fprintf('\n');
fprintf('MATLAB op_grad (imag): '); fprintf('%.6f ', imag(op_grad_matlab')); fprintf('\n');

% Python: 1j * k * kappa * exp(1j * k * dx/2)
op_grad_py = 1j * k_py .* kappa_py .* exp(1j * k_py * dx/2);
fprintf('Python op_grad (real): '); fprintf('%.6f ', real(op_grad_py)); fprintf('\n');
fprintf('Python op_grad (imag): '); fprintf('%.6f ', imag(op_grad_py)); fprintf('\n');

diff_op = max(abs(op_grad_matlab - op_grad_py'));
fprintf('Max diff in gradient operator: %.2e\n', diff_op);

%% Test with actual field
fprintf('\n=== GRADIENT APPLICATION TEST ===\n');

p_test = zeros(Nx, Ny);
p_test(Nx/2, Ny/2) = 1;

% MATLAB: bsxfun for broadcasting Nx x 1 operator over Nx x Ny
p_k = fft2(p_test);
% Full MATLAB operator includes kappa separately in some cases
grad_p_matlab = real(ifft2(bsxfun(@times, op_grad_matlab, p_k)));

% Python approach: op already includes kappa, reshape to Nx x 1
op_grad_2d = reshape(op_grad_py, [Nx, 1]);
grad_p_py = real(ifft2(op_grad_2d .* p_k));

fprintf('MATLAB grad at center row: '); fprintf('%.6f ', grad_p_matlab(Nx/2, :)); fprintf('\n');
fprintf('Python grad at center row: '); fprintf('%.6f ', grad_p_py(Nx/2, :)); fprintf('\n');

diff_grad = max(abs(grad_p_matlab(:) - grad_p_py(:)));
fprintf('Max diff in gradient field: %.2e\n', diff_grad);

fprintf('\n=== SUMMARY ===\n');
if diff_k < 1e-14 && diff_kappa < 1e-14 && diff_op < 1e-14 && diff_grad < 1e-14
    fprintf('All operators match within machine precision!\n');
else
    fprintf('Differences found - investigate further.\n');
end
