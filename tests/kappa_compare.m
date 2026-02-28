% Compare kappa values between MATLAB and Python
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'k-Wave'));

kgrid = kWaveGrid(32, 1e-4, 32, 1e-4);
kgrid.setTime(50, 1e-8);
c_ref = 1500;

% MATLAB kappa (in ifftshift/centered order first, then ifftshift to FFT order)
A = c_ref * kgrid.k * kgrid.dt / 2;
kappa_centered = sinc(A);
kappa_fft = ifftshift(kappa_centered);

% Print A values along first row (ky=0)
fprintf('--- sinc argument A = c_ref * k * dt / 2 ---\n');
fprintf('A(center) = %.6e\n', A(17, 17));  % center of ifftshifted grid = k=0
fprintf('A(1,1 centered) = %.6e  (Nyquist corner)\n', A(1, 1));
fprintf('A(17,1 centered) = %.6e  (kx=0, ky=Nyquist)\n', A(17, 1));
fprintf('A(18,17 centered) = %.6e  (first nonzero kx)\n', A(18, 17));

% kgrid.k values (centered order, k=0 at center)
fprintf('\n--- kgrid.k values ---\n');
fprintf('k(center=17,17) = %.6f\n', kgrid.k(17, 17));
fprintf('k(18,17) = %.6f  (first nonzero kx)\n', kgrid.k(18, 17));
fprintf('k(1,1) = %.6f  (corner/Nyquist)\n', kgrid.k(1, 1));

% Now compute what Python gets
py_kx = 2 * pi * (0:31) / (32 * 1e-4);
py_kx(17:end) = py_kx(17:end) - 2 * pi / 1e-4;  % wrap negatives
py_ky = py_kx;  % same for y
[KX, KY] = meshgrid(py_kx, py_ky);
KX = KX'; KY = KY';  % match indexing
py_k_mag = sqrt(KX.^2 + KY.^2);

fprintf('\n--- Python k values (FFT order) ---\n');
fprintf('k(1,1) = %.6f  (should be 0)\n', py_k_mag(1, 1));
fprintf('k(2,1) = %.6f  (first nonzero kx)\n', py_k_mag(2, 1));
fprintf('k(17,1) = %.6f  (Nyquist kx)\n', py_k_mag(17, 1));

% Python kappa (sinc(A/pi))
A_py = c_ref * py_k_mag * kgrid.dt / 2;
kappa_py = sinc(A_py / pi);  % Python convention: sinc(x/pi)

fprintf('\n--- kappa comparison (FFT order) ---\n');
fprintf('idx  | MATLAB kappa    | Python kappa    | diff\n');
fprintf('(1,1)  %.15f  %.15f  %.2e\n', kappa_fft(1,1), kappa_py(1,1), abs(kappa_fft(1,1)-kappa_py(1,1)));
fprintf('(2,1)  %.15f  %.15f  %.2e\n', kappa_fft(2,1), kappa_py(2,1), abs(kappa_fft(2,1)-kappa_py(2,1)));
fprintf('(17,1) %.15f  %.15f  %.2e\n', kappa_fft(17,1), kappa_py(17,1), abs(kappa_fft(17,1)-kappa_py(17,1)));

% Also get ACTUAL Python kappa via pyenv
py.sys.path().insert(int32(0), fullfile(fileparts(mfilename('fullpath')), '..', 'k-Wave', 'python'));
kWavePy = py.importlib.import_module('kWavePy');
py_kgrid = py.dict(pyargs('Nx', int32(32), 'dx', 1e-4, 'Ny', int32(32), 'dy', 1e-4));
py_kgrid{'Nt'} = int32(50);
py_kgrid{'dt'} = 1e-8;
py_medium = py.dict(pyargs('sound_speed', 1500.0, 'density', 1000.0));
py_source = py.dict(pyargs('p0', py.numpy.zeros(int32([32, 32]), pyargs('order', 'F'))));
py_sensor = py.dict(pyargs('mask', py.numpy.ones(int32([32, 32]), pyargs('order', 'F'))));
sim = kWavePy.Simulation(py_kgrid, py_medium, py_source, py_sensor);
sim.setup();

% Extract actual Python kappa
actual_py_kappa = double(sim.kappa.flatten(pyargs('order', 'F')));
actual_py_kappa = reshape(actual_py_kappa, [32, 32]);

fprintf('\n--- ACTUAL Python kappa from sim.kappa ---\n');
fprintf('(1,1)  %.15f  vs MATLAB %.15f  diff=%.2e\n', actual_py_kappa(1,1), kappa_fft(1,1), abs(actual_py_kappa(1,1)-kappa_fft(1,1)));
fprintf('(2,1)  %.15f  vs MATLAB %.15f  diff=%.2e\n', actual_py_kappa(2,1), kappa_fft(2,1), abs(actual_py_kappa(2,1)-kappa_fft(2,1)));
fprintf('(17,1) %.15f  vs MATLAB %.15f  diff=%.2e\n', actual_py_kappa(17,1), kappa_fft(17,1), abs(actual_py_kappa(17,1)-kappa_fft(17,1)));

max_kappa_diff = max(abs(actual_py_kappa(:) - kappa_fft(:)));
fprintf('\nMax kappa difference: %.2e\n', max_kappa_diff);
if max_kappa_diff > 1e-10
    fprintf('*** KAPPA VALUES DIFFER SIGNIFICANTLY - this is likely the source of the 2D parity gap ***\n');
else
    fprintf('Kappa values match to machine precision.\n');
end
