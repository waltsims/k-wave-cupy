%DIAGNOSE_2D_PARITY Step-by-step comparison of Python vs MATLAB for 2D.
%
% This script identifies WHERE the ~10% error between Python and MATLAB
% originates by comparing at increasing time steps.

clearvars;
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'k-Wave'));

% Setup Python
pyenv('Version', fullfile(project_root, '.venv310', 'bin', 'python'));

% Ensure kWavePy is reloaded
module_dir = fullfile(project_root, 'k-Wave', 'python');
if ~any(strcmp(cell(py.sys.path), module_dir))
    insert(py.sys.path, int32(0), module_dir);
end
kWavePy = py.importlib.import_module('kWavePy');
py.importlib.reload(kWavePy);

%% Setup identical simulation parameters
Nx = 32; Ny = 32;
dx = 0.1e-3; dy = 0.1e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1500;
medium.density = 1000;

% Simple point source at center
source.p0 = zeros(Nx, Ny);
source.p0(Nx/2, Ny/2) = 1;

sensor.mask = ones(Nx, Ny);

% Short simulation for diagnosis
kgrid.makeTime(medium.sound_speed);
Nt_test = 10;  % Just 10 steps for diagnosis
kgrid.setTime(Nt_test, kgrid.dt);

fprintf('Grid: %dx%d, dt=%.2e, Nt=%d\n', Nx, Ny, kgrid.dt, Nt_test);

%% Run MATLAB reference (full simulation, extract final pressure)
sensor_data_matlab = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PMLInside', false, 'PMLSize', 0, 'PlotSim', false, 'Smooth', false, 'DataCast', 'off');

p_matlab_final = reshape(sensor_data_matlab(:, end), [Nx, Ny]);
fprintf('MATLAB final pressure: max=%.6e, sum=%.6e\n', max(p_matlab_final(:)), sum(p_matlab_final(:)));

%% Run Python step-by-step using Simulation class
toNumpy = @(x) py.numpy.array(double(x), pyargs('order', 'F'));

k_py = py.dict(pyargs('Nx', int64(Nx), 'dx', dx, 'Ny', int64(Ny), 'dy', dy, ...
    'Nt', int64(Nt_test), 'dt', kgrid.dt));
m_py = py.dict(pyargs('sound_speed', toNumpy(medium.sound_speed), ...
    'density', toNumpy(medium.density)));
s_py = py.dict(pyargs('p0', toNumpy(source.p0)));
d_py = py.dict(pyargs('mask', toNumpy(sensor.mask)));

% Create Simulation object using helper function
sim = kWavePy.create_simulation(k_py, m_py, s_py, d_py, 'numpy');
sim.setup();

% Check setup values using py.getattr
fprintf('\n=== SETUP COMPARISON ===\n');
fprintf('Python c_ref: %.2f\n', double(py.getattr(sim, 'c_ref')));
fprintf('Python dt: %.2e\n', double(py.getattr(sim, 'dt')));
fprintf('Python ndim: %d\n', double(py.getattr(sim, 'ndim')));

% Get initial p0
p_py_t0 = double(py.numpy.array(py.getattr(sim, 'p')));
fprintf('Python p0 (before any steps): max=%.6e\n', max(p_py_t0(:)));

%% Step through and compare at each step
fprintf('\n=== STEP-BY-STEP COMPARISON ===\n');
fprintf('Step | Python max | MATLAB max | Rel Error\n');
fprintf('-----|------------|------------|----------\n');

for t = 1:Nt_test
    sim.step();

    p_py = double(py.numpy.array(py.getattr(sim, 'p')));
    p_matlab = reshape(sensor_data_matlab(:, t), [Nx, Ny]);

    % Compute relative error
    diff_abs = abs(p_py - p_matlab);
    max_diff = max(diff_abs(:));
    rel_err = max_diff / max(abs(p_matlab(:)) + 1e-15);

    fprintf('%4d | %10.4e | %10.4e | %8.2e\n', ...
        t, max(abs(p_py(:))), max(abs(p_matlab(:))), rel_err);
end

%% Visualize final comparison
figure('Position', [100 100 1200 400]);

subplot(1,4,1);
imagesc(p_matlab_final);
title('MATLAB Final');
colorbar; axis image;

p_py_final = double(py.numpy.array(py.getattr(sim, 'p')));
subplot(1,4,2);
imagesc(p_py_final);
title('Python Final');
colorbar; axis image;

subplot(1,4,3);
imagesc(p_matlab_final - p_py_final);
title('Difference');
colorbar; axis image;

subplot(1,4,4);
plot(p_matlab_final(Nx/2, :), 'b-', 'LineWidth', 1.5);
hold on;
plot(p_py_final(Nx/2, :), 'r--', 'LineWidth', 1.5);
legend('MATLAB', 'Python');
title('Cross-section at y=Ny/2');

plot_path = fullfile(project_root, 'tests', 'plots', 'diagnose_2d_parity.png');
saveas(gcf, plot_path);
fprintf('\nSaved: %s\n', plot_path);
