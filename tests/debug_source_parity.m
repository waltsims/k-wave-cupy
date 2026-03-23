% Debug source injection parity between MATLAB and Python
% Run from k-wave-cupy root: arch -arm64 matlab -batch "run('tests/debug_source_parity.m')"

repo_root = fileparts(fileparts(mfilename('fullpath')));
pyenv('Version', fullfile(repo_root, '.venv310', 'bin', 'python'));
addpath(fullfile(repo_root, 'k-Wave'));

% Small grid
Nx = 32; dx = 1;
kgrid = kWaveGrid(Nx, dx, Nx, dx);
kgrid.setTime(5, 1/(2*1500));

medium.sound_speed = 1500;
medium.density = 1000;

% Pressure source: single row, additive
source.p_mask = zeros(Nx, Nx);
source.p_mask(1, :) = 1;
source.p = ones(1, 5);  % constant unit pressure for simplicity
source.p_mode = 'additive';

sensor.mask = ones(Nx, Nx);
sensor.record = {'p'};

% --- MATLAB reference (no shims) ---
disp('=== MATLAB reference ===');
ref = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PMLInside', true, 'PMLSize', 10, 'PMLAlpha', 2, 'Smooth', false);
ref_p = ref.p;
disp(['  max |p|: ' num2str(max(abs(ref_p(:))), '%.6e')]);
disp(['  p(1,1): ' num2str(ref_p(1,1), '%.6e')]);
disp(['  p(1,5): ' num2str(ref_p(1,5), '%.6e')]);

% --- Python via simulate_from_dicts ---
disp('=== Python ===');
module = py.importlib.import_module('kwave.solvers.kspace_solver');

kgrid_dict = py.dict(pyargs('Nx', int32(Nx), 'Ny', int32(Nx), ...
    'dx', dx, 'dy', dx, ...
    'Nt', int32(5), 'dt', kgrid.dt, ...
    'pml_size_x', int32(10), 'pml_size_y', int32(10), ...
    'pml_alpha_x', 2.0, 'pml_alpha_y', 2.0));

medium_dict = py.dict(pyargs('sound_speed', 1500.0, 'density', 1000.0));

% Build source dict
p_mask_py = py.numpy.array(source.p_mask, pyargs('dtype', 'float64', 'order', 'F'));
p_signal_py = py.numpy.array(source.p, pyargs('dtype', 'float64'));
source_dict = py.dict(pyargs('p0', py.None, ...
    'p', p_signal_py, 'p_mask', p_mask_py, 'p_mode', 'additive', ...
    'ux', py.None, 'uy', py.None, 'uz', py.None, ...
    'u_mask', py.None, 'u_mode', 'additive'));

sensor_dict = py.dict(pyargs('mask', py.numpy.ones(int32([Nx Nx]), pyargs('dtype', 'bool')), ...
    'record', py.tuple({'p'}), 'record_start_index', int32(1)));

py_result = module.simulate_from_dicts(kgrid_dict, medium_dict, source_dict, sensor_dict, pyargs('device', 'cpu'));
py_p = double(py_result{'p'});
disp(['  max |p|: ' num2str(max(abs(py_p(:))), '%.6e')]);
disp(['  p(1,1): ' num2str(py_p(1,1), '%.6e')]);
disp(['  p(1,5): ' num2str(py_p(1,5), '%.6e')]);

% --- Compare ---
max_diff = max(abs(ref_p(:) - py_p(:)));
max_val = max(abs(ref_p(:)));
disp('=== Comparison ===');
disp(['  Max abs diff: ' num2str(max_diff, '%.6e')]);
disp(['  Rel error: ' num2str(max_diff / max_val, '%.6e')]);

save('tests/debug_source_parity_result.mat', 'ref_p', 'py_p', 'kgrid', 'source', 'medium');
