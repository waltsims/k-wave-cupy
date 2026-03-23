% Debug p0 source parity — minimal case
repo_root = fileparts(fileparts(mfilename('fullpath')));
pyenv('Version', fullfile(repo_root, '.venv310', 'bin', 'python'));
addpath(fullfile(repo_root, 'k-Wave'));

Nx = 32; dx = 1;
kgrid = kWaveGrid(Nx, dx, Nx, dx);
kgrid.makeTime(1500);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = zeros(Nx, Nx);
source.p0(Nx/2+1, Nx/2+1) = 1;

sensor.mask = ones(Nx, Nx);
sensor.record = {'p'};

% MATLAB reference (no shim)
disp('=== MATLAB (native) ===');
ref = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', true, 'PMLSize', 10, 'Smooth', false);
ref_p = ref.p;
disp(['  shape: ' num2str(size(ref_p))]);
disp(['  max: ' num2str(max(abs(ref_p(:))), '%.6e')]);
disp(['  t=0 max: ' num2str(max(abs(ref_p(:,1))), '%.6e')]);
disp(['  t=1 max: ' num2str(max(abs(ref_p(:,2))), '%.6e')]);

% Python via shim
disp('=== Python (via kspaceFirstOrderPy) ===');
py_res = kspaceFirstOrderPy(kgrid, medium, source, sensor, 'PMLInside', true, 'PMLSize', 10, 'Smooth', false);
py_p = py_res.p;
disp(['  shape: ' num2str(size(py_p))]);
disp(['  max: ' num2str(max(abs(py_p(:))), '%.6e')]);
disp(['  t=0 max: ' num2str(max(abs(py_p(:,1))), '%.6e')]);
disp(['  t=1 max: ' num2str(max(abs(py_p(:,2))), '%.6e')]);

disp('=== Comparison ===');
diff = max(abs(ref_p(:) - py_p(:)));
disp(['  Max diff: ' num2str(diff, '%.6e')]);
disp(['  Rel err: ' num2str(diff/max(abs(ref_p(:))), '%.6e')]);
