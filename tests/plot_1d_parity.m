%PLOT_1D_PARITY Visualise MATLAB vs Python 1D PSTD outputs.
%
% Saves a plot to tests/plots/1d_parity.png comparing sensor data from
% kspaceFirstOrder1D (MATLAB) and kspaceFirstOrderPy (Python backend).

repo = fileparts(fileparts(mfilename('fullpath')));
plots_dir = fullfile(repo, 'tests', 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

% Configure Python interpreter
try
    pyenv('Version', fullfile(repo, '.venv310', 'bin', 'python'));
catch
    % keep going if already loaded
end

addpath(fullfile(repo, 'k-Wave'));
addpath(fullfile(repo, 'tests'));
addpath(fullfile(repo, 'k-Wave', 'testing', 'unit'));

kgrid = kWaveGrid(16, 1e-3);
kgrid.setTime(20, 1e-7);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = zeros(kgrid.Nx, 1);
source.p0(8) = 1;

sensor.mask = ones(kgrid.Nx, 1);

py_data = kspaceFirstOrderPy(kgrid, medium, source, sensor);

opts = { ...
    'PlotSim', false, ...
    'PMLAlpha', 0, ...
    'PMLSize', 0, ...
    'PMLInside', true, ...
    'DataCast', 'off', ...
    'DisplayMask', 'off', ...
    'PlotPML', false, ...
    'Smooth', false, ...
    'RecordMovie', false ...
    };

mat_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, opts{:});
if isstruct(mat_data) && isfield(mat_data, 'p')
    mat_data = mat_data.p;
end

diff_data = py_data - mat_data;
max_diff = max(abs(diff_data(:)));
fprintf('Max abs diff: %g\n', max_diff);

f = figure('Name', '1D Parity', 'Visible', 'off');
subplot(3, 1, 1);
imagesc(py_data);
title('Python kspaceFirstOrderPy sensor data');
xlabel('time step');
ylabel('sensor index');
colorbar;

subplot(3, 1, 2);
imagesc(mat_data);
title('MATLAB kspaceFirstOrder1D sensor data');
xlabel('time step');
ylabel('sensor index');
colorbar;

subplot(3, 1, 3);
imagesc(diff_data);
title(sprintf('Difference (Python - MATLAB), max |diff| = %.3g', max_diff));
xlabel('time step');
ylabel('sensor index');
colorbar;

saveas(f, fullfile(plots_dir, '1d_parity.png'));
close(f);
