% BENCHMARK_2D  Wall-clock comparison of MATLAB vs Python/NumPy on 2D grids.
%
% Prints a table and saves results to tests/benchmark_2d_results.mat.
%
% Usage:
%   pyenv('Version', fullfile(pwd,'.venv310','bin','python'));
%   addpath('k-Wave'); addpath('tests');
%   benchmark_2d

% grid sizes and matching time steps
grid_sizes = [64, 128, 256, 512];
Nt_list    = round(grid_sizes * 1.5);

% fixed parameters
pml_size = 10;
cfl      = 0.1;
c0       = 1500;   % [m/s]
rho0     = 1000;   % [kg/m^3]

% preallocate results
n = numel(grid_sizes);
t_matlab = zeros(n, 1);
t_python = zeros(n, 1);

% check Python availability
try
    env = pyenv;
    if strcmp(env.Status, "NotLoaded")
        py.list();
    end
    python_ok = true;
catch
    warning('Python not available — only MATLAB will be benchmarked.');
    python_ok = false;
end

for ii = 1:n
    Nx = grid_sizes(ii);
    Ny = Nx;
    Nt = Nt_list(ii);
    dx = 1e-3;

    % create grid
    kgrid = kWaveGrid(Nx, dx, Ny, dx);
    kgrid.makeTime(c0, cfl, Nt * cfl * dx / c0);

    % medium
    medium.sound_speed = c0;
    medium.density     = rho0;

    % point source at centre
    source.p_mask = zeros(Nx, Ny);
    source.p_mask(round(Nx/2), round(Ny/2)) = 1;
    source.p = sin(2 * pi * 1e6 / c0 * kgrid.t_array);

    % single-point sensor
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(round(Nx/4), round(Ny/4)) = 1;

    % --- MATLAB backend ---
    tic;
    kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
        'PlotSim', false, 'PMLSize', pml_size);
    t_matlab(ii) = toc;

    % --- Python backend ---
    if python_ok
        tic;
        kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
            'PMLSize', pml_size);
        t_python(ii) = toc;
    end
end

% print table
fprintf('\n%-12s %-6s %-12s %-12s %-8s\n', ...
    'Grid', 'Nt', 'MATLAB(s)', 'Python(s)', 'Ratio');
fprintf('%s\n', repmat('-', 1, 54));
for ii = 1:n
    Nx = grid_sizes(ii);
    Nt = Nt_list(ii);
    if python_ok
        ratio = t_matlab(ii) / t_python(ii);
        fprintf('%-12s %-6d %-12.3f %-12.3f %-8.2f\n', ...
            sprintf('%dx%d', Nx, Nx), Nt, t_matlab(ii), t_python(ii), ratio);
    else
        fprintf('%-12s %-6d %-12.3f %-12s %-8s\n', ...
            sprintf('%dx%d', Nx, Nx), Nt, t_matlab(ii), 'N/A', 'N/A');
    end
end

% save results
results.grid_sizes = grid_sizes;
results.Nt_list    = Nt_list;
results.t_matlab   = t_matlab;
results.t_python   = t_python;
if python_ok
    results.ratio = t_matlab ./ t_python;
end

save(fullfile(fileparts(mfilename('fullpath')), 'benchmark_2d_results.mat'), ...
    '-struct', 'results');
fprintf('\nResults saved to tests/benchmark_2d_results.mat\n');
