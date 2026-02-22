function pass = test_1d_parity(plot_simulations, plot_comparisons)
%TEST_1D_PARITY Compare MATLAB and Python 1D PSTD stepper.

if nargin < 2, plot_comparisons = false; end
if nargin < 1, plot_simulations = false; end

pass = true;

% Check Python availability
try
    py.numpy.array(1);
catch
    warning('Python/NumPy not available. Skipping parity test.');
    return;
end

try
    % Setup
    kgrid = kWaveGrid(16, 1e-3);
    kgrid.setTime(20, 1e-7);
    
    medium.sound_speed = 1500;
    medium.density = 1000;
    
    source.p0 = zeros(kgrid.Nx, 1);
    source.p0(8) = 1;
    
    sensor.mask = ones(kgrid.Nx, 1);
    
    % Run Python
    py_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, ...
        'PMLSize', 0, 'PMLAlpha', 0, 'Smooth', false);

    % Run MATLAB
    opts = {'PlotSim', false, 'PMLAlpha', 0, 'PMLSize', 0, ...
            'PMLInside', true, 'DataCast', 'off', 'DisplayMask', 'off', ...
            'PlotPML', false, 'Smooth', false};
    mat_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, opts{:});
    if isstruct(mat_data), mat_data = mat_data.p; end
    
    % Compare
    diff = max(abs(py_data(:) - mat_data(:)));
    if diff > 1e-6
        pass = false;
        fprintf('Parity check failed. Max diff: %e\n', diff);
    end
    
    if plot_comparisons
        figure;
        subplot(3,1,1); plot(py_data); title('Python');
        subplot(3,1,2); plot(mat_data); title('MATLAB');
        subplot(3,1,3); plot(py_data - mat_data); title('Difference');
    end
    
catch ME
    disp(ME.message);
    pass = false;
end
end
