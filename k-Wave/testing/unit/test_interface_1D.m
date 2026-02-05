function test_pass = test_interface_1D(plot_simulations, plot_comparisons)
%TEST_INTERFACE_1D Smoke-test the Python 1D wrapper.

if nargin < 2, plot_comparisons = false; end
if nargin < 1, plot_simulations = false; end

test_pass = true;

% Check Python availability
try
    env = pyenv;
    % trigger lazy load if needed
    if strcmp(env.Status, "NotLoaded")
        py.list();
    end
catch ME
    warning('Python unavailable: %s. Skipping test.', ME.message);
    return;
end

try
    % Setup test
    kgrid = kWaveGrid(8, 1e-3);
    kgrid.setTime(6, 1e-7);

    medium.sound_speed = 1500;
    medium.density = 1000;

    source.p0 = zeros(kgrid.Nx, 1);
    source.p0(4) = 1;

    sensor.mask = ones(kgrid.Nx, 1);

    % Run simulation
    data = kspaceFirstOrderPy(kgrid, medium, source, sensor);

    % Verify results
    expected_size = [nnz(sensor.mask), kgrid.Nt];
    if ~isequal(size(data), expected_size)
        fprintf('FAIL: Expected size [%d, %d], got [%d, %d]\n', ...
            expected_size(1), expected_size(2), size(data, 1), size(data, 2));
        test_pass = false;
        return;
    end

    if any(isnan(data(:)))
        fprintf('FAIL: Python backend returned NaNs\n');
        test_pass = false;
        return;
    end

    fprintf('PASS: Python 1D interface test completed successfully\n');

catch ME
    fprintf('FAIL: %s\n', ME.message);
    test_pass = false;
end
end
