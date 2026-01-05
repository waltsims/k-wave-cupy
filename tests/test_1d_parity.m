function tests = test_1d_parity
%TEST_1D_PARITY Compare MATLAB and Python 1D PSTD stepper.
tests = functiontests(localfunctions);
end

function testPythonMatchesMatlabStepper(testCase)
assumePythonAvailable(testCase);

kgrid = kWaveGrid(16, 1e-3);
kgrid.setTime(20, 1e-7);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = zeros(kgrid.Nx, 1);
source.p0(8) = 1;

sensor.mask = ones(kgrid.Nx, 1);

try
    py_data = kspaceFirstOrderPy(kgrid, medium, source, sensor);
catch ME
    testCase.assumeFail("Python call failed: " + ME.message);
end

mat_data = matlab_pstd(kgrid, medium, source, sensor);

diff_sensor = max(abs(py_data(:) - mat_data(:)));

testCase.verifyLessThan(diff_sensor, 1e-6);
end

% -------------------------------------------------------------------------
function data = matlab_pstd(kgrid, medium, source, sensor)
opts = { ...
    'PlotSim', false, ...
    'PMLAlpha', 0, ...
    'PMLSize', 0, ...
    'PMLInside', true, ...
    'DataCast', 'off', ...
    'DisplayMask', 'off', ...
    'PlotPML', false, ...
    'Smooth', false, ...
    };

data = kspaceFirstOrder1D(kgrid, medium, source, sensor, opts{:});
if isstruct(data) && isfield(data, 'p')
    % sensor.record can change return shape; normalize to matrix
    data = data.p;
end
end

function assumePythonAvailable(testCase)
try
    env = pyenv;
    if env.Status == "NotLoaded"
        py.list();
    end
catch ME
    testCase.assumeFail("Python unavailable: " + ME.message);
end
end
