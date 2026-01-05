function tests = test_interface_1D
%TEST_INTERFACE_1D Smoke-test the Python 1D wrapper.
tests = functiontests(localfunctions);
end

function testSimpleCallReturnsData(testCase)
assumePythonAvailable(testCase);

kgrid = kWaveGrid(8, 1e-3);
kgrid.setTime(6, 1e-7);

medium.sound_speed = 1500;
medium.density = 1000;

source.p0 = zeros(kgrid.Nx, 1);
source.p0(4) = 1;

sensor.mask = ones(kgrid.Nx, 1);

data = kspaceFirstOrderPy(kgrid, medium, source, sensor);

testCase.verifySize(data, [nnz(sensor.mask), kgrid.Nt]);
testCase.verifyFalse(any(isnan(data(:))), 'Python backend returned NaNs');
end

function assumePythonAvailable(testCase)
try
    env = pyenv;
    % trigger lazy load if needed
    if strcmp(env.Status, "NotLoaded")
        py.list();
    end
catch ME
    testCase.assumeFail("Python unavailable: " + ME.message);
end
end
