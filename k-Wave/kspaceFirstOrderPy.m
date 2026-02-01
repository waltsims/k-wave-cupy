function sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDERPY Minimal Python-backed 1D k-Wave solver.

% 1. Validation & Parsing
if kgrid.dim ~= 1, error('Only 1D grids supported.'); end
if isempty(kgrid.dt), error('Time step dt not set.'); end
p = inputParser; addParameter(p, 'Backend', 'auto'); parse(p, varargin{:});

% 2. Setup Python Environment
module_dir = fullfile(fileparts(mfilename('fullpath')), 'python');
if ~any(strcmp(cell(py.sys.path), module_dir))
    insert(py.sys.path, int32(0), module_dir);
end
kWavePy = py.importlib.import_module('kWavePy');
py.importlib.reload(kWavePy);

% 3. Marshal Data
toPy = @(x) py.numpy.array(castForPy(x), pyargs('order', 'F'));
val  = @(s,n,d) toPy(getField(s, n, d));

k_py = py.dict(pyargs('Nx',int64(kgrid.Nx), 'dx',kgrid.dx, 'Nt',int64(kgrid.Nt), 'dt',kgrid.dt));
m_py = py.dict(pyargs('sound_speed', val(medium,{'sound_speed','c0'},[]), 'density', val(medium,{'density','rho0'},1000)));
s_py = py.dict(pyargs('p0', val(source,{'p0'},0)));
d_py = py.dict(pyargs('mask', val(sensor,{'mask'},1)));

% 4. Run Simulation
res = kWavePy.simulate(k_py, m_py, s_py, d_py, pyargs('backend', char(p.Results.Backend)));
sensor_data = double(res{'sensor_data'});
end

function y = castForPy(x)
    if isa(x, 'single')
        y = x;
    else
        y = double(x);
    end
end

function v = getField(s, names, default)
    v = default;
    for i=1:numel(names)
        if isprop(s, names{i}) || (isstruct(s) && isfield(s, names{i}))
            v = s.(names{i}); return;
        end
    end
end
