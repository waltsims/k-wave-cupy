function sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDERPY Minimal Python-backed 1D k-Wave solver.
%
% sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor)
%
% Only 1D grids are supported. Call kgrid.makeTime or setTime before
% invoking this wrapper.

backend = parseBackend(varargin{:});
assert1DGrid(kgrid);
assertTimeIsSpecified(kgrid);

module = importPythonModule();

kgrid_py = py.dict(pyargs( ...
    'Nx', int64(kgrid.Nx), ...
    'dx', kgrid.dx, ...
    'Nt', int64(kgrid.Nt), ...
    'dt', kgrid.dt));

medium_py = py.dict(pyargs( ...
    'sound_speed', toNumpy(fetchField(medium, {'sound_speed', 'c0'})), ...
    'density', toNumpy(fetchField(medium, {'density', 'rho0'}, 1000))));

source_py = py.dict(pyargs( ...
    'p0', toNumpy(fetchField(source, {'p0'}, zeros(kgrid.Nx, 1)))));

sensor_py = py.dict(pyargs( ...
    'mask', toNumpy(fetchField(sensor, {'mask'}, ones(kgrid.Nx, 1)))));

result = module.simulate(kgrid_py, medium_py, source_py, sensor_py, pyargs('backend', backend));

sensor_data = double(result{'sensor_data'});
end

% -------------------------------------------------------------------------
function module = importPythonModule()
module_dir = fullfile(fileparts(mfilename('fullpath')), 'python');

% ensure the python folder is on sys.path
paths = cell(py.sys.path);
if ~any(strcmp(paths, module_dir))
    insert(py.sys.path, int32(0), module_dir);
end

module = py.importlib.import_module('kWavePy');
py.importlib.reload(module);
end

function data = toNumpy(array)
data = py.numpy.array(double(array), pyargs('order', 'F'));
end

function value = fetchField(structLike, names, default)
if nargin < 3
    default = [];
end

value = default;
for idx = 1:numel(names)
    name = names{idx};
    if isstruct(structLike) && isfield(structLike, name)
        value = structLike.(name);
        return;
    end
    if isprop(structLike, name)
        value = structLike.(name);
        return;
    end
end
end

function assert1DGrid(kgrid)
if isprop(kgrid, 'dim') && kgrid.dim ~= 1
    error('kspaceFirstOrderPy:dim', 'Only 1D grids are supported by kspaceFirstOrderPy.');
end
if isprop(kgrid, 'Ny') && kgrid.Ny > 1
    error('kspaceFirstOrderPy:dim', 'Only 1D grids are supported by kspaceFirstOrderPy.');
end
if isprop(kgrid, 'Nz') && kgrid.Nz > 1
    error('kspaceFirstOrderPy:dim', 'Only 1D grids are supported by kspaceFirstOrderPy.');
end
end

function assertTimeIsSpecified(kgrid)
dt = kgrid.dt;
if ~isnumeric(dt) || isempty(dt) || any(dt <= 0) ...
        || ischar(dt) || (isstring(dt) && (strlength(dt) == 0 || dt == "auto"))
    error('kspaceFirstOrderPy:time', 'kgrid.dt is not set. Call kgrid.makeTime or setTime first.');
end
end

function backend = parseBackend(varargin)
backend = 'auto';
if isempty(varargin)
    return;
end

parser = inputParser;
addParameter(parser, 'Backend', backend);
parse(parser, varargin{:});

backend = parser.Results.Backend;
if isstring(backend)
    backend = char(backend);
end
end
