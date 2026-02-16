function sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDERPY Minimal Python-backed 1D k-Wave solver.

% Validate inputs (Python backend has stricter requirements than MATLAB)
if kgrid.dim ~= 1, error('Only 1D grids supported.'); end
if isempty(kgrid.dt) || (ischar(kgrid.dt) && strcmp(kgrid.dt, 'auto'))
    error('kspaceFirstOrderPy:TimeNotSet', ...
        ['Time array not set. The Python backend requires explicit time array.\n' ...
         'Use one of:\n' ...
         '  kgrid.makeTime(medium.sound_speed)      %% Auto-calculate based on CFL\n' ...
         '  kgrid.setTime(Nt, dt)                   %% Set explicitly']);
end
p = inputParser; p.KeepUnmatched = true; addParameter(p, 'Backend', 'auto'); parse(p, varargin{:});

% Load Python module once per session (persistent avoids repeated imports)
persistent kWavePy
if isempty(kWavePy)
    module_dir = fullfile(fileparts(mfilename('fullpath')), 'python');
    if ~any(strcmp(cell(py.sys.path), module_dir)), insert(py.sys.path, int32(0), module_dir); end
    kWavePy = py.importlib.import_module('kWavePy');
    py.importlib.reload(kWavePy);
end

% Marshal MATLAB structs to Python dicts (column-major order for MATLAB compatibility)
toNumpy = @(x) py.numpy.array(double(x), pyargs('order', 'F'));
getField = @(s, names, default) getFieldValue(s, names, default);

k_py = py.dict(pyargs('Nx', int64(kgrid.Nx), 'dx', kgrid.dx, 'Nt', int64(kgrid.Nt), 'dt', kgrid.dt));
m_py = py.dict(pyargs( ...
    'sound_speed', toNumpy(getField(medium, {'sound_speed','c0'}, [])), ...
    'density',     toNumpy(getField(medium, {'density','rho0'}, 1000)), ...
    'alpha_coeff', toNumpy(getField(medium, {'alpha_coeff'}, 0)), ...
    'alpha_power', toNumpy(getField(medium, {'alpha_power'}, 1.5)), ...
    'BonA',        toNumpy(getField(medium, {'BonA'}, 0))));
s_py = py.dict(pyargs( ...
    'p0',     toNumpy(getField(source, {'p0'}, 0)), ...
    'p_mask', toNumpy(getField(source, {'p_mask'}, 0)), ...
    'p',      toNumpy(getField(source, {'p'}, 0)), ...
    'p_mode', getField(source, {'p_mode'}, 'additive'), ...
    'u_mask', toNumpy(getField(source, {'u_mask'}, 0)), ...
    'ux',     toNumpy(getField(source, {'ux'}, 0)), ...
    'u_mode', getField(source, {'u_mode'}, 'additive')));
d_py = py.dict(pyargs('mask', toNumpy(getField(sensor, {'mask'}, 1))));

% Run simulation and convert result back to MATLAB double
res = kWavePy.simulate_from_dicts(k_py, m_py, s_py, d_py, pyargs('backend', char(p.Results.Backend)));
sensor_data = double(res{'sensor_data'});
end

function value = getFieldValue(s, names, default)
    % Try each field name in order (supports aliasing like 'sound_speed'/'c0')
    if ischar(names), names = {names}; end
    for i = 1:numel(names)
        if isprop(s, names{i}) || (isstruct(s) && isfield(s, names{i}))
            value = s.(names{i});
            return;
        end
    end
    value = default;
end
