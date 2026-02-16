function sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDERPY Minimal Python-backed N-D k-Wave solver (1D, 2D, 3D).

% Validate inputs (Python backend requires explicit time array)
if kgrid.dim > 3, error('kspaceFirstOrderPy:UnsupportedDimension', 'Only 1D, 2D, 3D supported.'); end
if isempty(kgrid.dt) || (ischar(kgrid.dt) && strcmp(kgrid.dt, 'auto'))
    error('kspaceFirstOrderPy:TimeNotSet', ...
        ['Time array not set. The Python backend requires explicit time array.\n' ...
         'Use one of:\n' ...
         '  kgrid.makeTime(medium.sound_speed)      %% Auto-calculate based on CFL\n' ...
         '  kgrid.setTime(Nt, dt)                   %% Set explicitly']);
end
p = inputParser; p.KeepUnmatched = true;
addParameter(p, 'Backend', 'auto');
addParameter(p, 'PMLSize', 20);
addParameter(p, 'PMLAlpha', 2);
parse(p, varargin{:});

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

% Build kgrid dict with conditional fields for each dimension
% Handle PMLSize as scalar or [x,y] or [x,y,z] array
pml_size = p.Results.PMLSize;
pml_alpha = p.Results.PMLAlpha;
if isscalar(pml_size), pml_size = repmat(pml_size, 1, kgrid.dim); end
if isscalar(pml_alpha), pml_alpha = repmat(pml_alpha, 1, kgrid.dim); end

kgrid_args = {'Nx', int64(kgrid.Nx), 'dx', kgrid.dx, 'Nt', int64(kgrid.Nt), 'dt', kgrid.dt, ...
    'pml_size_x', int64(pml_size(1)), 'pml_alpha_x', pml_alpha(1)};
if kgrid.dim >= 2
    kgrid_args = [kgrid_args, {'Ny', int64(kgrid.Ny), 'dy', kgrid.dy, ...
        'pml_size_y', int64(pml_size(2)), 'pml_alpha_y', pml_alpha(2)}];
end
if kgrid.dim >= 3
    kgrid_args = [kgrid_args, {'Nz', int64(kgrid.Nz), 'dz', kgrid.dz, ...
        'pml_size_z', int64(pml_size(3)), 'pml_alpha_z', pml_alpha(3)}];
end
k_py = py.dict(pyargs(kgrid_args{:}));

m_py = py.dict(pyargs( ...
    'sound_speed', toNumpy(getField(medium, {'sound_speed','c0'}, [])), ...
    'density',     toNumpy(getField(medium, {'density','rho0'}, 1000)), ...
    'alpha_coeff', toNumpy(getField(medium, {'alpha_coeff'}, 0)), ...
    'alpha_power', toNumpy(getField(medium, {'alpha_power'}, 1.5)), ...
    'BonA',        toNumpy(getField(medium, {'BonA'}, 0))));

% Build source dict with velocity components for each dimension
source_args = { ...
    'p0',     toNumpy(getField(source, {'p0'}, 0)), ...
    'p_mask', toNumpy(getField(source, {'p_mask'}, 0)), ...
    'p',      toNumpy(getField(source, {'p'}, 0)), ...
    'p_mode', getField(source, {'p_mode'}, 'additive'), ...
    'u_mask', toNumpy(getField(source, {'u_mask'}, 0)), ...
    'ux',     toNumpy(getField(source, {'ux'}, 0)), ...
    'u_mode', getField(source, {'u_mode'}, 'additive')};
if kgrid.dim >= 2
    source_args = [source_args, {'uy', toNumpy(getField(source, {'uy'}, 0))}];
end
if kgrid.dim >= 3
    source_args = [source_args, {'uz', toNumpy(getField(source, {'uz'}, 0))}];
end
s_py = py.dict(pyargs(source_args{:}));

d_py = py.dict(pyargs('mask', toNumpy(getField(sensor, {'mask'}, 1)), ...
    'record_start_index', int64(getField(sensor, {'record_start_index'}, 1))));

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
