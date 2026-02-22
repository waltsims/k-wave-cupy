function [pressure_max, pressure_time] = angularSpectrumPy(input_plane, dx, dt, z_pos, medium, varargin)
%ANGULARSPECTRUMPY Python-backed angular spectrum propagator.
%
% Mirrors the interface of angularSpectrum.m but delegates computation
% to the Python/NumPy backend in kWavePy.angular_spectrum.
%
% USAGE:
%     pressure_max = angularSpectrumPy(input_plane, dx, dt, z_pos, c0)
%     [pressure_max, pressure_time] = angularSpectrumPy(input_plane, dx, dt, z_pos, medium, ...)
%
% See also angularSpectrum

% Parse medium
if isstruct(medium)
    c0 = medium.sound_speed;
    alpha_coeff = getFieldDefault(medium, 'alpha_coeff', 0);
    alpha_power = getFieldDefault(medium, 'alpha_power', 1.5);
else
    c0 = medium;
    alpha_coeff = 0;
    alpha_power = 1.5;
end

% Parse optional arguments
p = inputParser;
addParameter(p, 'AngularRestriction', true);
addParameter(p, 'GridExpansion', 0);
addParameter(p, 'Reverse', false);
parse(p, varargin{:});

% Load Python module
persistent kWavePy
if isempty(kWavePy)
    module_dir = fullfile(fileparts(mfilename('fullpath')), 'python');
    if ~any(strcmp(cell(py.sys.path), module_dir)), insert(py.sys.path, int32(0), module_dir); end
    kWavePy = py.importlib.import_module('kWavePy');
    py.importlib.reload(kWavePy);
end

toNumpy = @(x) py.numpy.array(double(x), pyargs('order', 'F'));

% Call Python
res = kWavePy.angular_spectrum( ...
    toNumpy(input_plane), dx, dt, toNumpy(z_pos), c0, ...
    pyargs('alpha_coeff', alpha_coeff, 'alpha_power', alpha_power, ...
           'angular_restriction', logical(p.Results.AngularRestriction), ...
           'grid_expansion', int64(p.Results.GridExpansion), ...
           'reverse', logical(p.Results.Reverse)));

pressure_max = double(res{1});
if nargout >= 2
    pressure_time = double(res{2});
end
end

function val = getFieldDefault(s, name, default)
    if isfield(s, name)
        val = s.(name);
    else
        val = default;
    end
end
