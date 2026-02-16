function sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, varargin)
%KSPACEFIRSTORDER2D Shim that redirects to Python backend.
%
% This shim allows existing k-Wave tests to transparently use the
% Python/CuPy backend by replacing kspaceFirstOrder2D calls.
%
% USAGE:
%   Add tests/shims to the MATLAB path BEFORE k-Wave to enable.
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

% Redirect to Python backend
sensor_data = kspaceFirstOrderPy(kgrid, medium, source, sensor, varargin{:});

end
