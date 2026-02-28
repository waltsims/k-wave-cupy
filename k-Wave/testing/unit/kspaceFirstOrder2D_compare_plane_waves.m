function test_pass = kspaceFirstOrder2D_compare_plane_waves(plot_comparisons, plot_simulations)
% DESCRIPTION:
%     Unit test to validate 2D wave propagation with comprehensive
%     feature coverage including absorption, nonlinearity, and all source
%     types. This is the 2D version of
%     kspaceFirstOrder1D_compare_plane_waves.m
%
%     When run with standard MATLAB backend, generates reference data.
%     When run with Python backend (via shim), compares against reference.
%
%     90 tests are performed:
%
%         1.  linear + lossless + source.p0 + homogeneous
%         2.  linear + lossless + source.p0 + heterogeneous
%         3.  linear + lossless + source.p (additive) + homogeneous
%         4.  linear + lossless + source.p (additive) + heterogeneous
%         5.  linear + lossless + source.p (additive-no-correction) + homogeneous
%         6.  linear + lossless + source.p (additive-no-correction) + heterogeneous
%         7.  linear + lossless + source.p (dirichlet) + homogeneous
%         8.  linear + lossless + source.p (dirichlet) + heterogeneous
%         9.  linear + lossless + source.ux (additive) + homogeneous
%         10. linear + lossless + source.ux (additive) + heterogeneous
%         11. linear + lossless + source.ux (additive-no-correction) + homogeneous
%         12. linear + lossless + source.ux (additive-no-correction) + heterogeneous
%         13. linear + lossless + source.ux (dirichlet) + homogeneous
%         14. linear + lossless + source.ux (dirichlet) + heterogeneous
%
%         15. linear + lossy + source.p0 + homogeneous
%         16. linear + lossy + source.p0 + heterogeneous
%         17. linear + lossy + source.p (additive) + homogeneous
%         18. linear + lossy + source.p (additive) + heterogeneous
%         19. linear + lossy + source.p (additive-no-correction) + homogeneous
%         20. linear + lossy + source.p (additive-no-correction) + heterogeneous
%         21. linear + lossy + source.p (dirichlet) + homogeneous
%         22. linear + lossy + source.p (dirichlet) + heterogeneous
%         23. linear + lossy + source.ux (additive) + homogeneous
%         24. linear + lossy + source.ux (additive) + heterogeneous
%         25. linear + lossy + source.ux (additive-no-correction) + homogeneous
%         26. linear + lossy + source.ux (additive-no-correction) + heterogeneous
%         27. linear + lossy + source.ux (dirichlet) + homogeneous
%         28. linear + lossy + source.ux (dirichlet) + heterogeneous
%
%         29. linear + stokes + source.p0 + homogeneous
%         30. linear + stokes + source.p0 + heterogeneous
%         31. linear + stokes + source.p (additive) + homogeneous
%         32. linear + stokes + source.p (additive) + heterogeneous
%         33. linear + stokes + source.p (additive-no-correction) + homogeneous
%         34. linear + stokes + source.p (additive-no-correction) + heterogeneous
%         35. linear + stokes + source.p (dirichlet) + homogeneous
%         36. linear + stokes + source.p (dirichlet) + heterogeneous
%         37. linear + stokes + source.ux (additive) + homogeneous
%         38. linear + stokes + source.ux (additive) + heterogeneous
%         39. linear + stokes + source.ux (additive-no-correction) + homogeneous
%         40. linear + stokes + source.ux (additive-no-correction) + heterogeneous
%         41. linear + stokes + source.ux (dirichlet) + homogeneous
%         42. linear + stokes + source.ux (dirichlet) + heterogeneous
%
%         43. nonlinear + lossless + source.p0 + homogeneous
%         44. nonlinear + lossless + source.p0 + heterogeneous
%         45. nonlinear + lossless + source.p (additive) + homogeneous
%         46. nonlinear + lossless + source.p (additive) + heterogeneous
%         47. nonlinear + lossless + source.p (additive-no-correction) + homogeneous
%         48. nonlinear + lossless + source.p (additive-no-correction) + heterogeneous
%         49. nonlinear + lossless + source.p (dirichlet) + homogeneous
%         50. nonlinear + lossless + source.p (dirichlet) + heterogeneous
%         51. nonlinear + lossless + source.ux (additive) + homogeneous
%         52. nonlinear + lossless + source.ux (additive) + heterogeneous
%         53. nonlinear + lossless + source.ux (additive-no-correction) + homogeneous
%         54. nonlinear + lossless + source.ux (additive-no-correction) + heterogeneous
%         55. nonlinear + lossless + source.ux (dirichlet) + homogeneous
%         56. nonlinear + lossless + source.ux (dirichlet) + heterogeneous
%
%         57. nonlinear + lossy + source.p0 + homogeneous
%         58. nonlinear + lossy + source.p0 + heterogeneous
%         59. nonlinear + lossy + source.p (additive) + homogeneous
%         60. nonlinear + lossy + source.p (additive) + heterogeneous
%         61. nonlinear + lossy + source.p (additive-no-correction) + homogeneous
%         62. nonlinear + lossy + source.p (additive-no-correction) + heterogeneous
%         63. nonlinear + lossy + source.p (dirichlet) + homogeneous
%         64. nonlinear + lossy + source.p (dirichlet) + heterogeneous
%         65. nonlinear + lossy + source.ux (additive) + homogeneous
%         66. nonlinear + lossy + source.ux (additive) + heterogeneous
%         67. nonlinear + lossy + source.ux (additive-no-correction) + homogeneous
%         68. nonlinear + lossy + source.ux (additive-no-correction) + heterogeneous
%         69. nonlinear + lossy + source.ux (dirichlet) + homogeneous
%         70. nonlinear + lossy + source.ux (dirichlet) + heterogeneous
%
%         71. nonlinear + stokes + source.p0 + homogeneous
%         72. nonlinear + stokes + source.p0 + heterogeneous
%         73. nonlinear + stokes + source.p (additive) + homogeneous
%         74. nonlinear + stokes + source.p (additive) + heterogeneous
%         75. nonlinear + stokes + source.p (additive-no-correction) + homogeneous
%         76. nonlinear + stokes + source.p (additive-no-correction) + heterogeneous
%         77. nonlinear + stokes + source.p (dirichlet) + homogeneous
%         78. nonlinear + stokes + source.p (dirichlet) + heterogeneous
%         79. nonlinear + stokes + source.ux (additive) + homogeneous
%         80. nonlinear + stokes + source.ux (additive) + heterogeneous
%         81. nonlinear + stokes + source.ux (additive-no-correction) + homogeneous
%         82. nonlinear + stokes + source.ux (additive-no-correction) + heterogeneous
%         83. nonlinear + stokes + source.ux (dirichlet) + homogeneous
%         84. nonlinear + stokes + source.ux (dirichlet) + heterogeneous
%
%         85. linear + lossless + source.uy (additive) + homogeneous
%         86. linear + lossless + source.uy (additive) + heterogeneous
%         87. linear + lossless + source.uy (additive-no-correction) + homogeneous
%         88. linear + lossless + source.uy (additive-no-correction) + heterogeneous
%         89. linear + lossless + source.uy (dirichlet) + homogeneous
%         90. linear + lossless + source.uy (dirichlet) + heterogeneous
%
% ABOUT:
%     author      - Bradley Treeby (original), Claude Code (2D adaptation)
%     date        - 27th February 2026
%     last update - 27th February 2026
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2026- Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
%
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>.

%#ok<*NOPRT>

% check for plot inputs, and set to false if nargin is zero
if nargin == 0
    plot_comparisons = false;
    plot_simulations = false;
end

% =========================================================================
% DETECT BACKEND AND MODE
% =========================================================================

% Detect if Python backend is active by checking where kspaceFirstOrder2D resolves
% Note: Shim detection requires correct MATLAB path order. When using shims:
%   addpath('k-Wave');      % Add k-Wave first
%   addpath('tests/shims'); % Add shims last (goes to front of search path)
func_path = which('kspaceFirstOrder2D');
using_python_backend = contains(func_path, 'shims');

% Debug: Show which function is being used
disp(['DEBUG: kspaceFirstOrder2D resolves to: ' func_path]);
disp(['DEBUG: Using Python backend: ' mat2str(using_python_backend)]);

% Reference data file path
ref_data_file = fullfile(fileparts(mfilename('fullpath')), 'kspaceFirstOrder2D_reference_data.mat');

disp(['DEBUG: Reference data file: ' ref_data_file]);
disp(['DEBUG: Reference file exists: ' mat2str(exist(ref_data_file, 'file') == 2)]);

if using_python_backend
    % Load reference data for comparison
    if ~exist(ref_data_file, 'file')
        error('Reference data not found. Run test without Python shim first to generate reference data.');
    end
    disp('Python backend detected. Loading MATLAB reference data for comparison...');
    ref_data = load(ref_data_file);
    disp(['Reference data loaded: ' num2str(length(ref_data.reference_results)) ' test cases']);
else
    % Will generate reference data
    disp('MATLAB backend detected. Will generate reference data...');
    reference_results = cell(90, 1);
end

% set additional literals
STAGGERED_GRID      = true;
USE_KSPACE          = true;
USE_PML             = true;
SMOOTH_P0_SOURCE    = false;

% =========================================================================
% SIMULATION PARAMETERS
% =========================================================================

% define grid size
Nx = 32;
dx = 1;
Ny = 32;
dy = 1;

% define PML properties
PML_size    = 10;
if USE_PML
    PML_alpha_def   = 2;
else
    PML_alpha_def   = 0;
end

% define medium properties
c0      = 1500;
rho0    = 1000;
BonA0   = 10;
alpha0  = 5;
y       = 1.5;

% define medium properties for heterogeneous case
c1      = 2000;
rho1    = 1200;
BonA1   = 5;
alpha1  = 2;

% position of the heterogeneous interface
interface_position = Nx / 2;

% define time array
cfl = 0.1;
Nt  = 50;
dt  = cfl * dx / c0;
t_array = 0:dt:(Nt - 1) * dt;

% define optional inputs
input_args = {'PMLSize', PML_size, 'Smooth', false, ...
    'UsekSpace', USE_KSPACE, 'UseSG', STAGGERED_GRID, ...
    'PlotSim', plot_simulations};

% define source properties
source_strength = 5e6;
source_position_x = PML_size - 6;
source_position_y = Ny / 2;
source_freq     = 200;
source_signal   = source_strength * sin(2 * pi * source_freq * t_array);

% define sensor properties
sensor_position_x = Nx - PML_size;
sensor_position_y = Ny / 2;

% set pass variable
test_pass = true;

% test names
test_names = {...
    'linear + lossless + source.p0 + homogeneous', ...
    'linear + lossless + source.p0 + heterogeneous', ...
    'linear + lossless + source.p (additive) + homogeneous', ...
    'linear + lossless + source.p (additive) + heterogeneous', ...
    'linear + lossless + source.p (additive-no-correction) + homogeneous', ...
    'linear + lossless + source.p (additive-no-correction) + heterogeneous', ...
    'linear + lossless + source.p (dirichlet) + homogeneous', ...
    'linear + lossless + source.p (dirichlet) + heterogeneous', ...
    'linear + lossless + source.ux (additive) + homogeneous', ...
    'linear + lossless + source.ux (additive) + heterogeneous', ...
    'linear + lossless + source.ux (additive-no-correction) + homogeneous', ...
    'linear + lossless + source.ux (additive-no-correction) + heterogeneous', ...
    'linear + lossless + source.ux (dirichlet) + homogeneous', ...
    'linear + lossless + source.ux (dirichlet) + heterogeneous', ...
    'linear + lossy + source.p0 + homogeneous', ...
    'linear + lossy + source.p0 + heterogeneous', ...
    'linear + lossy + source.p (additive) + homogeneous', ...
    'linear + lossy + source.p (additive) + heterogeneous', ...
    'linear + lossy + source.p (additive-no-correction) + homogeneous', ...
    'linear + lossy + source.p (additive-no-correction) + heterogeneous', ...
    'linear + lossy + source.p (dirichlet) + homogeneous', ...
    'linear + lossy + source.p (dirichlet) + heterogeneous', ...
    'linear + lossy + source.ux (additive) + homogeneous', ...
    'linear + lossy + source.ux (additive) + heterogeneous', ...
    'linear + lossy + source.ux (additive-no-correction) + homogeneous', ...
    'linear + lossy + source.ux (additive-no-correction) + heterogeneous', ...
    'linear + lossy + source.ux (dirichlet) + homogeneous', ...
    'linear + lossy + source.ux (dirichlet) + heterogeneous', ...
    'linear + stokes + source.p0 + homogeneous', ...
    'linear + stokes + source.p0 + heterogeneous', ...
    'linear + stokes + source.p (additive) + homogeneous', ...
    'linear + stokes + source.p (additive) + heterogeneous', ...
    'linear + stokes + source.p (additive-no-correction) + homogeneous', ...
    'linear + stokes + source.p (additive-no-correction) + heterogeneous', ...
    'linear + stokes + source.p (dirichlet) + homogeneous', ...
    'linear + stokes + source.p (dirichlet) + heterogeneous', ...
    'linear + stokes + source.ux (additive) + homogeneous', ...
    'linear + stokes + source.ux (additive) + heterogeneous', ...
    'linear + stokes + source.ux (additive-no-correction) + homogeneous', ...
    'linear + stokes + source.ux (additive-no-correction) + heterogeneous', ...
    'linear + stokes + source.ux (dirichlet) + homogeneous', ...
    'linear + stokes + source.ux (dirichlet) + heterogeneous', ...
    'nonlinear + lossless + source.p0 + homogeneous', ...
    'nonlinear + lossless + source.p0 + heterogeneous', ...
    'nonlinear + lossless + source.p (additive) + homogeneous', ...
    'nonlinear + lossless + source.p (additive) + heterogeneous', ...
    'nonlinear + lossless + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + lossless + source.p (additive-no-correction) + heterogeneous', ...
    'nonlinear + lossless + source.p (dirichlet) + homogeneous', ...
    'nonlinear + lossless + source.p (dirichlet) + heterogeneous', ...
    'nonlinear + lossless + source.ux (additive) + homogeneous', ...
    'nonlinear + lossless + source.ux (additive) + heterogeneous', ...
    'nonlinear + lossless + source.ux (additive-no-correction) + homogeneous', ...
    'nonlinear + lossless + source.ux (additive-no-correction) + heterogeneous', ...
    'nonlinear + lossless + source.ux (dirichlet) + homogeneous', ...
    'nonlinear + lossless + source.ux (dirichlet) + heterogeneous', ...
    'nonlinear + lossy + source.p0 + homogeneous', ...
    'nonlinear + lossy + source.p0 + heterogeneous', ...
    'nonlinear + lossy + source.p (additive) + homogeneous', ...
    'nonlinear + lossy + source.p (additive) + heterogeneous', ...
    'nonlinear + lossy + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + lossy + source.p (additive-no-correction) + heterogeneous', ...
    'nonlinear + lossy + source.p (dirichlet) + homogeneous', ...
    'nonlinear + lossy + source.p (dirichlet) + heterogeneous', ...
    'nonlinear + lossy + source.ux (additive) + homogeneous', ...
    'nonlinear + lossy + source.ux (additive) + heterogeneous', ...
    'nonlinear + lossy + source.ux (additive-no-correction) + homogeneous', ...
    'nonlinear + lossy + source.ux (additive-no-correction) + heterogeneous', ...
    'nonlinear + lossy + source.ux (dirichlet) + homogeneous', ...
    'nonlinear + lossy + source.ux (dirichlet) + heterogeneous', ...
    'nonlinear + stokes + source.p0 + homogeneous', ...
    'nonlinear + stokes + source.p0 + heterogeneous', ...
    'nonlinear + stokes + source.p (additive) + homogeneous', ...
    'nonlinear + stokes + source.p (additive) + heterogeneous', ...
    'nonlinear + stokes + source.p (additive-no-correction) + homogeneous', ...
    'nonlinear + stokes + source.p (additive-no-correction) + heterogeneous', ...
    'nonlinear + stokes + source.p (dirichlet) + homogeneous', ...
    'nonlinear + stokes + source.p (dirichlet) + heterogeneous', ...
    'nonlinear + stokes + source.ux (additive) + homogeneous', ...
    'nonlinear + stokes + source.ux (additive) + heterogeneous', ...
    'nonlinear + stokes + source.ux (additive-no-correction) + homogeneous', ...
    'nonlinear + stokes + source.ux (additive-no-correction) + heterogeneous', ...
    'nonlinear + stokes + source.ux (dirichlet) + homogeneous', ...
    'nonlinear + stokes + source.ux (dirichlet) + heterogeneous', ...
    'linear + lossless + source.uy (additive) + homogeneous', ...
    'linear + lossless + source.uy (additive) + heterogeneous', ...
    'linear + lossless + source.uy (additive-no-correction) + homogeneous', ...
    'linear + lossless + source.uy (additive-no-correction) + heterogeneous', ...
    'linear + lossless + source.uy (dirichlet) + homogeneous', ...
    'linear + lossless + source.uy (dirichlet) + heterogeneous', ...
    };

% lists used to set properties
p0_tests = [1, 2, 15, 16, 29, 30, 43, 44, 57, 58, 71, 72];
p_tests  = [3:8, 17:22, 31:36, 45:50, 59:64, 73:78];
u_tests  = [9:14, 23:28, 37:42, 51:56, 65:70, 79:84];
additive_tests = [3,  4,  9,  10, 17, 18, 23, 24, 31, 32, 37, 38, ...
                  45, 46, 51, 52, 59, 60, 65, 66, 73, 74, 79, 80];
additive_no_correction_tests = additive_tests + 2;
dirichlet_tests = additive_tests + 4;

uy_tests = 85:90;  % source.uy scenarios (2D only)
% uy scenarios follow same additive/no-correction/dirichlet pattern
uy_additive_tests = [85, 86];
uy_additive_no_correction_tests = [87, 88];
uy_dirichlet_tests = [89, 90];

% =========================================================================
% SIMULATIONS
% =========================================================================

% loop through tests
for test_num = 1:90

    % clear structures
    clear source medium

    % update command line
    disp(['Running Test ' num2str(test_num) ' of 90: ' test_names{test_num}]);

    % assign medium properties
    medium.sound_speed  = c0;
    medium.density      = rho0;
    PML_alpha = PML_alpha_def;
    switch test_num
        case {1,2,3,4,5,6,7,8,9,10,11,12,13,14}

            % linear + lossless

        case {15,16,17,18,19,20,21,22,23,24,25,26,27,28}

            % linear + lossy
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = y;

        case {29,30,31,32,33,34,35,36,37,38,39,40,41,42}

            % linear + stokes
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = 2;
            medium.alpha_mode   = 'stokes';

        case {43,44,45,46,47,48,49,50,51,52,53,54,55,56}

            % nonlinear + lossless
            medium.BonA         = BonA0;

        case {57,58,59,60,61,62,63,64,65,66,67,68,69,70}

            % nonlinear + lossy
            medium.BonA         = BonA0;
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = y;

        case {71,72,73,74,75,76,77,78,79,80,81,82,83,84}

            % nonlinear + stokes
            medium.BonA         = BonA0;
            medium.alpha_coeff  = alpha0;
            medium.alpha_power  = 2;
            medium.alpha_mode   = 'stokes';

        case {85,86,87,88,89,90}

            % linear + lossless (source.uy scenarios)

        otherwise

            error('Unknown test number.');

    end

    % ----------------
    % 2D SIMULATION
    % ----------------

    % create computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    kgrid.t_array = t_array;

    % heterogeneous medium properties
    if ~rem(test_num, 2)
        setMaterialProperties();
    end

    % source
    clear source;
    if any(p0_tests == test_num)
        source.p0 = zeros(Nx, Ny);
        source.p0(source_position_x, source_position_y) = source_strength;
        if SMOOTH_P0_SOURCE
             source.p0 = smooth(source.p0, true);
        end
    elseif any(p_tests == test_num)
        source.p_mask = zeros(Nx, Ny);
        source.p_mask(source_position_x, source_position_y) = 1;
        source.p = source_signal;
        if any(dirichlet_tests == test_num)
            source.p_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.p_mode = 'additive';
            PML_alpha = 0;
        elseif any(additive_no_correction_tests == test_num)
            source.p_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end
    elseif any(u_tests == test_num)
        source.u_mask = zeros(Nx, Ny);
        source.u_mask(source_position_x, source_position_y) = 1;
        source.ux = source_signal ./ (c0 * rho0);
        if any(dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end
    elseif any(uy_tests == test_num)
        source.u_mask = zeros(Nx, Ny);
        source.u_mask(source_position_x, source_position_y) = 1;
        source.uy = source_signal ./ (c0 * rho0);
        if any(uy_dirichlet_tests == test_num)
            source.u_mode = 'dirichlet';
        elseif any(uy_additive_tests == test_num)
            source.u_mode = 'additive';
        elseif any(uy_additive_no_correction_tests == test_num)
            source.u_mode = 'additive-no-correction';
        else
            error('Unknown source mode.');
        end
    else
        error('Unknown source condition.');
    end

    % sensor
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(sensor_position_x, sensor_position_y) = 1;

    % run simulation
    sensor_data_2D = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:}, 'PMLAlpha', PML_alpha);

    % -------------
    % VALIDATION
    % -------------

    if using_python_backend
        % Compare against MATLAB reference
        ref_sensor_data = ref_data.reference_results{test_num};

        % Compute relative error
        ref_max = max(abs(ref_sensor_data(:)));
        if ref_max > 0
            rel_error = max(abs(sensor_data_2D(:) - ref_sensor_data(:))) / ref_max;
        else
            rel_error = max(abs(sensor_data_2D(:)));
        end

        % Check threshold
        COMPARISON_THRESH = 1e-14;
        if rel_error > COMPARISON_THRESH
            disp(['  TEST FAILED: Relative error = ' num2str(rel_error, '%.2e') ' > ' num2str(COMPARISON_THRESH, '%.2e')]);
            test_pass = false;
        else
            disp(['  Passed (rel_error = ' num2str(rel_error, '%.2e') ')']);
        end
    else
        % Store reference data for MATLAB baseline
        reference_results{test_num} = sensor_data_2D;

        % Basic sanity check
        if isempty(sensor_data_2D) || any(isnan(sensor_data_2D(:))) || any(isinf(sensor_data_2D(:)))
            disp(['  TEST FAILED: Invalid sensor data']);
            test_pass = false;
        end
    end

    % -------------
    % PLOTTING
    % -------------

    if plot_comparisons
        figure;
        plot(sensor_data_2D);
        title(test_names{test_num});
        xlabel('Time Step');
        ylabel('Pressure');
        drawnow;
    end

end

% =========================================================================
% SAVE REFERENCE DATA (MATLAB BASELINE MODE)
% =========================================================================

if ~using_python_backend
    % Save reference data for future Python backend comparisons
    save(ref_data_file, 'reference_results', '-v7.3');
    disp(['Reference data saved to: ' ref_data_file]);
    disp(['File size: ' num2str(dir(ref_data_file).bytes / 1024 / 1024, '%.2f') ' MB']);
end

% =========================================================================
% SUB FUNCTIONS
% =========================================================================

function setMaterialProperties()

    % sound speed and density (heterogeneous step in x-direction, uniform in y)
    medium.sound_speed = c0 * ones(Nx, Ny);
    medium.density = rho0 * ones(Nx, Ny);
    medium.sound_speed(interface_position:end, :) = c1;
    medium.density(interface_position:end, :) = rho1;

    % absorption (if defined)
    if isfield(medium, 'alpha_coeff')
        medium.alpha_coeff = alpha0 * ones(Nx, Ny);
        medium.alpha_coeff(interface_position:end, :) = alpha1;
    end

    % nonlinearity (if defined)
    if isfield(medium, 'BonA')
        medium.BonA = BonA0 * ones(Nx, Ny);
        medium.BonA(interface_position:end, :) = BonA1;
    end

end

end
