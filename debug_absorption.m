% Debug absorption operators by printing intermediate values
clearvars;

% Test parameters
Nx = 64;
dx = 1.0;
dt = 0.1 * dx / 2000;
Nt = 100;

% Medium parameters (linear + lossy + homogeneous)
c0 = 1500.0;
rho0 = 1000.0;
alpha_coeff = 5.0;  % dB/(MHz^y cm)
alpha_power = 1.5;

fprintf('=== Debugging Absorption Operators (MATLAB) ===\n\n');

% Step 1: dB to Nepers conversion
alpha_np = db2neper(alpha_coeff, alpha_power);
fprintf('1. alpha_np (Nepers): %.10e\n', alpha_np);

% Step 2: Absorption coefficients
absorb_tau = -2 * alpha_np * c0^(alpha_power - 1);
absorb_eta = 2 * alpha_np * c0^alpha_power * tan(pi * alpha_power / 2);
fprintf('2. absorb_tau: %.10e\n', absorb_tau);
fprintf('3. absorb_eta: %.10e\n', absorb_eta);

% Step 3: Create k-space grid
kgrid = kWaveGrid(Nx, dx);
kgrid.setTime(Nt, dt);

% Step 4: Fractional Laplacian operators
absorb_nabla1 = (kgrid.k).^(alpha_power - 2);
absorb_nabla1(isinf(absorb_nabla1)) = 0;
absorb_nabla1_shifted = ifftshift(absorb_nabla1);

absorb_nabla2 = (kgrid.k).^(alpha_power - 1);
absorb_nabla2(isinf(absorb_nabla2)) = 0;
absorb_nabla2_shifted = ifftshift(absorb_nabla2);

fprintf('\n4. absorb_nabla1 (before ifftshift):\n');
fprintf('   Shape: [%d, %d]\n', size(absorb_nabla1));
fprintf('   First 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla1(1:5));
fprintf('   Last 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla1(end-4:end));
fprintf('   Max: %.6e\n', max(absorb_nabla1));

fprintf('\n5. absorb_nabla1 (after ifftshift):\n');
fprintf('   First 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla1_shifted(1:5));
fprintf('   Last 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla1_shifted(end-4:end));

fprintf('\n6. absorb_nabla2 (before ifftshift):\n');
fprintf('   First 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla2(1:5));
fprintf('   Max: %.6e\n', max(absorb_nabla2));

fprintf('\n7. absorb_nabla2 (after ifftshift):\n');
fprintf('   First 5: [%.6e %.6e %.6e %.6e %.6e]\n', absorb_nabla2_shifted(1:5));

% Step 5: Test operators
fprintf('\n=== Testing Operators ===\n');
test_field = sin(2 * pi * (0:Nx-1)' / Nx);

diff_k = @(f, nabla_op) real(ifft(nabla_op .* fft(f)));

result_nabla1 = diff_k(test_field, absorb_nabla1_shifted);
result_nabla2 = diff_k(test_field, absorb_nabla2_shifted);

fprintf('8. Test field max: %.6e\n', max(abs(test_field)));
fprintf('9. Result nabla1 max: %.6e\n', max(abs(result_nabla1)));
fprintf('10. Result nabla2 max: %.6e\n', max(abs(result_nabla2)));

% Step 6: Run full simulation
fprintf('\n=== Running Full Simulation ===\n');
medium.sound_speed = c0;
medium.density = rho0;
medium.alpha_coeff = alpha_coeff;
medium.alpha_power = alpha_power;

p0 = zeros(Nx, 1);
p0(11) = 5e6;
source.p0 = p0;

sensor.mask = ones(Nx, 1);

sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor);
fprintf('11. Max sensor data: %.6e\n', max(abs(sensor_data(:))));
fprintf('12. Sensor data shape: [%d, %d]\n', size(sensor_data));
fprintf('13. Final pressure would need p_final to compute\n');

% Additional debug: print kgrid.k values
fprintf('\n=== k-space vector ===\n');
fprintf('kgrid.k first 5: [%.6e %.6e %.6e %.6e %.6e]\n', kgrid.k(1:5));
fprintf('kgrid.k last 5: [%.6e %.6e %.6e %.6e %.6e]\n', kgrid.k(end-4:end));
fprintf('kgrid.k shape: [%d, %d]\n', size(kgrid.k));
