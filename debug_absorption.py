#!/usr/bin/env python3
"""Debug absorption operator by comparing intermediate values with MATLAB."""
import sys
sys.path.insert(0, 'k-Wave/python')
import kWavePy
import numpy as np

# Test parameters matching MATLAB test
Nx = 64
dx = 1.0
dt = 0.1 * dx / 2000  # CFL with max c
Nt = 100

# Medium parameters (linear + lossy + homogeneous)
c0 = 1500.0
rho0 = 1000.0
alpha_coeff = 5.0  # dB/(MHz^y cm)
alpha_power = 1.5

# Initial pressure
p0 = np.zeros(Nx)
p0[10] = 5e6

# Build k-space vector
k = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)

print("=== Debugging Absorption Operators ===\n")

# Step 1: dB to Nepers conversion
alpha_np = 100 * alpha_coeff * (1e-6 / (2 * np.pi))**alpha_power / (20 * np.log10(np.e))
print(f"1. alpha_np (Nepers): {alpha_np:.10e}")

# Step 2: Absorption coefficients
absorb_tau = -2 * alpha_np * c0**(alpha_power - 1)
absorb_eta = 2 * alpha_np * c0**alpha_power * np.tan(np.pi * alpha_power / 2)
print(f"2. absorb_tau: {absorb_tau:.10e}")
print(f"3. absorb_eta: {absorb_eta:.10e}")

# Step 3: Fractional Laplacian operators (matching kWavePy.py fix)
k_shifted = np.fft.fftshift(k)
k_mag = np.abs(k_shifted)
absorb_nabla1 = k_mag**(alpha_power - 2)
absorb_nabla1[np.isinf(absorb_nabla1)] = 0
absorb_nabla1_shifted = np.fft.ifftshift(absorb_nabla1)

print(f"\n3a. k_shifted first 5: {k_shifted[:5]}")
print(f"3b. k_mag first 5: {k_mag[:5]}")

absorb_nabla2 = k_mag**(alpha_power - 1)
absorb_nabla2[np.isinf(absorb_nabla2)] = 0
absorb_nabla2_shifted = np.fft.ifftshift(absorb_nabla2)

print(f"\n4. absorb_nabla1 (before ifftshift):")
print(f"   Shape: {absorb_nabla1.shape}")
print(f"   First 5: {absorb_nabla1[:5]}")
print(f"   Last 5: {absorb_nabla1[-5:]}")
print(f"   Max: {np.max(absorb_nabla1):.6e}")

print(f"\n5. absorb_nabla1 (after ifftshift):")
print(f"   First 5: {absorb_nabla1_shifted[:5]}")
print(f"   Last 5: {absorb_nabla1_shifted[-5:]}")

print(f"\n6. absorb_nabla2 (before ifftshift):")
print(f"   First 5: {absorb_nabla2[:5]}")
print(f"   Max: {np.max(absorb_nabla2):.6e}")

print(f"\n7. absorb_nabla2 (after ifftshift):")
print(f"   First 5: {absorb_nabla2_shifted[:5]}")

# Step 4: Test operators
print("\n=== Testing Operators ===")
test_field = np.sin(2 * np.pi * np.arange(Nx) / Nx)

def diff_k(f, nabla_op):
    return np.real(np.fft.ifft(nabla_op * np.fft.fft(f)))

result_nabla1 = diff_k(test_field, absorb_nabla1_shifted)
result_nabla2 = diff_k(test_field, absorb_nabla2_shifted)

print(f"8. Test field max: {np.max(np.abs(test_field)):.6e}")
print(f"9. Result nabla1 max: {np.max(np.abs(result_nabla1)):.6e}")
print(f"10. Result nabla2 max: {np.max(np.abs(result_nabla2)):.6e}")

# Step 5: Run full simulation
print("\n=== Running Full Simulation ===")
kgrid = {"Nx": Nx, "dx": dx, "dt": dt, "Nt": Nt}
medium = {"c0": c0, "rho0": rho0, "alpha_coeff": alpha_coeff, "alpha_power": alpha_power}
source = {"p0": p0}
sensor = {"mask": np.ones(Nx, dtype=bool)}

result = kWavePy.simulate(kgrid, medium, source, sensor)
print(f"11. Max sensor data: {np.max(np.abs(result['sensor_data'])):.6e}")
print(f"12. Sensor data shape: {result['sensor_data'].shape}")
print(f"13. Final pressure max: {np.max(np.abs(result['pressure'])):.6e}")

print("\n=== Compare with MATLAB reference ===")
print("Run this in MATLAB to get reference values:")
print("""
Nx = 64; dx = 1.0; dt = 0.1 * dx / 2000; Nt = 100;
kgrid = kWaveGrid(Nx, dx);
kgrid.setTime(Nt, dt);
medium.sound_speed = 1500;
medium.density = 1000;
medium.alpha_coeff = 5.0;
medium.alpha_power = 1.5;
p0 = zeros(Nx, 1); p0(11) = 5e6;
source.p0 = p0;
sensor.mask = ones(Nx, 1);
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor);
fprintf('MATLAB max sensor data: %.6e\\n', max(abs(sensor_data(:))));
""")
