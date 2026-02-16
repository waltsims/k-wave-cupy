"""
Minimal k-Wave Python Backend.
"""
import numpy as np
try: import cupy as cp
except ImportError: cp = None

def simulate(kgrid, medium, source, sensor, backend="auto"):
    """
    1D k-Space Pseudospectral Wave Propagator.
    """
    # 1. Setup Backend (CPU/GPU)
    xp = cp if cp and backend in ("auto", "gpu") else np
    
    # 2. Extract Simulation Parameters
    Nx, dx = int(kgrid["Nx"]), float(kgrid["dx"])
    Nt, dt = int(kgrid["Nt"]), float(kgrid["dt"])

    # 3. Load Physics (supports scalar expansion)
    def load(obj, keys, default=None, dtype=float):
        val = next((obj.get(k) for k in keys if obj.get(k) is not None), default)
        if val is None: raise ValueError(f"Missing required parameter: {keys[0]}")
        arr = xp.array(val, dtype=dtype).flatten(order="F")
        return xp.full(Nx, arr[0], dtype=dtype) if arr.size == 1 else arr

    c0   = load(medium, ["sound_speed", "c0"])
    rho0 = load(medium, ["density", "rho0"], 1000.0)
    p    = load(source, ["p0"], 0.0)
    mask = load(sensor, ["mask"], True, dtype=bool)

    # Interpolate density to staggered grid (x + dx/2) for heterogeneous media
    # Velocity is defined at x + dx/2, so we need density there for momentum equation
    if rho0.size > 1:
        # Staggered grid: rho_sgx[i] = 0.5 * (rho0[i] + rho0[i+1])
        # Last point uses extrapolation (same as MATLAB's NaN handling)
        rho0_sgx = xp.zeros_like(rho0)
        rho0_sgx[:-1] = 0.5 * (rho0[:-1] + rho0[1:])
        rho0_sgx[-1] = rho0[-1]  # Boundary: use same value
    else:
        rho0_sgx = rho0

    # CFL Check
    cfl = float(xp.max(c0) * dt / dx)
    if cfl > 1.0: print(f"Warning: Unstable CFL condition: {cfl:.2f} > 1.0")

    # 4. Initialize State
    u = xp.zeros_like(p)
    rho = p / c0**2
    sensor_data = xp.zeros((int(xp.sum(mask)), Nt), dtype=p.dtype)

    # 5. Precompute k-space Operators
    k = 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    kappa = xp.sinc((xp.max(c0) * k * dt / 2) / np.pi)
    op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)
    op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)

    def diff(f, op): return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # 6. Build Physics Operators (composable, return 0 if inactive)
    absorption, dispersion, nonlinearity, nonlinear_factor, source_p, source_u = \
        _build_physics_ops(medium, source, k, rho0, xp)

    # 7. Time Loop (Leapfrog with Staggered Grid)
    # Initialize u at t = -dt/2 (backward half-step)
    u += (dt / (2 * rho0_sgx)) * diff(p, op_grad)

    for t in range(Nt):
        sensor_data[:, t] = p[mask]

        # Momentum equation with velocity sources
        grad_p = diff(p, op_grad)
        u -= (dt / rho0_sgx) * (grad_p + source_u(t, u))

        # Mass conservation with nonlinearity
        duxdx = diff(u, op_div)
        rho -= (dt * rho0) * duxdx * nonlinear_factor(rho)

        # Equation of state (all physics terms compose additively)
        p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))

        # Pressure sources (additive or dirichlet)
        p = source_p(t, p)

    return {"sensor_data": _cpu(sensor_data), "pressure": _cpu(p)}

def _build_physics_ops(medium, source, k, rho0, xp):
    """
    Build composable physics operators.
    Each returns 0 (or identity) if feature not present.
    """
    # --- Absorption (power-law attenuation) ---
    absorption, dispersion = _build_absorption_ops(medium, k, rho0, xp)

    # --- Nonlinearity (BonA parameter) ---
    BonA = medium.get("BonA", 0)
    is_nonlinear = not (np.all(BonA == 0) if hasattr(BonA, '__len__') else BonA == 0)
    if is_nonlinear:
        nonlinearity = lambda rho: BonA * rho**2 / (2 * rho0)
        nonlinear_factor = lambda rho: (2*rho + rho0) / rho0
    else:
        nonlinearity = lambda rho: 0
        nonlinear_factor = lambda rho: 1.0

    # --- Pressure Sources (time-varying) ---
    # Not implemented yet - returns identity
    source_p = lambda t, p: p

    # --- Velocity Sources (time-varying) ---
    # Not implemented yet - returns 0
    source_u = lambda t, u: 0

    return absorption, dispersion, nonlinearity, nonlinear_factor, source_p, source_u

def _build_absorption_ops(medium, k, rho0, xp):
    """Build absorption and dispersion operators for power-law attenuation."""
    alpha_coeff = medium.get("alpha_coeff", 0)
    alpha_power = medium.get("alpha_power", 1.5)

    # Handle scalar or array alpha_coeff
    is_zero = np.all(alpha_coeff == 0) if hasattr(alpha_coeff, '__len__') else alpha_coeff == 0
    if is_zero:
        return lambda duxdx: 0, lambda rho: 0

    # Convert dB/(MHz^y cm) to Nepers/((rad/s)^y m)
    alpha_np = 100 * alpha_coeff * (1e-6 / (2 * np.pi))**alpha_power / (20 * np.log10(np.e))

    # Get sound speed (can be array for heterogeneous media)
    c0 = medium.get("sound_speed", medium.get("c0"))
    if not hasattr(c0, '__len__'):
        c0 = xp.array([c0])
    else:
        c0 = xp.array(c0)

    # Get reference sound speed for k-space operators (use max for heterogeneous)
    c_ref = float(xp.max(c0))

    # Stokes absorption (alpha_power = 2) - special case, no fractional Laplacian
    if abs(alpha_power - 2.0) < 1e-10:
        absorb_tau = -2 * alpha_np * c0
        return lambda duxdx: absorb_tau * rho0 * duxdx, lambda rho: 0

    # General power-law absorption - uses fractional Laplacian
    absorb_tau = -2 * alpha_np * c0**(alpha_power - 1)
    absorb_eta = 2 * alpha_np * c0**alpha_power * xp.tan(np.pi * alpha_power / 2)

    # Fractional Laplacian operators in k-space
    # MATLAB uses fftshift(k) for computing, so match that:
    # 1. Shift k to centered format
    # 2. Compute power
    # 3. Shift back to FFT ordering
    k_shifted = xp.fft.fftshift(k)
    k_mag = xp.abs(k_shifted)

    absorb_nabla1 = k_mag**(alpha_power - 2)
    absorb_nabla1 = xp.where(xp.isinf(absorb_nabla1), 0, absorb_nabla1)
    absorb_nabla1 = xp.fft.ifftshift(absorb_nabla1)

    absorb_nabla2 = k_mag**(alpha_power - 1)
    absorb_nabla2 = xp.where(xp.isinf(absorb_nabla2), 0, absorb_nabla2)
    absorb_nabla2 = xp.fft.ifftshift(absorb_nabla2)

    def diff_k(f, nabla_op):
        return xp.real(xp.fft.ifft(nabla_op * xp.fft.fft(f)))

    absorption = lambda duxdx: absorb_tau * diff_k(rho0 * duxdx, absorb_nabla1)
    dispersion = lambda rho: absorb_eta * diff_k(rho, absorb_nabla2)

    return absorption, dispersion

def interop_sanity(arr):
    """Verify MATLAB/Python data layout."""
    a = np.array(arr, copy=True)
    a[0, 1] = 99
    return a

def _cpu(x): return x.get() if hasattr(x, "get") else x
