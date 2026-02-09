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
        _build_physics_ops(medium, source, rho0, xp)

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

def _build_physics_ops(medium, source, rho0, xp):
    """
    Build composable physics operators.
    Each returns 0 (or identity) if feature not present.
    """
    # --- Absorption (power-law attenuation) ---
    # Not implemented yet - returns 0
    absorption = lambda duxdx: 0
    dispersion = lambda rho: 0

    # --- Nonlinearity (BonA parameter) ---
    # Not implemented yet - returns 0 for p term, 1.0 for mass conservation
    BonA = medium.get("BonA", 0)
    if BonA != 0:
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

def interop_sanity(arr):
    """Verify MATLAB/Python data layout."""
    a = np.array(arr, copy=True)
    a[0, 1] = 99
    return a

def _cpu(x): return x.get() if hasattr(x, "get") else x
