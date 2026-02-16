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
    mask = load(sensor, ["mask"], True, dtype=bool)

    # Initialize pressure to zero (source.p0 is applied at t=0 like MATLAB)
    p = xp.zeros(Nx, dtype=float)

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
    c_ref = float(xp.max(c0))
    k = 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    kappa = xp.sinc((c_ref * k * dt / 2) / np.pi)
    source_kappa = xp.cos(c_ref * k * dt / 2)  # For additive source correction
    op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)
    op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)

    def diff(f, op): return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # 6. Build Physics Operators (composable, return 0 if inactive)
    absorption, dispersion, nonlinearity, nonlinear_factor = \
        _build_physics_ops(medium, source, k, rho0, xp)

    # 7. Build Source Operators
    source_p_op, source_u_op = _build_source_ops(source, c0, rho0, dt, dx, Nx, Nt, source_kappa, xp)

    # 8. Time Loop (Leapfrog with Staggered Grid)
    # Check if using initial pressure source (p0) vs time-varying (p)
    p0_raw = source.get("p0", 0)
    has_p0 = not (np.all(p0_raw == 0) if hasattr(p0_raw, '__len__') else p0_raw == 0)
    p0_initial = xp.array(p0_raw, dtype=float).flatten(order="F") if has_p0 else None

    # Initialize u at t = -dt/2 (backward half-step)
    # Note: For source.p0, p starts at 0, so this gives u=0 initially (like MATLAB)
    u += (dt / (2 * rho0_sgx)) * diff(p, op_grad)

    for t in range(Nt):
        # Momentum equation
        grad_p = diff(p, op_grad)
        u -= (dt / rho0_sgx) * grad_p

        # Velocity sources (applied after momentum equation)
        u = source_u_op(t, u)

        # Mass conservation with nonlinearity
        duxdx = diff(u, op_div)
        rho -= (dt * rho0) * duxdx * nonlinear_factor(rho)

        # Pressure sources (applied to rho before equation of state)
        rho = source_p_op(t, rho)

        # Equation of state (all physics terms compose additively)
        p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))

        # For source.p0: override at t=0 with initial pressure and velocity (MATLAB behavior)
        if t == 0 and has_p0:
            p = p0_initial.copy()
            rho = p / c0**2
            # Set initial velocity: u(t=-dt/2) such that u(t=0) = 0
            # This is u = (dt/2) * (1/rho0_sgx) * grad(p)
            u = (dt / (2 * rho0_sgx)) * diff(p, op_grad)

        # Record sensor data at END of time step (like MATLAB)
        sensor_data[:, t] = p[mask]

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

    return absorption, dispersion, nonlinearity, nonlinear_factor

def _build_source_ops(source, c0, rho0, dt, dx, Nx, Nt, source_kappa, xp):
    """Build time-varying source operators for pressure and velocity."""

    def apply_kspace_correction(source_mat):
        return xp.real(xp.fft.ifft(source_kappa * xp.fft.fft(source_mat)))

    # --- Pressure Source ---
    p_mask_raw = source.get("p_mask", 0)
    p_signal_raw = source.get("p", 0)
    p_mode = source.get("p_mode", "additive")

    has_p_mask = not (np.all(p_mask_raw == 0) if hasattr(p_mask_raw, '__len__') else p_mask_raw == 0)
    has_p_signal = not (np.all(p_signal_raw == 0) if hasattr(p_signal_raw, '__len__') else p_signal_raw == 0)

    if has_p_mask and has_p_signal:
        p_mask = xp.array(p_mask_raw, dtype=bool).flatten(order="F")
        p_signal = xp.array(p_signal_raw, dtype=float).flatten(order="F")

        # Scale source signal (MATLAB: kspaceFirstOrder_scaleSourceTerms.m)
        # N = 1 for 1D
        c0_at_source = c0[p_mask] if c0.size > 1 else c0
        if p_mode == "dirichlet":
            # Dirichlet: scale by 1/(N * c0^2) = 1/c0^2
            p_scaled = p_signal / (c0_at_source**2)
        else:
            # Additive: scale by 2*dt/(N * c0 * dx) = 2*dt/(c0 * dx)
            p_scaled = p_signal * (2 * dt / (c0_at_source * dx))

        if p_mode == "dirichlet":
            def source_p_op(t, rho):
                if t < len(p_scaled):
                    rho[p_mask] = p_scaled[t] if p_scaled.ndim == 1 else p_scaled[:, t]
                return rho
        elif p_mode == "additive":
            def source_p_op(t, rho):
                if t < len(p_scaled):
                    source_mat = xp.zeros(Nx, dtype=rho.dtype)
                    source_mat[p_mask] = p_scaled[t] if p_scaled.ndim == 1 else p_scaled[:, t]
                    rho = rho + apply_kspace_correction(source_mat)
                return rho
        else:  # additive-no-correction
            def source_p_op(t, rho):
                if t < len(p_scaled):
                    rho[p_mask] = rho[p_mask] + (p_scaled[t] if p_scaled.ndim == 1 else p_scaled[:, t])
                return rho
    else:
        source_p_op = lambda t, rho: rho

    # --- Velocity Source ---
    u_mask_raw = source.get("u_mask", 0)
    u_signal_raw = source.get("ux", 0)
    u_mode = source.get("u_mode", "additive")

    has_u_mask = not (np.all(u_mask_raw == 0) if hasattr(u_mask_raw, '__len__') else u_mask_raw == 0)
    has_u_signal = not (np.all(u_signal_raw == 0) if hasattr(u_signal_raw, '__len__') else u_signal_raw == 0)

    if has_u_mask and has_u_signal:
        u_mask = xp.array(u_mask_raw, dtype=bool).flatten(order="F")
        u_signal = xp.array(u_signal_raw, dtype=float).flatten(order="F")

        # Scale velocity source (MATLAB: kspaceFirstOrder_scaleSourceTerms.m)
        # Dirichlet: no scaling (values used directly)
        # Additive: scale by 2*c0*dt/dx
        c0_at_source = c0[u_mask] if c0.size > 1 else c0
        if u_mode == "dirichlet":
            u_scaled = u_signal  # No scaling for dirichlet
        else:
            u_scaled = u_signal * (2 * c0_at_source * dt / dx)

        if u_mode == "dirichlet":
            def source_u_op(t, u):
                if t < len(u_scaled):
                    u[u_mask] = u_scaled[t] if u_scaled.ndim == 1 else u_scaled[:, t]
                return u
        elif u_mode == "additive":
            def source_u_op(t, u):
                if t < len(u_scaled):
                    source_mat = xp.zeros(Nx, dtype=u.dtype)
                    source_mat[u_mask] = u_scaled[t] if u_scaled.ndim == 1 else u_scaled[:, t]
                    u = u + apply_kspace_correction(source_mat)
                return u
        else:  # additive-no-correction
            def source_u_op(t, u):
                if t < len(u_scaled):
                    u[u_mask] = u[u_mask] + (u_scaled[t] if u_scaled.ndim == 1 else u_scaled[:, t])
                return u
    else:
        source_u_op = lambda t, u: u

    return source_p_op, source_u_op


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
