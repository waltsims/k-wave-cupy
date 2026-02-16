"""
Minimal k-Wave Python Backend.
"""
import numpy as np
try: import cupy as cp
except ImportError: cp = None

def _is_nonzero(x):
    """Check if x is nonzero (handles scalars and arrays uniformly)."""
    return not (np.all(x == 0) if hasattr(x, '__len__') else x == 0)

def _to_cpu(x):
    """Move CuPy array to CPU if needed."""
    return x.get() if hasattr(x, "get") else x

def _expand_to_grid(val, Nx, xp, name="parameter"):
    """Ensure medium parameters are 1D arrays matching the grid length."""
    arr = xp.array(val, dtype=float).flatten(order="F")
    if arr.size == 1: return xp.full(Nx, float(arr[0]), dtype=float)
    if arr.size == Nx: return arr
    raise ValueError(f"{name} length {arr.size} incompatible with grid size {Nx}")

def simulate(kgrid, medium, source, sensor, backend="auto"):
    """1D k-Space Pseudospectral Wave Propagator."""
    # Use GPU if available and requested
    xp = cp if cp and backend in ("auto", "gpu") else np

    Nx, dx = int(kgrid["Nx"]), float(kgrid["dx"])
    Nt, dt = int(kgrid["Nt"]), float(kgrid["dt"])

    # Load with scalar expansion: single values become uniform arrays
    def load_field(obj, keys, default=None, dtype=float):
        val = next((obj.get(k) for k in keys if obj.get(k) is not None), default)
        if val is None: raise ValueError(f"Missing required parameter: {keys[0]}")
        arr = xp.array(val, dtype=dtype).flatten(order="F")
        return xp.full(Nx, arr[0], dtype=dtype) if arr.size == 1 else arr

    c0   = load_field(medium, ["sound_speed", "c0"])
    rho0 = load_field(medium, ["density", "rho0"], 1000.0)
    mask_raw = sensor.get("mask")
    if mask_raw is None:
        mask = xp.ones(Nx, dtype=bool)  # Record all points by default
    else:
        # Force copy to avoid memory aliasing with MATLAB-created arrays
        mask = xp.array(mask_raw, dtype=bool, copy=True).flatten(order="F")
        if mask.size == 1:
            mask = xp.full(Nx, bool(mask[0]), dtype=bool)

    # Validate mask shape
    if mask.shape != (Nx,):
        raise ValueError(f"Sensor mask shape {mask.shape} doesn't match grid size {Nx}")

    p = xp.zeros(Nx, dtype=float)

    # Density at staggered grid points (x + dx/2) for velocity update
    rho0_sgx = xp.concatenate([0.5 * (rho0[:-1] + rho0[1:]), rho0[-1:]]) if rho0.size > 1 else rho0

    # Warn if simulation will be unstable
    cfl = float(xp.max(c0) * dt / dx)
    if cfl > 1.0: print(f"Warning: Unstable CFL condition: {cfl:.2f} > 1.0")

    # k-space operators with time-staggering correction (sinc factor)
    c_ref, k = float(xp.max(c0)), 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    kappa, source_kappa = xp.sinc((c_ref * k * dt / 2) / np.pi), xp.cos(c_ref * k * dt / 2)
    op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)  # Gradient with forward shift
    op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)  # Divergence with backward shift
    spectral_diff = lambda f, op: xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # Physics operators return 0 or identity when feature is disabled
    absorption, dispersion, nonlinearity, nonlinear_factor = _build_physics_ops(medium, k, rho0, Nx, xp)
    source_p_op, source_u_op = _build_source_ops(source, c0, dt, dx, Nx, source_kappa, xp)

    # Initial pressure source (applied at t=0, overriding computed values)
    p0_raw = source.get("p0", 0)
    p0_initial = xp.array(p0_raw, dtype=float).flatten(order="F") if _is_nonzero(p0_raw) else None

    u, rho = xp.zeros_like(p), xp.zeros_like(p)
    sensor_data = xp.zeros((int(xp.sum(mask)), Nt), dtype=p.dtype)

    # Initialize velocity at t=-dt/2 for leapfrog staggering
    u += (dt / (2 * rho0_sgx)) * spectral_diff(p, op_grad)

    for t in range(Nt):
        # Momentum equation: du/dt = -grad(p)/rho
        u -= (dt / rho0_sgx) * spectral_diff(p, op_grad)
        u = source_u_op(t, u)

        # Mass conservation: drho/dt = -rho0 * div(u) * nonlinear_factor
        duxdx = spectral_diff(u, op_div)
        rho -= (dt * rho0) * duxdx * nonlinear_factor(rho)
        rho = source_p_op(t, rho)

        # Equation of state with absorption/dispersion/nonlinearity
        p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))

        # For source.p0: override computed values at t=0 (MATLAB convention)
        if t == 0 and p0_initial is not None:
            p, rho = p0_initial.copy(), p0_initial / c0**2
            u = (dt / (2 * rho0_sgx)) * spectral_diff(p, op_grad)

        sensor_data[:, t] = p[mask]

    return {"sensor_data": _to_cpu(sensor_data), "pressure": _to_cpu(p)}

def _build_physics_ops(medium, k, rho0, Nx, xp):
    """Build composable physics operators (return 0/identity when disabled)."""
    absorption, dispersion = _build_absorption_ops(medium, k, rho0, Nx, xp)

    BonA_raw = medium.get("BonA", 0)
    if not _is_nonzero(BonA_raw):
        return absorption, dispersion, lambda rho: 0, lambda rho: 1.0

    # Nonlinear acoustics: pressure depends on rho^2
    BonA = _expand_to_grid(BonA_raw, Nx, xp, "BonA")
    return (absorption, dispersion,
            lambda rho: BonA * rho**2 / (2 * rho0),
            lambda rho: (2*rho + rho0) / rho0)

def _build_source_ops(source, c0, dt, dx, Nx, source_kappa, xp):
    """Build time-varying source operators for pressure and velocity."""

    def build_source_operator(mask_raw, signal_raw, mode, dirichlet_scale, additive_scale):
        """Generic factory for source operators (pressure or velocity)."""
        if not (_is_nonzero(mask_raw) and _is_nonzero(signal_raw)):
            return lambda t, field: field

        mask = xp.array(mask_raw, dtype=bool).flatten(order="F")
        signal = xp.array(signal_raw, dtype=float).flatten(order="F")
        c0_src = c0[mask] if c0.size > 1 else c0
        scaled = signal / dirichlet_scale(c0_src) if mode == "dirichlet" else signal * additive_scale(c0_src)
        get_value = lambda t: scaled[t] if scaled.ndim == 1 else scaled[:, t]

        if mode == "dirichlet":
            return lambda t, field: (field.__setitem__(mask, get_value(t)), field)[1] if t < len(scaled) else field
        if mode == "additive":
            def apply_additive(t, field):
                if t >= len(scaled): return field
                src = xp.zeros(Nx, dtype=field.dtype)
                src[mask] = get_value(t)
                return field + xp.real(xp.fft.ifft(source_kappa * xp.fft.fft(src)))
            return apply_additive
        # additive-no-correction
        return lambda t, field: (field.__setitem__(mask, field[mask] + get_value(t)), field)[1] if t < len(scaled) else field

    # Pressure source: scale by 1/c0^2 (dirichlet) or 2*dt/(c0*dx) (additive)
    source_p_op = build_source_operator(
        source.get("p_mask", 0), source.get("p", 0), source.get("p_mode", "additive"),
        lambda c: c**2, lambda c: 2*dt/(c*dx))

    # Velocity source: no scaling (dirichlet) or 2*c0*dt/dx (additive)
    source_u_op = build_source_operator(
        source.get("u_mask", 0), source.get("ux", 0), source.get("u_mode", "additive"),
        lambda c: 1, lambda c: 2*c*dt/dx)

    return source_p_op, source_u_op

def _build_absorption_ops(medium, k, rho0, Nx, xp):
    """Build absorption and dispersion operators for power-law attenuation."""
    alpha_coeff_raw = medium.get("alpha_coeff", 0)
    if not _is_nonzero(alpha_coeff_raw):
        return lambda duxdx: 0, lambda rho: 0

    alpha_coeff = _expand_to_grid(alpha_coeff_raw, Nx, xp, "alpha_coeff")
    c0 = _expand_to_grid(medium.get("sound_speed", medium.get("c0")), Nx, xp, "sound_speed")
    alpha_power = float(xp.array(medium.get("alpha_power", 1.5), dtype=float).flatten(order="F")[0])

    # Convert from dB/(MHz^y cm) to Nepers/((rad/s)^y m)
    alpha_np = 100 * alpha_coeff * (1e-6 / (2 * np.pi))**alpha_power / (20 * np.log10(np.e))
    spectral_diff = lambda f, op: xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # Stokes absorption (y=2): simple viscous damping, no fractional Laplacian needed
    if abs(alpha_power - 2.0) < 1e-10:
        return lambda duxdx: -2 * alpha_np * c0 * rho0 * duxdx, lambda rho: 0

    # General power-law: requires fractional Laplacian for causality
    def fractional_laplacian(power):
        k_mag = xp.abs(xp.fft.fftshift(k))
        return xp.fft.ifftshift(xp.where(k_mag == 0, 0, k_mag**power))

    tau, eta = -2 * alpha_np * c0**(alpha_power - 1), 2 * alpha_np * c0**alpha_power * xp.tan(np.pi * alpha_power / 2)
    nabla1, nabla2 = fractional_laplacian(alpha_power - 2), fractional_laplacian(alpha_power - 1)

    return (lambda duxdx: tau * spectral_diff(rho0 * duxdx, nabla1),
            lambda rho: eta * spectral_diff(rho, nabla2))

def interop_sanity(arr):
    """Verify MATLAB/Python data layout (for testing)."""
    a = np.array(arr, copy=True)
    a[0, 1] = 99
    return a
