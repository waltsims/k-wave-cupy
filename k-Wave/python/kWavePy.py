"""
Minimal k-Wave Python Backend.
"""
import numpy as np
from types import SimpleNamespace
try: import cupy as cp
except ImportError: cp = None

def _attr(obj, name, default=None):
    """Get attribute with default (works for objects and SimpleNamespace)."""
    return getattr(obj, name, default)

def _is_enabled(x):
    """Check if parameter is nonzero/active (handles scalars and arrays)."""
    return x is not None and not (np.all(x == 0) if hasattr(x, '__len__') else x == 0)

def _to_cpu(x):
    """Move CuPy array to CPU if needed."""
    return x.get() if hasattr(x, "get") else x

def _expand_to_grid(val, Nx, xp, name="parameter"):
    """Expand scalar or array to match grid length Nx."""
    if val is None: raise ValueError(f"Missing required parameter: {name}")
    arr = xp.array(val, dtype=float).flatten(order="F")
    if arr.size == 1: return xp.full(Nx, float(arr[0]), dtype=float)
    if arr.size == Nx: return arr
    raise ValueError(f"{name} length {arr.size} incompatible with grid size {Nx}")

def _spectral_diff(f, op, xp):
    """Apply spectral operator: F^-1[op * F[f]] (core k-space method)."""
    return xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

def _fractional_laplacian(k, power, xp):
    """Fractional Laplacian |k|^power in k-space (zero at DC to avoid singularity)."""
    k_mag = xp.abs(xp.fft.fftshift(k))
    return xp.fft.ifftshift(xp.where(k_mag == 0, 0, k_mag**power))

def _parse_inputs(kgrid, medium, sensor, xp):
    """Parse and validate simulation inputs, expanding scalars to grid arrays."""
    Nx, dx = int(kgrid.Nx), float(kgrid.dx)
    Nt, dt = int(kgrid.Nt), float(kgrid.dt)

    c0   = _expand_to_grid(medium.sound_speed, Nx, xp, "sound_speed")
    rho0 = _expand_to_grid(_attr(medium, 'density', 1000.0), Nx, xp, "density")

    mask_raw = _attr(sensor, 'mask', None)
    if mask_raw is None:
        mask = xp.ones(Nx, dtype=bool)
    else:
        mask = xp.array(mask_raw, dtype=bool, copy=True).flatten(order="F")
        if mask.size == 1:
            mask = xp.full(Nx, bool(mask[0]), dtype=bool)
    if mask.shape != (Nx,):
        raise ValueError(f"Sensor mask shape {mask.shape} doesn't match grid size {Nx}")

    return Nx, dx, Nt, dt, c0, rho0, mask

def simulate(kgrid, medium, source, sensor, backend="auto"):
    """1D k-Space Pseudospectral Wave Propagator.

    Accepts objects with attributes (kWaveGrid, kWaveMedium, etc.) or SimpleNamespace.
    For dict-based input (MATLAB interop), use simulate_from_dicts().
    """
    xp = cp if cp and backend in ("auto", "gpu") else np
    Nx, dx, Nt, dt, c0, rho0, mask = _parse_inputs(kgrid, medium, sensor, xp)
    diff = lambda f, op: _spectral_diff(f, op, xp)

    p = xp.zeros(Nx, dtype=float)

    # Density at staggered grid points (x + dx/2) for velocity update
    rho0_staggered = xp.concatenate([0.5 * (rho0[:-1] + rho0[1:]), rho0[-1:]]) if rho0.size > 1 else rho0

    # Warn if simulation will be unstable
    cfl = float(xp.max(c0) * dt / dx)
    if cfl > 1.0: print(f"Warning: Unstable CFL condition: {cfl:.2f} > 1.0")

    # k-space operators with time-staggering correction (sinc factor)
    c_ref, k = float(xp.max(c0)), 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    kappa, source_kappa = xp.sinc((c_ref * k * dt / 2) / np.pi), xp.cos(c_ref * k * dt / 2)
    op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)  # Gradient with forward shift
    op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)  # Divergence with backward shift

    # Physics operators return 0 or identity when feature is disabled
    absorption, dispersion, nonlinearity, nonlinear_factor = _build_physics_ops(medium, k, rho0, Nx, xp)
    source_p_op, source_u_op = _build_source_ops(source, c0, dt, dx, Nx, source_kappa, xp)

    # Initial pressure source (applied at t=0, overriding computed values)
    p0_raw = _attr(source, 'p0', 0)
    p0_initial = xp.array(p0_raw, dtype=float).flatten(order="F") if _is_enabled(p0_raw) else None

    u, rho = xp.zeros_like(p), xp.zeros_like(p)
    sensor_data = xp.zeros((int(xp.sum(mask)), Nt), dtype=p.dtype)

    # Initialize velocity at t=-dt/2 for leapfrog staggering
    u += (dt / (2 * rho0_staggered)) * diff(p, op_grad)

    for t in range(Nt):
        # Momentum equation: du/dt = -grad(p)/rho
        u -= (dt / rho0_staggered) * diff(p, op_grad)
        u = source_u_op(t, u)

        # Mass conservation: drho/dt = -rho0 * div(u) * nonlinear_factor
        duxdx = diff(u, op_div)
        rho -= (dt * rho0) * duxdx * nonlinear_factor(rho)
        rho = source_p_op(t, rho)

        # Equation of state with absorption/dispersion/nonlinearity
        p = c0**2 * (rho + absorption(duxdx) - dispersion(rho) + nonlinearity(rho))

        # For source.p0: override computed values at t=0 (MATLAB convention)
        if t == 0 and p0_initial is not None:
            p, rho = p0_initial.copy(), p0_initial / c0**2
            u = (dt / (2 * rho0_staggered)) * diff(p, op_grad)

        sensor_data[:, t] = p[mask]

    return {"sensor_data": _to_cpu(sensor_data), "pressure": _to_cpu(p)}

def _build_physics_ops(medium, k, rho0, Nx, xp):
    """Build composable physics operators (return 0/identity when disabled)."""
    absorption, dispersion = _build_absorption_ops(medium, k, rho0, Nx, xp)

    BonA_raw = _attr(medium, 'BonA', 0)
    if not _is_enabled(BonA_raw):
        # (absorption, dispersion, nonlinearity, nonlinear_factor)
        return absorption, dispersion, lambda rho: 0, lambda rho: 1.0

    # Nonlinear acoustics: pressure depends on rho^2
    BonA = _expand_to_grid(BonA_raw, Nx, xp, "BonA")
    # (absorption, dispersion, nonlinearity, nonlinear_factor)
    return (absorption, dispersion,
            lambda rho: BonA * rho**2 / (2 * rho0),
            lambda rho: (2*rho + rho0) / rho0)

def _build_source_ops(source, c0, dt, dx, Nx, source_kappa, xp):
    """Build time-varying source operators for pressure and velocity."""

    def build_op(mask_raw, signal_raw, mode, scale_dirichlet, scale_additive):
        if not (_is_enabled(mask_raw) and _is_enabled(signal_raw)):
            return lambda t, field: field

        mask = xp.array(mask_raw, dtype=bool).flatten(order="F")
        signal = xp.array(signal_raw, dtype=float).flatten(order="F")
        c0_src = c0[mask] if c0.size > 1 else c0
        scaled = signal / scale_dirichlet(c0_src) if mode == "dirichlet" else signal * scale_additive(c0_src)
        get_val = lambda t: scaled[t] if scaled.ndim == 1 else scaled[:, t]

        def dirichlet(t, field):
            if t < len(scaled): field[mask] = get_val(t)
            return field

        def additive(t, field):
            if t >= len(scaled): return field
            src = xp.zeros(Nx, dtype=field.dtype)
            src[mask] = get_val(t)
            return field + xp.real(xp.fft.ifft(source_kappa * xp.fft.fft(src)))

        def additive_raw(t, field):
            if t < len(scaled): field[mask] += get_val(t)
            return field

        return {"dirichlet": dirichlet, "additive": additive}.get(mode, additive_raw)

    # Pressure: dirichlet divides by c^2, additive multiplies by 2*dt/(c*dx)
    source_p_op = build_op(
        _attr(source, 'p_mask', 0), _attr(source, 'p', 0), _attr(source, 'p_mode', 'additive'),
        lambda c: c**2, lambda c: 2*dt/(c*dx))

    # Velocity: dirichlet unchanged, additive multiplies by 2*c*dt/dx
    source_u_op = build_op(
        _attr(source, 'u_mask', 0), _attr(source, 'ux', 0), _attr(source, 'u_mode', 'additive'),
        lambda c: 1, lambda c: 2*c*dt/dx)

    return source_p_op, source_u_op

def _build_absorption_ops(medium, k, rho0, Nx, xp):
    """Build absorption and dispersion operators for power-law attenuation."""
    alpha_coeff_raw = _attr(medium, 'alpha_coeff', 0)
    if not _is_enabled(alpha_coeff_raw):
        return lambda duxdx: 0, lambda rho: 0

    alpha_coeff = _expand_to_grid(alpha_coeff_raw, Nx, xp, "alpha_coeff")
    c0 = _expand_to_grid(medium.sound_speed, Nx, xp, "sound_speed")
    alpha_power = float(xp.array(_attr(medium, 'alpha_power', 1.5), dtype=float).flatten(order="F")[0])

    # Convert from dB/(MHz^y cm) to Nepers/((rad/s)^y m)
    alpha_np = 100 * alpha_coeff * (1e-6 / (2 * np.pi))**alpha_power / (20 * np.log10(np.e))
    diff = lambda f, op: _spectral_diff(f, op, xp)

    # Stokes absorption (y=2): simple viscous damping, no fractional Laplacian needed
    is_stokes_absorption = abs(alpha_power - 2.0) < 1e-10
    if is_stokes_absorption:
        return lambda duxdx: -2 * alpha_np * c0 * rho0 * duxdx, lambda rho: 0

    # General power-law: requires fractional Laplacian for causality
    tau, eta = -2 * alpha_np * c0**(alpha_power - 1), 2 * alpha_np * c0**alpha_power * xp.tan(np.pi * alpha_power / 2)
    nabla1 = _fractional_laplacian(k, alpha_power - 2, xp)
    nabla2 = _fractional_laplacian(k, alpha_power - 1, xp)

    return (lambda duxdx: tau * diff(rho0 * duxdx, nabla1),
            lambda rho: eta * diff(rho, nabla2))

# =============================================================================
# MATLAB Interop
# =============================================================================

def simulate_from_dicts(kgrid, medium, source, sensor, backend="auto"):
    """MATLAB interop entry point - converts dicts to namespace objects.

    Handles field name aliases (c0 -> sound_speed, rho0 -> density) and
    provides defaults for optional parameters.
    """
    def to_namespace(d):
        return SimpleNamespace(**dict(d))

    def normalize_medium(m):
        """Normalize field names and apply defaults."""
        d = dict(m)
        # Aliases: legacy k-Wave names -> canonical names
        if 'c0' in d and 'sound_speed' not in d:
            d['sound_speed'] = d.pop('c0')
        if 'rho0' in d and 'density' not in d:
            d['density'] = d.pop('rho0')
        return d

    return simulate(
        to_namespace(kgrid),
        to_namespace(normalize_medium(medium)),
        to_namespace(source),
        to_namespace(sensor),
        backend
    )

def interop_sanity(arr):
    """Verify MATLAB/Python data layout (for testing)."""
    a = np.array(arr, copy=True)
    a[0, 1] = 99
    return a
