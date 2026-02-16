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
    c_ref = float(xp.max(c0))
    k = 2 * np.pi * xp.fft.fftfreq(Nx, d=dx)
    kappa = xp.sinc((c_ref * k * dt / 2) / np.pi)
    source_kappa = xp.cos(c_ref * k * dt / 2)  # Correction for additive sources
    op_grad = 1j * k * kappa * xp.exp( 1j * k * dx/2)  # Gradient with forward shift
    op_div  = 1j * k * kappa * xp.exp(-1j * k * dx/2)  # Divergence with backward shift
    spectral_diff = lambda f, op: xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # Physics operators return 0 or identity when feature is disabled
    absorption, dispersion, nonlinearity, nonlinear_factor = _build_physics_ops(medium, k, rho0, xp)
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

        # Debug CI issue: detailed mask analysis
        p_masked = p[mask]
        if p_masked.shape[0] != sensor_data.shape[0]:
            raise ValueError(
                f"Shape mismatch at t={t}: p[mask].shape={p_masked.shape}, "
                f"sensor_data[:,t].shape={sensor_data[:,t].shape}, "
                f"sum(mask)={int(xp.sum(mask))}, mask.dtype={mask.dtype}, "
                f"mask.shape={mask.shape}, p.shape={p.shape}, "
                f"mask[:10]={mask[:10]}, mask[-10:]={mask[-10:]}"
            )
        sensor_data[:, t] = p_masked

    return {"sensor_data": _to_cpu(sensor_data), "pressure": _to_cpu(p)}

def _build_physics_ops(medium, k, rho0, xp):
    """Build composable physics operators (return 0/identity when disabled)."""
    absorption, dispersion = _build_absorption_ops(medium, k, rho0, xp)

    BonA = medium.get("BonA", 0)
    if _is_nonzero(BonA):
        # Nonlinear acoustics: pressure depends on rho^2
        nonlinearity = lambda rho: BonA * rho**2 / (2 * rho0)
        nonlinear_factor = lambda rho: (2*rho + rho0) / rho0
    else:
        nonlinearity = lambda rho: 0
        nonlinear_factor = lambda rho: 1.0

    return absorption, dispersion, nonlinearity, nonlinear_factor

def _build_source_ops(source, c0, dt, dx, Nx, source_kappa, xp):
    """Build time-varying source operators for pressure and velocity."""

    def build_source_operator(mask_raw, signal_raw, mode, scale_dirichlet, scale_additive):
        """Generic factory for source operators (pressure or velocity)."""
        if not (_is_nonzero(mask_raw) and _is_nonzero(signal_raw)):
            return lambda t, field: field

        mask = xp.array(mask_raw, dtype=bool).flatten(order="F")
        signal = xp.array(signal_raw, dtype=float).flatten(order="F")
        c0_at_source = c0[mask] if c0.size > 1 else c0

        # Scale factors from MATLAB's kspaceFirstOrder_scaleSourceTerms.m
        if mode == "dirichlet":
            scaled = signal / scale_dirichlet(c0_at_source)
        else:
            scaled = signal * scale_additive(c0_at_source)

        get_value = lambda t: scaled[t] if scaled.ndim == 1 else scaled[:, t]

        if mode == "dirichlet":
            # Replace field values at source points
            def apply_source(t, field):
                if t < len(scaled): field[mask] = get_value(t)
                return field
        elif mode == "additive":
            # Add with k-space correction for numerical stability
            def apply_source(t, field):
                if t < len(scaled):
                    src = xp.zeros(Nx, dtype=field.dtype)
                    src[mask] = get_value(t)
                    return field + xp.real(xp.fft.ifft(source_kappa * xp.fft.fft(src)))
                return field
        else:  # additive-no-correction
            # Add without k-space correction (for testing/comparison)
            def apply_source(t, field):
                if t < len(scaled): field[mask] = field[mask] + get_value(t)
                return field

        return apply_source

    # Pressure source: scale by 1/c0^2 (dirichlet) or 2*dt/(c0*dx) (additive)
    source_p_op = build_source_operator(
        source.get("p_mask", 0), source.get("p", 0), source.get("p_mode", "additive"),
        scale_dirichlet=lambda c: c**2,
        scale_additive=lambda c: 2*dt/(c*dx))

    # Velocity source: no scaling (dirichlet) or 2*c0*dt/dx (additive)
    source_u_op = build_source_operator(
        source.get("u_mask", 0), source.get("ux", 0), source.get("u_mode", "additive"),
        scale_dirichlet=lambda c: 1,
        scale_additive=lambda c: 2*c*dt/dx)

    return source_p_op, source_u_op

def _build_absorption_ops(medium, k, rho0, xp):
    """Build absorption and dispersion operators for power-law attenuation."""
    alpha_coeff = medium.get("alpha_coeff", 0)
    alpha_power = medium.get("alpha_power", 1.5)

    if not _is_nonzero(alpha_coeff):
        return lambda duxdx: 0, lambda rho: 0

    # Convert from dB/(MHz^y cm) to Nepers/((rad/s)^y m)
    alpha_np = 100 * alpha_coeff * (1e-6 / (2 * np.pi))**alpha_power / (20 * np.log10(np.e))
    c0 = xp.atleast_1d(xp.array(medium.get("sound_speed", medium.get("c0"))))

    spectral_diff = lambda f, op: xp.real(xp.fft.ifft(op * xp.fft.fft(f)))

    # Stokes absorption (y=2): simple viscous damping, no fractional Laplacian needed
    if abs(alpha_power - 2.0) < 1e-10:
        tau = -2 * alpha_np * c0
        return lambda duxdx: tau * rho0 * duxdx, lambda rho: 0

    # General power-law: requires fractional Laplacian for causality
    tau = -2 * alpha_np * c0**(alpha_power - 1)
    eta = 2 * alpha_np * c0**alpha_power * xp.tan(np.pi * alpha_power / 2)

    def fractional_laplacian(power):
        """Compute |k|^power with proper handling of k=0 singularity."""
        k_mag = xp.abs(xp.fft.fftshift(k))
        op = xp.where(xp.isinf(k_mag**power), 0, k_mag**power)
        return xp.fft.ifftshift(op)

    nabla1 = fractional_laplacian(alpha_power - 2)
    nabla2 = fractional_laplacian(alpha_power - 1)

    absorption = lambda duxdx: tau * spectral_diff(rho0 * duxdx, nabla1)
    dispersion = lambda rho: eta * spectral_diff(rho, nabla2)

    return absorption, dispersion

def interop_sanity(arr):
    """Verify MATLAB/Python data layout (for testing)."""
    a = np.array(arr, copy=True)
    a[0, 1] = 99
    return a
