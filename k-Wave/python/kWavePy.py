"""
Minimal CuPy/NumPy backend for k-Wave (1D steel thread).
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np

try:
    import cupy as cp

    _CUPY_AVAILABLE = True
except Exception:  # pragma: no cover - optional dependency
    cp = None
    _CUPY_AVAILABLE = False


def _select_xp(backend: str):
    prefer_gpu = backend in ("auto", "gpu")
    if prefer_gpu and _CUPY_AVAILABLE:
        return cp
    return np


def _to_numpy(array):
    if _CUPY_AVAILABLE and isinstance(array, cp.ndarray):
        return cp.asnumpy(array)
    return np.asarray(array)


def _flatten_f(array, xp):
    arr = xp.asarray(array)
    return arr.reshape(-1, order="F")


def _prepare_field(values, xp, length: int, dtype):
    arr = _flatten_f(values, xp)
    if arr.size == 1:
        return xp.full(length, float(arr.reshape(-1)[0]), dtype=dtype)
    if arr.size != length:
        raise ValueError(f"Expected {length} samples, got {arr.size}")
    return arr.astype(dtype, copy=False)


def _spectral_derivative(field, dx: float, xp):
    k = xp.fft.fftfreq(field.size, d=dx) * (2 * xp.pi)
    return xp.real(xp.fft.ifft(1j * k * xp.fft.fft(field)))


def interop_sanity(array):
    """
    Mutate a Python array to check MATLAB/Python index ordering.
    """
    arr = np.array(array, copy=True)
    arr[0, 1] = 99
    return arr


def simulate(
    kgrid: Dict[str, Any],
    medium: Dict[str, Any],
    source: Dict[str, Any],
    sensor: Dict[str, Any],
    backend: str = "auto",
):
    """
    Minimal 1D k-space pseudospectral leapfrog solver.
    """
    xp = _select_xp(backend)

    Nx = int(kgrid["Nx"])
    dx = float(kgrid["dx"])
    Nt = int(kgrid["Nt"])
    dt = float(kgrid["dt"])

    c0 = medium.get("sound_speed", medium.get("c0", 1500.0))
    rho0 = medium.get("density", medium.get("rho0", 1000.0))

    p0 = source.get("p0", xp.zeros(Nx))
    mask = sensor.get("mask", xp.ones(Nx, dtype=bool))

    field_dtype = xp.float64
    p = _prepare_field(p0, xp, Nx, dtype=field_dtype)
    c0_arr = _prepare_field(c0, xp, Nx, dtype=field_dtype)
    rho0_arr = _prepare_field(rho0, xp, Nx, dtype=field_dtype)
    mask_arr = _prepare_field(mask, xp, Nx, dtype=bool)

    ux = xp.zeros_like(p)
    bulk_modulus = rho0_arr * (c0_arr ** 2)
    inv_rho0 = 1.0 / rho0_arr

    sensor_idx = xp.nonzero(mask_arr)[0]
    sensor_data = xp.zeros((sensor_idx.size, Nt), dtype=field_dtype)

    for step in range(Nt):
        dpdx = _spectral_derivative(p, dx, xp)
        ux = ux - dt * inv_rho0 * dpdx

        dudx = _spectral_derivative(ux, dx, xp)
        p = p - dt * bulk_modulus * dudx

        sensor_data[:, step] = p[sensor_idx]

    return {
        "sensor_data": _to_numpy(sensor_data),
        "pressure": _to_numpy(p),
    }
