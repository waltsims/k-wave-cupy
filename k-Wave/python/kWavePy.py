"""
Thin shim re-exporting from k-wave-python.

The canonical solver lives in the k-wave-python package
(kwave.solvers.kspace_solver). This module re-exports the public API
so that the MATLAB wrapper (kspaceFirstOrderPy.m) continues to work
unchanged via py.importlib.import_module('kWavePy').
"""
from kwave.solvers.kspace_solver import (  # noqa: F401
    Simulation,
    create_simulation,
    simulate_from_dicts,
    interop_sanity,
)
