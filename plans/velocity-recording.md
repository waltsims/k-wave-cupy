# Add Velocity Recording to kWavePy.py

## Context

k-Wave examples frequently record particle velocity (`sensor.record = {'p', 'ux', 'uy'}`). The Python backend currently only records pressure. Adding velocity recording completes a key feature gap and enables parity testing of velocity outputs.

## API Design

```python
# Input — sensor.record is a tuple of field names
sensor.record = ('p',)              # default: pressure only (backward compat)
sensor.record = ('p', 'ux')         # pressure + x-velocity
sensor.record = ('p', 'u')          # shorthand: pressure + all velocity components

# Output — dict keys match requested fields
result['p']            # pressure time series at sensor points (n_sensors x Nt)
result['ux']           # x-velocity time series (if requested)
result['uy']           # y-velocity time series (if requested, 2D/3D only)
result['uz']           # z-velocity time series (if requested, 3D only)
result['pressure']     # final pressure field (always included)
result['sensor_data']  # backward-compat alias for result['p']
```

`'u'` is a shorthand that expands to `('ux',)`, `('ux', 'uy')`, or `('ux', 'uy', 'uz')` based on `self.ndim`.

## File to Modify

`k-Wave/python/kWavePy.py` (lines ~111-127, ~327-340, ~399-405)

## Changes

### A. Parse `sensor.record` in `_setup_sensor_mask()` (~4 lines)

```python
record = _attr(self.sensor, 'record', ('p',))
if isinstance(record, str): record = (record,)
self.record = set(record)
if 'u' in self.record:
    self.record.update(['ux', 'uy', 'uz'][:self.ndim])
```

### B. Allocate velocity storage in `_setup_fields()` (~3 lines)

```python
self.sensor_data_u = {}
for i, v in enumerate(['ux', 'uy', 'uz'][:self.ndim]):
    if v in self.record:
        self.sensor_data_u[v] = xp.zeros((self.n_sensor_points, self.num_recorded_time_points))
```

### C. Record velocity in `step()` (~3 lines)

After existing pressure recording (line 395):
```python
for i, v in enumerate(['ux', 'uy', 'uz'][:self.ndim]):
    if v in self.sensor_data_u:
        self.sensor_data_u[v][:, file_index] = self.u[i][self.mask]
```

### D. Return results in `run()` (~4 lines)

```python
p_data = _to_cpu(self.sensor_data)
results = {"p": p_data, "sensor_data": p_data, "pressure": _to_cpu(self.p)}
for v, data in self.sensor_data_u.items():
    results[v] = _to_cpu(data)
return results
```

**Total: ~14 lines changed/added.** Backward compatible — `sensor_data` remains as alias for `p`.

## Parity Test Integration

After this plan is implemented, velocity recording scenarios will be added to `tests/test_parity.py` (see parity-test-suite.md):
- 1D: `sensor.record = ('p', 'ux')`
- 2D: `sensor.record = ('p', 'ux', 'uy')`
- 3D: `sensor.record = ('p', 'u')`

Compare velocity outputs across backends (NumPy vs CuPy, and vs MATLAB if reference data includes velocity).

## Verification

```bash
python -c "
import sys; sys.path.insert(0, 'k-Wave/python')
import kWavePy
result = kWavePy.simulate_from_dicts(
    {'Nx': 64, 'dx': 1.0, 'Nt': 100, 'dt': 6.67e-5, 'pml_size_x': 0, 'pml_alpha_x': 0},
    {'sound_speed': 1500, 'density': 1000, 'alpha_coeff': 0, 'alpha_power': 1.5, 'BonA': 0},
    {'p0': [0]*31 + [5e6] + [0]*32, 'p_mask': 0, 'p': 0, 'p_mode': 'additive',
     'u_mask': 0, 'ux': 0, 'u_mode': 'additive'},
    {'mask': [True]*64, 'record': ('p', 'ux')})
print('Keys:', sorted(result.keys()))
print('ux shape:', result['ux'].shape)
"
# Expected: Keys: ['p', 'pressure', 'sensor_data', 'ux']
#           ux shape: (64, 100)
```
