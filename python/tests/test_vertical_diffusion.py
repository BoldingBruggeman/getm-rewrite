import sys

import numpy as np
import numpy.random

import pygetm

extent = 50000
domain = pygetm.domain.create_cartesian(np.linspace(0, extent, 50), np.linspace(0, extent, 52), 25, f=0, H=50)
domain.mask[...] = numpy.random.random_sample(domain.mask.shape) > 0.5   # randomly mask half of the domain
sim = pygetm.Simulation(domain, runtype=2)

assert (domain.T.ho.all_values == domain.T.hn.all_values).all(), 'ho and hn are not identical'

dt = 600
nstep = 100
cnpar = 1.

nuh = domain.T.array(fill=1e-2, z=pygetm.INTERFACES)

# Set diffusivity at all masked points to NaN
nuh.all_values[:, domain.T.mask.all_values != 1] = np.nan

# Set diffusivity at the very surface and bottom to NaN,
# so we can later check that this value has not been propagated (used)
nuh.all_values[0, ...] = np.nan
nuh.all_values[-1, ...] = np.nan

tracer = domain.T.array(z=pygetm.CENTERS)
tracer.values[...] = 0
tracer.values[0, ...] = 1

mask = np.broadcast_to(domain.T.mask.values != 1, tracer.shape)
valid_tracer = np.ma.array(tracer.values, mask=mask)

ini_min = tracer.values[...].min()
ini_max = tracer.values[...].max()
vdif = pygetm.operators.VerticalDiffusion(tracer.grid, cnpar=cnpar)

tolerance = 1e-14

for _ in range(nstep):
    vdif(nuh, dt, tracer)

if not np.isfinite(tracer)[...].all():
    print('ERROR: tracer contains non-finite values after diffusion')
    sys.exit(1)

col_min = valid_tracer.min(axis=1).min(axis=1)
col_max = valid_tracer.max(axis=1).max(axis=1)

col_range = col_max - col_min
rel_col_range = 2 * col_range / (col_min + col_max)
print('Absolute range in profiles: %s' % (col_range,))
print('Relative range in profiles: %s' % (rel_col_range,))
max_rel_col_range = rel_col_range.max()
global_min, global_max = col_min.min(), col_max.max()

if max_rel_col_range > tolerance:
    print('ERROR: maximum range across domain %s exceeds tolerance %s' % (max_rel_col_range, tolerance))
    sys.exit(1)
if global_max > ini_max:
    print('ERROR: final maximum value exceeds %s initial maximum %s' % (global_max, ini_max))
    sys.exit(1)
if global_min < ini_min:
    print('ERROR: final global minimum value %s below initial minimum %s' % (global_min, ini_min))
    sys.exit(1)

delta = np.abs(valid_tracer.sum(axis=0) - 1)
max_delta = delta.max()
if max_delta > tolerance:
    print('ERROR: difference between initial and final depth integral %s exceeds tolerance %s' % (max_delta, tolerance))
    sys.exit(1)

# Now try without spatial gradient and source term only
tracer.values[...] = 0
sources = domain.T.array(fill=1. / dt, z=pygetm.CENTERS) * dt * domain.T.hn  # note that sources should be time- and layer-integrated!
for _ in range(nstep):
    vdif(nuh, dt, tracer, ea4=sources)
delta = valid_tracer / nstep - 1
error = np.abs(delta).max()
if error > tolerance:
    print('ERROR: error %s in tracer after using vertical diffusion solver to integrate sources exceeds tolerance %s' % (error, tolerance))
    sys.exit(1)

# Now try without spatial gradient and relative [linear] source term only
tolerance = 1e-13
tracer.values[...] = 1
r = 0.1   # relative rate of increase
rel_sources = domain.T.array(fill=r / dt, z=pygetm.CENTERS) * dt * domain.T.hn   # note that sources should be time- and layer-integrated!
for _ in range(nstep):
    vdif(nuh, dt, tracer, ea2=rel_sources)
expected = 1. / (1. - r)**nstep
delta = valid_tracer - expected
rel_delta = delta / expected
rel_delta_min = rel_delta.min(axis=(1, 2))
rel_delta_max = rel_delta.max(axis=(1, 2))
if (rel_delta_min - rel_delta_max).any():
    print('ERROR: relative error in tracer varies horizontally after using vertical diffusion solver to integrate relative sources: %s vs %s' % (rel_delta_min, rel_delta_max))
rel_error = np.abs(rel_delta).max()
if rel_error > tolerance:
    print('ERROR: maximum relative error %s in tracer after using vertical diffusion solver to integrate relative sources exceeds tolerance %s' % (rel_error, tolerance))
    sys.exit(1)

# Now try without spatial gradient in tracer, but with variable diffusivity
tolerance = 1e-11
rng = np.random.default_rng()
nuh.all_values[:, :, ] = 10.**rng.uniform(-6., 0., nuh.all_values.shape)

# Set diffusivity at all masked points to NaN
nuh.all_values[:, domain.T.mask.all_values != 1] = np.nan

# Set diffusivity at the very surface and bottom to NaN,
# so we can later check that this value has not been propagated (used)
nuh.all_values[0, ...] = np.nan
nuh.all_values[-1, ...] = np.nan

constant_value = 35.
tracer.all_values[...] = np.nan
tracer.values[:, domain.T.mask.values == 1] = constant_value
for _ in range(nstep):
    vdif(nuh, dt, tracer)
error = np.abs(tracer.ma / constant_value - 1.).max()
if error > tolerance:
    print('ERROR: error %s in tracer after using vertical diffusion solver on tracer without gradient exceeds tolerance %s' % (error, tolerance))
    sys.exit(1)