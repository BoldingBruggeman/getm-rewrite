import sys

import numpy

import pygetm

extent = 50000
domain = pygetm.domain.create_cartesian(numpy.linspace(0, extent, 50), numpy.linspace(0, extent, 52), 25, f=0, H=50)
domain.mask[...] = numpy.random.random_sample(domain.mask.shape) > 0.5   # randomly mask half of the domain
sim = pygetm.Simulation(domain, runtype=2)
domain.do_vertical()

assert (domain.T.ho.all_values == domain.T.hn.all_values).all(), 'ho and hn are not identical'

dt = 600
nstep = 100
cnpar = 1.

nuh = domain.W.array(fill=1e-2, is_3d=True)

# Set diffusivity at all masked points to NaN
nuh.all_values[:, domain.W.mask.all_values != 1] = numpy.nan

# Set diffusivity at the very surface and bottom to NaN,
# so we can later check that this value has not been propagated (used)
nuh.all_values[0, ...] = numpy.nan
nuh.all_values[-1, ...] = numpy.nan

tracer = domain.T.array(is_3d=True)
tracer.values[...] = 0
tracer.values[0, ...] = 1

ini_min = tracer.values[...].min()
ini_max = tracer.values[...].max()
vdif = pygetm.VerticalDiffusion(tracer.grid, cnpar=cnpar)

for _ in range(nstep):
    vdif(nuh, dt, tracer)

if not numpy.isfinite(tracer)[...].all():
    print('ERROR: tracer contains non-finite values after diffusion')
    sys.exit(1)
mask = numpy.broadcast_to(domain.T.mask.values != 1, tracer.shape)
valid_tracer = numpy.ma.array(tracer.values, mask=mask)

col_min = valid_tracer.min(axis=1).min(axis=1)
col_max = valid_tracer.max(axis=1).max(axis=1)

col_range = col_max - col_min
rel_col_range = 2 * col_range / (col_min + col_max)
print('Absolute range in profiles: %s' % (col_range,))
print('Relative range in profiles: %s' % (rel_col_range,))
max_rel_col_range = rel_col_range.max()
global_min, global_max = col_min.min(), col_max.max()
tolerance = 1e-14
if max_rel_col_range > tolerance:
    print('ERROR: maximum range across domain %s exceeds tolerance %s' % (max_rel_col_range, tolerance))
    sys.exit(1)
if global_max > ini_max:
    print('ERROR: final maximum value exceeds %s initial maximum %s' % (global_max, ini_max))
    sys.exit(1)
if global_min < ini_min:
    print('ERROR: final global minimum value %s below initial minimum %s' % (global_min, ini_min))
    sys.exit(1)

delta = numpy.abs(valid_tracer.sum(axis=0) - 1)
max_delta = delta.max()
if max_delta > tolerance:
    print('ERROR: difference between initial and final depth integral %s exceeds tolerance %s' % (max_delta, tolerance))
    sys.exit(1)

# Now try without spatial gradient and source term only
tracer.values[...] = 0
sources = domain.T.array(fill=1. / dt, is_3d=True) * dt * domain.T.hn  # note that sources should be time- and layer-integrated!
for _ in range(nstep):
    vdif(nuh, dt, tracer, ea4=sources)
delta = valid_tracer / nstep - 1
error = numpy.abs(delta).max()
if error > tolerance:
    print('ERROR: error %s in tracer after using vertical diffusion solver to integrate sources exceeds tolerance %s' % (error, tolerance))
    sys.exit(1)

# Now try without spatial gradient and source term only
tracer.values[...] = 1
sources = domain.T.array(fill=-.1 / dt, is_3d=True) * dt * domain.T.hn   # note that sources should be time- and layer-integrated!
for _ in range(1):
    vdif(nuh, dt, tracer, ea2=sources)
print(valid_tracer)
