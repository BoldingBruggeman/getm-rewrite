import sys

import numpy

import pygetm

extent = 50000
domain = pygetm.domain.Domain.create_cartesian(numpy.linspace(0, extent, 50), numpy.linspace(0, extent, 52), 25, f=0, H=50)
domain.mask[...] = numpy.random.random_sample(domain.mask.shape) > 0.5
domain.H[...] = 50
sim = pygetm.Simulation(domain, runtype=2)
domain.do_vertical()

assert (domain.T.hn.all_values == domain.T.hn.all_values).all()

nuh = domain.W.array(fill=1e-2, is_3d=True)

# Set diffusivity at all masked points to NaN
nuh.all_values[:, domain.T.mask.all_values != 1] = numpy.nan

# Set diffusivity at the very surface and bottom to NaN,
# so we can later check that this value has not been propagated (used)
nuh.all_values[0, ...] = numpy.nan
nuh.all_values[-1, ...] = numpy.nan

tracer = domain.T.array(is_3d=True)
tracer.values[...] = 0
tracer.values[0, ...] = 1

ini_min = tracer.values[...].min()
ini_max = tracer.values[...].max()
vdif = pygetm.VerticalDiffusion(tracer.grid, cnpar=1.)

for _ in range(1):
    vdif.calculate(nuh, 600, tracer)

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
tol = 1e-14
if max_rel_col_range > 1e-14:
    print('ERROR: maximum range across domain %s exceeds tolerance %s' % (max_rel_col_range, tol))
    sys.exit(1)
if global_max > ini_max:
    print('ERROR: final maximum value exceeds %s initial maximum %s' % (global_max, ini_max))
    sys.exit(1)
if global_min < ini_min:
    print('ERROR: final global minimum value %s below initial minimum %s' % (global_min, ini_min))
    sys.exit(1)
    