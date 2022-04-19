import sys
import argparse

import numpy

import pygetm

def check_range(name, values, rtol=1e-12, atol=1e-12, target_value=None):
    values = numpy.asarray(values)
    print('  %s... ' % name, end='', flush=True)
    absrange = values.max() - values.min()
    mean = values.mean()
    crit = rtol * abs(mean) + atol
    valid = absrange < crit
    result = 'OK' if valid else 'FAILED because field is not homogeneous'
    relrange = 0 if absrange == 0 else absrange / abs(mean)
    if valid and target_value is not None:
        valid = abs(mean - target_value) < rtol * abs(target_value) + atol
        if not valid:
            result = 'FAILED to match target value %.6g' % target_value
    print('%s (mean %.6g, spread %.6g, relative spread %.6g)' % (result, mean, absrange, relrange))
    if not valid:
        plot('%s.png' % name, values)
    return valid

def plot(name, values):
    import matplotlib.pyplot
    fig, ax = matplotlib.pyplot.subplots()
    pc = ax.pcolormesh(values)
    fig.colorbar(pc)
    fig.savefig(name, dpi=300)

def test(name, periodic_x: bool=False, periodic_y: bool=False, tau_x: float=0., tau_y: float=0., timestep: float=10., ntime: int=360, apply_bottom_friction: bool=False):
    print('%s, tau_x = %s, tau_y = %s...' % (name, tau_x, tau_y), flush=True)

    # Set up rectangular domain (all points unmasked)
    domain = pygetm.domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 1, f=0, H=50, periodic_x=periodic_x, periodic_y=periodic_y)
    sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_2D, advection_scheme=pygetm.AdvectionScheme.HSIMT, apply_bottom_friction=apply_bottom_friction)

    # Idealized surface forcing
    tausx = domain.U.array(fill=tau_x)
    tausy = domain.V.array(fill=tau_y)
    sp = domain.T.array(fill=0.)

    sim.U.update_halos()
    sim.V.update_halos()
    for istep in range(ntime):
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.update_2d_momentum(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.update_sealevel(timestep, sim.U, sim.V, sim.fwf)
        sim.update_depth()

    success = check_range('U', sim.U, target_value=None if apply_bottom_friction else ntime * timestep * tau_x / pygetm.RHO0)
    success = check_range('V', sim.V, target_value=None if apply_bottom_friction else ntime * timestep * tau_y / pygetm.RHO0) and success
    success = check_range('z', domain.T.z, target_value=0) and success
    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply_bottom_friction', action='store_true')
    args = parser.parse_args()

    success = test('Periodic in x', periodic_x=True, tau_x=0.01, apply_bottom_friction=args.apply_bottom_friction)
    success = test('Periodic in x', periodic_x=True, tau_x=-0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in y', periodic_y=True, tau_y=0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in y', periodic_y=True, tau_y=-0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in x and y', periodic_x=True, periodic_y=True, tau_x=0.01, tau_y=0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in x and y', periodic_x=True, periodic_y=True, tau_x=0.01, tau_y=-0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in x and y', periodic_x=True, periodic_y=True, tau_x=-0.01, tau_y=0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in x and y', periodic_x=True, periodic_y=True, tau_x=-0.01, tau_y=-0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    if not success:
        sys.exit(1)

