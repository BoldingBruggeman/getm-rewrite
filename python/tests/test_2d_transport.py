import sys

import numpy

import pygetm
import pygetm.domain

def check_range(name, values, rtol=1e-12, atol=1e-12):
    absrange = values.max() - values.min()
    mean = values.mean()
    crit = rtol * mean + atol
    valid = absrange < crit
    if not valid:
        relrange = absrange / mean
        print('  Range check failed on %s:' % name)
        print('    absolute range: %s' % (absrange,))
        print('    relative range: %s' % (relrange,))
    return valid

def test(name, periodic_x=False, periodic_y=False, tau_x=0., tau_y=0., timestep=10., ntime=360):
    print('%s:' % name, flush=True)

    # Set up rectangular domain with outer points masked
    domain = pygetm.domain.Domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 1, f=0, H=50, periodic_x=periodic_x, periodic_y=periodic_y)
    xslice = slice(None) if periodic_x else slice(3, -3)
    yslice = slice(None) if periodic_y else slice(3, -3)
    domain.mask[yslice, xslice] = 1
    sim = pygetm.Simulation(domain, runtype=1, advection_scheme=1)

    # Idealized surface forcing
    tausx, tausx_ = domain.T.array(fill=tau_x)
    tausy, tausy_ = domain.T.array(fill=tau_y)
    sp, sp_ = domain.T.array(fill=0.)
    
    dist_U = domain.distribute(sim.momentum.U_)
    dist_V = domain.distribute(sim.momentum.V_)
    dist_U.update_halos()
    dist_V.update_halos()
    for istep in range(ntime):
        sim.pressure.surface(domain.T.z_, sp_)
        sim.momentum.uv_momentum_2d(timestep, tausx_, tausy_, sim.pressure.dpdx_, sim.pressure.dpdy_)
        dist_U.update_halos()
        dist_V.update_halos()
        sim.sealevel.update(timestep, sim.momentum.U_, sim.momentum.V_)
        domain.depth_update()

    xslice = slice(None) if periodic_x else slice(1, -1)
    yslice = slice(None) if periodic_y else slice(1, -1)
    success = check_range('U', sim.momentum.U[yslice, xslice])
    success = check_range('V', sim.momentum.V[yslice, xslice]) and success
    success = check_range('z', domain.T.z[yslice, xslice]) and success
    return success

if __name__ == '__main__':
    success = test('Periodic in x', periodic_x=True, tau_x=0.01)
    success = test('Periodic in y', periodic_y=True, tau_y=0.01) and success
    success = test('Periodic in x and y', periodic_x=True, periodic_y=True, tau_x=0.01, tau_y=0.01) and success
    if not success:
        sys.exit(1)