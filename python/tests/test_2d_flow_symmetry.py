import sys
import argparse

import numpy

import pygetm


def check_symmetry(name, field, rtol=1e-12, atol=1e-12):
    asym = (field + field[::-1, :]).max() - (field + field[::-1, :]).min()
    ran = field.max() - field.min()
    print('  Asymmetry in %s: %s' % (name, asym / ran))
    return asym < ran * rtol + atol

def test(name, periodic='X', tau=0., timestep=10., ntime=360, apply_bottom_friction=False):
    print('%s, tau = %s...' % (name, tau), flush=True)

    # Set up rectangular domain (all points unmasked)
    domain = pygetm.domain.Domain.create_cartesian(numpy.linspace(0, 50000, 100), numpy.linspace(0, 50000, 102), 1, f=0, H=50, periodic_x=periodic == 'X', periodic_y=periodic == 'Y')
    distance_from_center = numpy.sqrt((domain.x - 25000)**2 + (domain.y - 25000)**2)
    domain.mask[distance_from_center < 7500] = 0
    sim = pygetm.Simulation(domain, runtype=1, advection_scheme=1, apply_bottom_friction=apply_bottom_friction)
    assert timestep < domain.maxdt

    # Idealized surface forcing
    tausx, tausx_ = domain.U.array(fill=tau if periodic == 'X' else 0.)
    tausy, tausy_ = domain.V.array(fill=tau if periodic == 'Y' else 0.)
    sp, sp_ = domain.T.array(fill=0.)

    dist_U = domain.distribute(sim.U_)
    dist_V = domain.distribute(sim.V_)
    dist_U.update_halos()
    dist_V.update_halos()
    V = sim.V[:-1,:] if periodic == 'X' else sim.U[:,:-1].T
    dp = sim.dpdy[:-1,:] if periodic == 'X' else sim.dpdx[:,:-1].T
    for istep in range(ntime):
        sim.update_surface_pressure_gradient(domain.T.z_, sp_)
        sim.uv_momentum_2d(timestep, tausx_, tausy_, sim.dpdx_, sim.dpdy_)
        dist_U.update_halos()
        dist_V.update_halos()
        sim.update_sealevel(timestep, sim.U_, sim.V_)
        sim.update_depth()
    
    success = True
    success = check_symmetry('surface pressure gradient', dp) and success
    success = check_symmetry('transport', V) and success

    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply_bottom_friction', action='store_true')
    args = parser.parse_args()

    success = True
    success = test('Periodic in x', periodic='X', tau=0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in y', periodic='Y', tau=0.01, apply_bottom_friction=args.apply_bottom_friction) and success

    if not success:
        sys.exit(1)

