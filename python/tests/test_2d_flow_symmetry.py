import sys
import argparse
from typing import Optional

import numpy

import pygetm


def check_symmetry(field, name: Optional[str]=None, rtol: float=1e-12, atol: float=1e-12, mirrored=True) -> bool:
    scale = -1. if mirrored else 1.
    diff = field - scale * field[::-1, :]
    asym = diff.max() - diff.min()
    ran = field.max() - field.min()
    success = asym < ran * rtol + atol
    if name is not None:
        print('  Asymmetry in %s: %s' % (name, asym / ran), 'OK' if success else 'FAILED')
    return success

def compare(message: str, value: float, reference: float, rtol: float=1e-12, atol: float=1e-12) -> bool:
    success = abs(value - reference) < rtol * reference + atol
    print(message, 'OK' if success else 'FAILED')
    return success

def test(name: str, periodic_x: bool=False, periodic_y: bool=False, tau_x: float=0., tau_y: float=0., timestep: float=10., ntime: int=360, apply_bottom_friction: bool=False) -> bool:
    assert tau_x == 0. or tau_y == 0.
    print('%s, tau_x = %s, tau_y = %s...' % (name, tau_x, tau_y), flush=True)

    # Set up rectangular domain (all points unmasked)
    extent = 50000
    domain = pygetm.domain.Domain.create_cartesian(numpy.linspace(0, extent, 50), numpy.linspace(0, extent, 52), 1, f=0, H=50, periodic_x=periodic_x, periodic_y=periodic_y)
    distance_from_center = numpy.sqrt((domain.x - 0.5 * extent)**2 + (domain.y - 0.5 * extent)**2)

    domain.mask[distance_from_center < extent * (1. / 6. + 1e-12)] = 0
    sim = pygetm.Simulation(domain, runtype=1, advection_scheme=1, apply_bottom_friction=apply_bottom_friction)
    assert timestep < domain.maxdt, 'Request time step %s exceeds maxdt=%.5f s' % (timestep, domain.maxdt)

    assert check_symmetry(domain.mask, mirrored=False), 'Mask is not symmetric'

    # Idealized surface forcing
    tausx = domain.U.array(fill=tau_x)
    tausy = domain.V.array(fill=tau_y)
    sp = domain.T.array(fill=0.)

    sim.U.update_halos()
    sim.V.update_halos()

    V = sim.V[:-1,:] if tau_x > 0 else sim.U[:,:-1].T
    dp = sim.dpdy[:-1,:] if tau_x > 0 else sim.dpdx[:,:-1].T
    E_input, ke = 0., 0.

    # Compute initial velocities on tracer grid
    u_T = sim.U.interp(domain.T) / domain.T.D
    v_T = sim.V.interp(domain.T) / domain.T.D

    for istep in range(ntime):
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.U.update_halos()
        sim.V.update_halos()

        # Compute updated velocities on tracer grid
        u_T_old, v_T_old = u_T, v_T
        u_T = sim.U.interp(domain.T) / domain.T.D
        v_T = sim.V.interp(domain.T) / domain.T.D

        # Energy input due to wind stress (per unit area!)
        E_input += (tau_x * (u_T_old + u_T) + tau_y * (v_T_old + v_T)) * 0.5 * timestep

        sim.update_sealevel(timestep, sim.U, sim.V)
        sim.update_depth()

    E_input = (E_input * domain.T.area).ma.sum()

    # Compute total kinetic energy
    ke = sim.Ekin.ma.sum()

    success = True
    success = check_symmetry(dp, name='surface pressure gradient') and success
    success = check_symmetry(V, name='transport') and success
    success = compare('  Kinetic energy in domain vs. input by wind: %.4e J vs %.4e J' % (ke, E_input), ke, E_input, rtol=0.01) and success
    meanz = numpy.mean(domain.T.z)
    success = compare('  Mean sea level: %s m' % (meanz,), meanz, 0.)

    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply_bottom_friction', action='store_true')
    args = parser.parse_args()

    success = True
    success = test('Periodic in x', periodic_x=True, tau_x=0.01, apply_bottom_friction=args.apply_bottom_friction) and success
    success = test('Periodic in y', periodic_y=True, tau_y=0.01, apply_bottom_friction=args.apply_bottom_friction) and success

    if not success:
        sys.exit(1)

