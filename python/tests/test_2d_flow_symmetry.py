import sys
import argparse
from typing import Optional

import numpy

import pygetm

class Debug:
    def check_symmetry(self, field: pygetm.core.Array, name: Optional[str]=None, rtol: float=1e-12, atol: float=1e-12, mirrored: bool=True, axis: int=0) -> bool:
        field = field.gather()
        if field is None:
            # Not root node
            return True

        slicespec = [slice(None)] * field.ndim
        if field.grid.type == pygetm._pygetm.UGRID:
            slicespec[-1] = slice(None, -1)
        elif field.grid.type == pygetm._pygetm.VGRID:
            slicespec[-2] = slice(None, -1)
        field = field[tuple(slicespec)]
        flipspec = [slice(None)] * field.ndim
        flipspec[axis] = slice(None, None, -1)

        scale = -1. if mirrored else 1.
        diff = field - scale * field[tuple(flipspec)]
        asym = diff.max() - diff.min()
        ran = field.max() - field.min()
        success = asym < ran * rtol + atol
        if name is not None:
            print('  Asymmetry in %s: %s' % (name, asym / ran), 'OK' if success else 'FAILED')
        return success

    def compare(self, message: str, value: float, reference: float, rtol: float=1e-12, atol: float=1e-12) -> bool:
        success = abs(value - reference) < rtol * reference + atol
        print(message, 'OK' if success else 'FAILED')
        return success

def test(name: str, periodic_x: bool=False, periodic_y: bool=False, tau_x: float=0., tau_y: float=0., timestep: float=10., ntime: int=360, apply_bottom_friction: bool=False, save: bool= False) -> bool:
    assert tau_x == 0. or tau_y == 0.
    print('%s, tau_x = %s, tau_y = %s...' % (name, tau_x, tau_y), flush=True)

    # Set up rectangular domain (all points unmasked)
    extent = 50000
    domain = pygetm.domain.Domain.create_cartesian(numpy.linspace(0, extent, 50), numpy.linspace(0, extent, 52), 1, f=0, H=50, periodic_x=periodic_x, periodic_y=periodic_y)
    distance_from_center = numpy.sqrt((domain.x - 0.5 * extent)**2 + (domain.y - 0.5 * extent)**2)

    domain.mask[distance_from_center < extent * (1. / 6. + 1e-12)] = 0
    sim = pygetm.Simulation(domain, runtype=1, advection_scheme=1, apply_bottom_friction=apply_bottom_friction)
    assert timestep < domain.maxdt, 'Request time step %s exceeds maxdt=%.5f s' % (timestep, domain.maxdt)

    if save:
        f = sim.output_manager.add_netcdf_file('island - %s.nc' % name)
        f.request('zt')

    deb = Debug()
    #assert deb.check_symmetry(domain.mask, mirrored=False), 'Mask is not symmetric'

    # Idealized surface forcing
    tausx = domain.U.array(fill=tau_x)
    tausy = domain.V.array(fill=tau_y)
    sp = domain.T.array(fill=0.)

    sim.U.update_halos()
    sim.V.update_halos()

    symmetry_axis = -2 if tau_x > 0 else -1
    V = sim.V if symmetry_axis == -2 else sim.U
    dp = sim.dpdy if symmetry_axis == -2 else sim.dpdx

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
        sim.output_manager.save()

    sim.output_manager.close()

    E_input = (E_input * domain.T.area).global_sum()

    # Compute total kinetic energy
    ke = sim.Ekin.global_sum()

    success = True
    success = deb.check_symmetry(dp, name='surface pressure gradient', axis=symmetry_axis) and success
    success = deb.check_symmetry(V, name='transport', axis=symmetry_axis) and success
    if E_input is not None:
        success = deb.compare('  Kinetic energy in domain vs. input by wind: %.4e J vs %.4e J' % (ke, E_input), ke, E_input, rtol=0.01) and success
    meanz = domain.T.z.global_mean(reproducible=True, where=domain.T.mask == 1)
    if meanz is not None:
        success = deb.compare('  Mean sea level: %s m' % (meanz,), meanz, 0.)

    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply_bottom_friction', action='store_true')
    parser.add_argument('--save', action='store_true')
    args = parser.parse_args()

    success = True
    success = test('Periodic in x', periodic_x=True, tau_x=0.01, apply_bottom_friction=args.apply_bottom_friction, save=args.save) and success
    success = test('Periodic in y', periodic_y=True, tau_y=0.01, apply_bottom_friction=args.apply_bottom_friction, save=args.save) and success

    if not success:
        sys.exit(1)

