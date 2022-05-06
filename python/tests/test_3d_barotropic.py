from distutils.log import debug
import sys
import argparse
import datetime

import numpy

import pygetm
import pygetm.domain
import pygetm.debug

def test(tau_x: float=0., tau_y: float=0., timestep: float=10., ntime: int=360, apply_bottom_friction: bool=False):
    domain = pygetm.domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 50, f=0, H=50)
    sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D, advection_scheme=1)
    sim.logger.info('Starting 3d barotropic test with tau_x=%s, tau_y=%s, apply_bottom_friction=%s' % (tau_x, tau_y, apply_bottom_friction))

    # Idealized surface forcing
    tausx = domain.U.array(fill=tau_x)
    tausy = domain.V.array(fill=tau_y)
    sp = domain.T.array(fill=0.)

    idpdx = domain.U.array(fill=0., z=pygetm.CENTERS)
    idpdy = domain.V.array(fill=0., z=pygetm.CENTERS)
    viscosity = domain.T.array(fill=0., z=pygetm.INTERFACES)

    t = domain.T.array(name='tracer', z=pygetm.CENTERS)
    rng = numpy.random.default_rng()
    rng.random(t.all_values.shape, out=t.all_values)
    adv = pygetm.operators.Advection(t.grid, scheme=1)

    times = timestep * numpy.arange(ntime)
    mode_split = 10
    domain.T.zio.all_values[...] = 0
    domain.T.zin.all_values[...] = 0
    sim.start_3d()
    z_sum_ini = domain.T.z.ma.sum()
    pre_tot = (t * domain.T.hn).values.sum()
    for istep, time in enumerate(times):
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.update_2d_momentum(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.update_sealevel(timestep, sim.U, sim.V, sim.fwf)
        sim.domain.update_depth()

        if istep % mode_split == 0:
            sim.Ui.all_values[...] /= mode_split
            sim.Vi.all_values[...] /= mode_split
            sim.start_3d()
            sim.update_surface_pressure_gradient(domain.T.zio, sp)

            sim.update_3d_momentum(timestep * mode_split, tausx, tausy, sim.dpdx, sim.dpdy, idpdx, idpdy, viscosity)

            div = numpy.zeros(domain.T.hn.shape)
            U1 = (sim.pk * domain.U.dy).all_values[:, 2:-2, 1:-3]
            U2 = (sim.pk * domain.U.dy).all_values[:, 2:-2, 2:-2]
            V1 = (sim.qk * domain.V.dx).all_values[:, 1:-3, 2:-2]
            V2 = (sim.qk * domain.V.dx).all_values[:, 2:-2, 2:-2]
            W1 = (sim.ww * domain.T.area).all_values[:-1, 2:-2, 2:-2]
            W2 = (sim.ww * domain.T.area).all_values[1:, 2:-2, 2:-2]
            dH = ((domain.T.hn - domain.T.ho) * domain.T.area).all_values[:, 2:-2, 2:-2] / (timestep * mode_split)
            div = U1 - U2 + V1 - V2 + W1 - W2 - dH
            maxtp = numpy.zeros_like(div)
            for ar in [U1, U2, V1, V2, W1, W2, dH]:
                maxtp = numpy.maximum(maxtp, numpy.abs(ar))
            reldiv = div / numpy.where(maxtp > 0., maxtp, 1.)
            if not pygetm.debug.check_zero('maximum divergence (as missing vertical velocity in m s-1)', div / domain.T.area.values):
                return False
            adv.apply_3d(sim.uk, sim.vk, sim.ww, timestep * mode_split, t)
            new_tot = (t * domain.T.hn).values.sum()
            if not pygetm.debug.check_equal('layer thicknesses', adv.h[:,2:-2,2:-2], domain.T.hn.values, rtol=1e-14, atol=1e-14):
                return False
            sim.Ui.all_values[...] = 0
            sim.Vi.all_values[...] = 0
    return pygetm.debug.check_equal('tracer total before and after simulation', new_tot, pre_tot) and pygetm.debug.check_equal('total volume before and after simulation', z_sum_ini, domain.T.z.ma.sum(), atol=1e-14, rtol=1e-14)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tau_x', type=float, default=0.)
    parser.add_argument('--tau_y', type=float, default=0.)
    parser.add_argument('--apply_bottom_friction', action='store_true')
    args = parser.parse_args()

    if not test(tau_x=args.tau_x, tau_y=args.tau_y):
        sys.exit(1)
