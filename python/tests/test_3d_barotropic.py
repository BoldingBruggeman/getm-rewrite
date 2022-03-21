import sys
import argparse
import datetime

import numpy

import pygetm
import pygetm.domain

def test(tau_x: float=0., tau_y: float=0., timestep: float=10., ntime: int=360, apply_bottom_friction: bool=False):
    domain = pygetm.domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 50, f=0, H=50)
    sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D, advection_scheme=1)

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
    adv = pygetm.Advection(t.grid, scheme=1)

    times = timestep * numpy.arange(ntime)
    mode_split = 10
    domain.T.zio.all_values[...] = 0
    domain.T.zin.all_values[...] = 0
    sim.start(datetime.datetime(2000, 1, 1), timestep, mode_split)
    pre_tot = (t * domain.T.hn * domain.T.area).values.sum()
    for istep, time in enumerate(times):
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
        sim.update_sealevel(timestep, sim.U, sim.V)
        sim.update_depth()

        if istep % mode_split == 0:
            sim.Ui.all_values[...] /= mode_split
            sim.Vi.all_values[...] /= mode_split
            sim.start_3d()
            domain.do_vertical()
            sim.update_surface_pressure_gradient(domain.T.zio, sp)

            sim.uvw_momentum_3d(timestep * mode_split, tausx, tausy, sim.dpdx, sim.dpdy, idpdx, idpdy, viscosity)

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
            missing_w = div / domain.T.area.values
            max_missing_w = numpy.abs(missing_w).max()
            print('maximum divergence (as missing vertical velocity in m s-1): %.6e ' % max_missing_w, end='')
            success = max_missing_w < 1e-15
            print('OK' if success else 'FAILED')
            if not success:
                return False
            #for k in range(div.shape[0]):
            #    print(reldiv[k, ...].min(), reldiv[k, ...].max())
            #print(reldiv.min(), reldiv.max(), div.min(), div.max())
            #print(reldiv[-1, ...].min(), reldiv[-1, ...].max())
            #pre_tot = (t * domain.T.hn).values.sum()
            #adv.apply_3d(sim.uk, sim.vk, sim.ww, timestep * mode_split, t)
            #print((t * adv.h).values.sum() / pre_tot - 1)
            sim.Ui.all_values[...] = 0
            sim.Vi.all_values[...] = 0
        #print((t * domain.T.hn * domain.T.area).values.sum() / pre_tot - 1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--apply_bottom_friction', action='store_true')
    args = parser.parse_args()

    if not test(tau_x=0.01):
        sys.exit(1)

