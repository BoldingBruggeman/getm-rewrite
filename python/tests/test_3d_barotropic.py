import sys
import argparse

import numpy as np

import pygetm
import pygetm.domain
import pygetm.debug


def test(
    tau_x: float = 0.0,
    tau_y: float = 0.0,
    timestep: float = 10.0,
    ntime: int = 360,
    apply_bottom_friction: bool = False,
):
    domain = pygetm.domain.create_cartesian(
        500.0 * np.arange(100), 500.0 * np.arange(30), 50, f=0, H=50
    )
    sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D, advection_scheme=1)
    sim.logger.info(
        "Starting 3d barotropic test with tau_x=%s, tau_y=%s, apply_bottom_friction=%s"
        % (tau_x, tau_y, apply_bottom_friction)
    )

    # Idealized surface forcing
    tausx = domain.U.array(fill=tau_x)
    tausy = domain.V.array(fill=tau_y)
    sp = domain.T.array(fill=0.0)

    idpdx = domain.U.array(fill=0.0, z=pygetm.CENTERS)
    idpdy = domain.V.array(fill=0.0, z=pygetm.CENTERS)
    viscosity = domain.T.array(fill=0.0, z=pygetm.INTERFACES)

    t = domain.T.array(name="tracer", z=pygetm.CENTERS)
    rng = np.random.default_rng()
    rng.random(t.all_values.shape, out=t.all_values)
    adv = pygetm.operators.Advection(t.grid, scheme=1)

    times = timestep * np.arange(ntime)
    mode_split = 10
    z_sum_ini = domain.T.z.ma.sum()
    pre_tot = (t * domain.T.hn).values.sum()
    for istep, time in enumerate(times):
        sim.update_surface_pressure_gradient(domain.T.z, sp)
        sim.momentum.advance_depth_integrated(
            timestep, tausx, tausy, sim.dpdx, sim.dpdy
        )
        sim.advance_surface_elevation(timestep, sim.momentum.U, sim.momentum.V, sim.fwf)
        sim.domain.update_depth()

        if istep % mode_split == 0:
            sim.domain.update_depth(True)
            sim.update_surface_pressure_gradient(domain.T.zio, sp)

            sim.momentum.advance(
                timestep * mode_split,
                mode_split,
                tausx,
                tausy,
                sim.dpdx,
                sim.dpdy,
                idpdx,
                idpdy,
                viscosity,
            )

            div = np.zeros(domain.T.hn.shape)
            U1 = (sim.momentum.pk * domain.U.dy).all_values[:, 2:-2, 1:-3]
            U2 = (sim.momentum.pk * domain.U.dy).all_values[:, 2:-2, 2:-2]
            V1 = (sim.momentum.qk * domain.V.dx).all_values[:, 1:-3, 2:-2]
            V2 = (sim.momentum.qk * domain.V.dx).all_values[:, 2:-2, 2:-2]
            W1 = (sim.momentum.ww * domain.T.area).all_values[:-1, 2:-2, 2:-2]
            W2 = (sim.momentum.ww * domain.T.area).all_values[1:, 2:-2, 2:-2]
            dH = ((domain.T.hn - domain.T.ho) * domain.T.area).all_values[
                :, 2:-2, 2:-2
            ] / (timestep * mode_split)
            div = U1 - U2 + V1 - V2 + W1 - W2 - dH
            maxtp = np.zeros_like(div)
            for ar in [U1, U2, V1, V2, W1, W2, dH]:
                maxtp = np.maximum(maxtp, np.abs(ar))
            reldiv = div / np.where(maxtp > 0.0, maxtp, 1.0)
            if not pygetm.debug.check_zero(
                "maximum divergence (as missing vertical velocity in m s-1)",
                div / domain.T.area.values,
            ):
                return False
            adv.apply_3d(
                sim.momentum.uk,
                sim.momentum.vk,
                sim.momentum.ww,
                timestep * mode_split,
                t,
            )
            new_tot = (t * domain.T.hn).values.sum()
            if not pygetm.debug.check_equal(
                "layer thicknesses",
                adv.h[:, 2:-2, 2:-2],
                domain.T.hn.values,
                rtol=1e-14,
                atol=1e-14,
            ):
                return False
    return pygetm.debug.check_equal(
        "tracer total before and after simulation", new_tot, pre_tot
    ) and pygetm.debug.check_equal(
        "total volume before and after simulation",
        z_sum_ini,
        domain.T.z.ma.sum(),
        atol=1e-14,
        rtol=1e-14,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tau_x", type=float, default=0.0)
    parser.add_argument("--tau_y", type=float, default=0.0)
    parser.add_argument("--apply_bottom_friction", action="store_true")
    args = parser.parse_args()

    if not test(tau_x=args.tau_x, tau_y=args.tau_y):
        sys.exit(1)
