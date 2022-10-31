import sys
import argparse
from typing import Optional
import unittest

import numpy as np

import pygetm


class TestFlowSymmetry(unittest.TestCase):
    def check_symmetry(
        self,
        field: pygetm.core.Array,
        name: Optional[str] = None,
        rtol: float = 1e-12,
        atol: float = 1e-12,
        mirrored: bool = True,
        axis: int = 0,
    ) -> bool:
        field = field.gather()
        if field is None:
            # Not root node
            return True

        slicespec = [slice(None)] * field.ndim
        if field.grid.type == pygetm._pygetm.UGRID:
            slicespec[-1] = slice(None, -1)
        elif field.grid.type == pygetm._pygetm.VGRID:
            slicespec[-2] = slice(None, -1)
        field = field.ma[tuple(slicespec)]
        flipspec = [slice(None)] * field.ndim
        flipspec[axis] = slice(None, None, -1)

        scale = -1.0 if mirrored else 1.0
        diff = field - scale * field[tuple(flipspec)]
        asym = diff.max() - diff.min()
        ran = field.max() - field.min()
        self.assertLess(
            asym, ran * rtol + atol, "Asymmetry in %s: %s" % (name, asym / ran)
        )

    def compare(
        self,
        message: str,
        value: float,
        reference: float,
        rtol: float = 1e-12,
        atol: float = 1e-12,
    ) -> bool:
        self.assertLess(abs(value - reference), rtol * reference + atol, message)

    def _test(
        self,
        name: str,
        periodic_x: bool = False,
        periodic_y: bool = False,
        tau_x: float = 0.0,
        tau_y: float = 0.0,
        timestep: float = 10.0,
        ntime: int = 360,
        apply_bottom_friction: bool = False,
        save: bool = False,
    ) -> bool:
        assert tau_x == 0.0 or tau_y == 0.0
        # print("%s, tau_x = %s, tau_y = %s..." % (name, tau_x, tau_y), flush=True)

        # Set up rectangular domain (all points unmasked)
        extent = 50000
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, extent, 50),
            np.linspace(0, extent, 52),
            1,
            f=0,
            H=50,
            periodic_x=periodic_x,
            periodic_y=periodic_y,
            z0=0.01 if apply_bottom_friction else 0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        # Add an island in the center of the domain
        distance_from_center = np.sqrt(
            (domain.x - 0.5 * extent) ** 2 + (domain.y - 0.5 * extent) ** 2
        )
        domain.mask[distance_from_center < extent * (1.0 / 6.0 + 1e-12)] = 0

        sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC)
        assert timestep < domain.maxdt, "Request time step %s exceeds maxdt=%.5f s" % (
            timestep,
            domain.maxdt,
        )

        if save:
            f = sim.output_manager.add_netcdf_file("island - %s.nc" % name)
            f.request("zt")

        # Idealized surface forcing
        tausx = domain.U.array(fill=tau_x)
        tausy = domain.V.array(fill=tau_y)
        sp = domain.T.array(fill=0.0)

        symmetry_axis = -2 if tau_x != 0.0 else -1
        V = sim.momentum.V if symmetry_axis == -2 else sim.momentum.U
        dp = sim.dpdy if symmetry_axis == -2 else sim.dpdx

        E_input, ke = 0.0, 0.0

        # Compute initial velocities on tracer grid
        u_T = sim.momentum.U.interp(domain.T) / domain.T.D
        v_T = sim.momentum.V.interp(domain.T) / domain.T.D

        for istep in range(ntime):
            sim.update_surface_pressure_gradient(domain.T.z, sp)
            sim.momentum.advance_depth_integrated(
                timestep, tausx, tausy, sim.dpdx, sim.dpdy
            )

            # Compute updated velocities on tracer grid
            u_T_old, v_T_old = u_T, v_T
            u_T = sim.momentum.U.interp(domain.T) / domain.T.D
            v_T = sim.momentum.V.interp(domain.T) / domain.T.D

            # Energy input due to wind stress (per unit area!)
            E_input += (
                (tau_x * (u_T_old + u_T) + tau_y * (v_T_old + v_T)) * 0.5 * timestep
            )

            sim.advance_surface_elevation(
                timestep, sim.momentum.U, sim.momentum.V, sim.fwf
            )
            sim.domain.update_depth()
            sim.output_manager.save(istep * timestep, istep)

        sim.output_manager.close(istep * timestep)

        E_input = (E_input * domain.T.area).global_sum(where=domain.T.mask != 0)

        # Compute total kinetic energy
        ke = sim.Ekin.global_sum(where=domain.T.mask != 0)

        self.check_symmetry(dp, name="surface pressure gradient", axis=symmetry_axis)
        self.check_symmetry(V, name="transport", axis=symmetry_axis)
        if E_input is not None:
            self.compare(
                "Kinetic energy in domain vs. input by wind: %.4e J vs %.4e J"
                % (ke, E_input),
                ke,
                E_input,
                rtol=0.01,
            )
        meanz = domain.T.z.global_mean(reproducible=True, where=domain.T.mask == 1)
        if meanz is not None:
            self.compare("Mean sea level: %s m" % (meanz,), meanz, 0.0)

    def test_periodic_x(self):
        for tau in (-0.01, 0.01):
            for apply_bottom_friction in (True, False):
                with self.subTest(tau=tau, apply_bottom_friction=apply_bottom_friction):
                    self._test(
                        "Periodic in x",
                        periodic_x=True,
                        tau_x=tau,
                        apply_bottom_friction=apply_bottom_friction,
                    )

    def test_periodic_y(self):
        for tau in (-0.01, 0.01):
            for apply_bottom_friction in (True, False):
                with self.subTest(tau=tau, apply_bottom_friction=apply_bottom_friction):
                    self._test(
                        "Periodic in y",
                        periodic_y=True,
                        tau_y=tau,
                        apply_bottom_friction=apply_bottom_friction,
                    )


if __name__ == "__main__":
    unittest.main()
