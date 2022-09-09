import unittest

import numpy as np
import pygetm

GRAVITY = 9.81
RHO0 = 1025.0


class Test(unittest.TestCase):
    def test_blumberg_mellor(self):
        self._test(pygetm.InternalPressure.BLUMBERG_MELLOR)

    def test_shchepetkin_mcwilliams(self):
        self._test(pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS)

    def _test(self, method: pygetm.InternalPressure, H=100.0, nz=30):
        rho_min = 1020.0
        rho_max = 1025.0

        x = np.linspace(0, 50000, 101)
        y = np.linspace(0, 100000, 100)
        domain = pygetm.domain.create_cartesian(
            x,
            y,
            nz,
            interfaces=True,
            f=0.0,
            H=H,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        sim = pygetm.Simulation(
            domain, pygetm.BAROCLINIC, internal_pressure_method=method,
        )

        # lock exchange density in x direction
        sim.rho.values[:, :, :50] = rho_min
        sim.rho.values[:, :, 50:] = rho_max
        sim.buoy.all_values[...] = (-GRAVITY / RHO0) * (sim.rho.all_values - RHO0)
        sim.update_internal_pressure_gradient(
            sim.buoy, sim.momentum.SxB, sim.momentum.SyB
        )
        self.assertTrue((sim.idpdy.ma == 0.0).all())
        dP_dx = (
            -domain.U.zc.values[:, 0, 0]
            * GRAVITY
            * (rho_max - rho_min)
            / domain.U.dx.values[0, 0]
        )
        acceleration = -dP_dx / RHO0
        dU = acceleration * domain.U.hn.values[:, 0, 0]
        tol = 1e-15
        self.assertTrue((np.abs(dU - sim.idpdx.values[:, 0, 49]) < tol).all())

        # linearly increasing density in x direction
        sim.rho.values[:, :, :] = rho_min + (rho_max - rho_min) * domain.T.x / 100000
        sim.buoy.all_values[...] = (-GRAVITY / RHO0) * (sim.rho.all_values - RHO0)
        sim.update_internal_pressure_gradient(
            sim.buoy, sim.momentum.SxB, sim.momentum.SyB
        )
        self.assertTrue((sim.idpdy.ma == 0.0).all())
        self.assertTrue(
            (np.abs(sim.idpdx.ma[:, 0, :] - sim.idpdx.values[:, 0, :1]) < tol).all()
        )

        # lock exchange density in y direction
        sim.rho.values[:, :50, :] = rho_min
        sim.rho.values[:, 50:, :] = rho_max
        sim.buoy.all_values[...] = (-GRAVITY / RHO0) * (sim.rho.all_values - RHO0)
        sim.update_internal_pressure_gradient(
            sim.buoy, sim.momentum.SxB, sim.momentum.SyB
        )
        self.assertTrue((sim.idpdx.ma == 0.0).all())
        dP_dy = (
            -domain.V.zc.values[:, 0, 0]
            * GRAVITY
            * (rho_max - rho_min)
            / domain.V.dy.values[0, 0]
        )
        acceleration = -dP_dy / RHO0
        dV = acceleration * domain.V.hn.values[:, 0, 0]
        tol = 1e-15
        self.assertTrue((np.abs(dV - sim.idpdy.values[:, 49, 0]) < tol).all())

        # linearly increasing density in y direction
        sim.rho.values[:, :, :] = rho_min + (rho_max - rho_min) * domain.T.y / 100000
        sim.buoy.all_values[...] = (-GRAVITY / RHO0) * (sim.rho.all_values - RHO0)
        sim.update_internal_pressure_gradient(
            sim.buoy, sim.momentum.SxB, sim.momentum.SyB
        )
        self.assertTrue((sim.idpdx.ma == 0.0).all())
        self.assertTrue(
            (np.abs(sim.idpdy.ma[:, :, 0] - sim.idpdy.values[:, :1, 0]) < tol).all()
        )


if __name__ == "__main__":
    unittest.main()
