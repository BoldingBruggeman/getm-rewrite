import unittest

import numpy as np
import pygetm

GRAVITY = 9.81
RHO0 = 1025.0


class TestInternalPressure(unittest.TestCase):
    def test_blumberg_mellor(self):
        for ddu in (0.0, 1.0, 2.0):
            with self.subTest(ddu=ddu):
                self._test(pygetm.internal_pressure.BlumbergMellor(), ddu=ddu)

    def test_shchepetkin_mcwilliams(self):
        for ddu in (0.0, 1.0, 2.0):
            with self.subTest(ddu=ddu):
                self._test(pygetm.internal_pressure.ShchepetkinMcwilliams(), ddu=ddu)

    def _test(self, ip: pygetm.internal_pressure.Base, H=100.0, nz=30, ddu=0.0):
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
            ddu=ddu,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        domain.initialize(pygetm.BAROCLINIC)
        domain.update_depth(_3d=True)
        domain.update_depth(_3d=True)
        ip.initialize(domain)
        rho = domain.T.array(z=pygetm.CENTERS)
        buoy = domain.T.array(z=pygetm.CENTERS)

        # lock exchange density in x-direction
        rho.values[:, :, :50] = rho_min
        rho.values[:, :, 50:] = rho_max
        buoy.all_values[...] = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdy.ma == 0.0).all())
        dP_dx = (
            -domain.U.zc.values[:, 0, 0]
            * GRAVITY
            * (rho_max - rho_min)
            / domain.U.dx.values[0, 0]
        )
        acceleration = -dP_dx / RHO0
        dU = acceleration * domain.U.hn.values[:, 0, 0]
        tol = 1e-14
        diff = dU - ip.idpdx.values[:, 0, 49]
        self.assertLess(np.abs(diff).max(), tol)

        # linearly increasing density in x-direction
        rho.values[:, :, :] = rho_min + (rho_max - rho_min) * domain.T.x / 100000
        buoy.all_values[...] = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdy.ma == 0.0).all())
        diff = ip.idpdx.ma[:, 0, :] - ip.idpdx.values[:, 0, :1]
        self.assertLess(np.abs(diff).max(), tol)

        # lock exchange density in y-direction
        rho.values[:, :50, :] = rho_min
        rho.values[:, 50:, :] = rho_max
        buoy.all_values[...] = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdx.ma == 0.0).all())
        dP_dy = (
            -domain.V.zc.values[:, 0, 0]
            * GRAVITY
            * (rho_max - rho_min)
            / domain.V.dy.values[0, 0]
        )
        acceleration = -dP_dy / RHO0
        dV = acceleration * domain.V.hn.values[:, 0, 0]
        tol = 1e-14
        self.assertLess(np.abs(dV - ip.idpdy.values[:, 49, 0]).max(), tol)

        # linearly increasing density in y-direction
        rho.values[:, :, :] = rho_min + (rho_max - rho_min) * domain.T.y / 100000
        buoy.all_values[...] = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdx.ma == 0.0).all())
        diff = ip.idpdy.ma[:, :, 0] - ip.idpdy.values[:, :1, 0]
        self.assertLess(np.abs(diff).max(), tol)

        # x = np.linspace(0, 50000, 101)
        # y = np.linspace(0, 100000, 100)
        # domain = pygetm.domain.create_cartesian(
        #     x,
        #     y,
        #     nz,
        #     interfaces=True,
        #     f=0.0,
        #     H=H,
        #     ddu=ddu,
        #     logger=pygetm.parallel.get_logger(level="ERROR"),
        # )
        # dH_dx = 50 / 50000
        # dH_dy = 0.0  # 50 / 100000
        # domain.H_[...] = domain.H_ + domain.x_ * dH_dx + domain.y_ * dH_dy

        # sim = pygetm.Simulation(
        #     domain, pygetm.BAROCLINIC, internal_pressure_method=method,
        # )
        # # ((D2 * g * rho) - (D1 * g * rho))/dx

        # #
        # # stratification wih same profile everywhere (but varying water depth!)
        # sim.rho.values[:, :, :] = rho_min
        # sim.buoy.all_values[...] = (-GRAVITY / RHO0) * (sim.rho.all_values - RHO0)
        # sim.update_internal_pressure_gradient(sim.buoy)
        # self.assertTrue((sim.idpdy.ma == 0.0).all())
        # zpos = (
        #     domain.V.zf.values[-1, :, :] - domain.V.zc.values[:, :, :]
        # ) / domain.V.D.values
        # dz_dx = zpos * dH_dy
        # rho_if = rho_min + zpos * (rho_max - rho_min)
        # dP_dx = rho_if * dz_dx * GRAVITY
        # print(dP_dx, sim.idpdx.values[:, :, :])
        # acceleration = -dP_dx / RHO0
        # dU = acceleration * domain.U.hn.values[:, 0, 0]
        # tol = 1e-14
        # diff = dU - sim.idpdx.values[:, 0, 49]
        # print(diff.min(), diff.max())
        # # self.assertTrue((np.abs(diff) < tol).all())


if __name__ == "__main__":
    unittest.main()
