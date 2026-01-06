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
            interfaces=True,
            f=0.0,
            H=H,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        vc = pygetm.vertical_coordinates.Sigma(nz, ddu=ddu)
        T = domain.create_grids(vc.nz, halox=2, haloy=2, velocity_grids=1)
        U, V = T.ugrid, T.vgrid
        vc.initialize(
            T, U, V, logger=domain.root_logger.getChild("vertical_coordinates")
        )
        vc.update(None)
        for g in (T, U, V):
            pygetm._pygetm.thickness2vertical_coordinates(g.mask, g.H, g.hn, g.zc, g.zf)

        ip.initialize(U, V)
        rho = T.array(z=pygetm.CENTERS)
        buoy = T.array(z=pygetm.CENTERS)

        # lock exchange density in x-direction
        rho.values[:, :, :50] = rho_min
        rho.values[:, :, 50:] = rho_max
        buoy.all_values = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdy.ma == 0.0).all())
        dP_dx = (
            -U.zc.values[:, 0, 0] * GRAVITY * (rho_max - rho_min) / U.dx.values[0, 0]
        )
        acceleration = -dP_dx / RHO0
        dU = acceleration * U.hn.values[:, 0, 0]
        tol = 1e-14
        diff = dU - ip.idpdx.values[:, 0, 49]
        self.assertLess(np.abs(diff).max(), tol)

        # linearly increasing density in x-direction
        rho.values[:, :, :] = rho_min + (rho_max - rho_min) * T.x / 100000
        buoy.all_values = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdy.ma == 0.0).all())
        diff = ip.idpdx.ma[:, 0, :] - ip.idpdx.values[:, 0, :1]
        self.assertLess(np.abs(diff).max(), tol)

        # lock exchange density in y-direction
        rho.values[:, :50, :] = rho_min
        rho.values[:, 50:, :] = rho_max
        buoy.all_values = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdx.ma == 0.0).all())
        dP_dy = (
            -V.zc.values[:, 0, 0] * GRAVITY * (rho_max - rho_min) / V.dy.values[0, 0]
        )
        acceleration = -dP_dy / RHO0
        dV = acceleration * V.hn.values[:, 0, 0]
        tol = 1e-14
        self.assertLess(np.abs(dV - ip.idpdy.values[:, 49, 0]).max(), tol)

        # linearly increasing density in y-direction
        rho.values[:, :, :] = rho_min + (rho_max - rho_min) * T.y / 100000
        buoy.all_values = (-GRAVITY / RHO0) * (rho.all_values - RHO0)
        ip(buoy)
        self.assertTrue((ip.idpdx.ma == 0.0).all())
        diff = ip.idpdy.ma[:, :, 0] - ip.idpdy.values[:, :1, 0]
        self.assertLess(np.abs(diff).max(), tol)

        # elevation gradient in x and constant buoyancy
        # analytical solution: idpdx / hu = buoy * d(elev)/dx
        lin = np.linspace(0.5, 1.5, x.size - 1)
        h_bck = T.hn.values[-1, ...].copy()
        T.hn.values[-1, ...] += lin
        T.hn.interp(U.hn)
        T.hn.interp(V.hn)
        pygetm._pygetm.thickness2vertical_coordinates(T.mask, T.H, T.hn, T.zc, T.zf)

        for const_buoy in (0.01, 0.0, -0.01):
            with self.subTest(const_buoy=const_buoy):
                buoy.all_values = const_buoy
                ip(buoy)
                idp_per_h = ip.idpdx.ma / U.hn.ma
                d_elev_dx = (lin[1] - lin[0]) / U.dx.values[0, 0]
                target = const_buoy * d_elev_dx
                self.assertLessEqual(
                    np.abs(idp_per_h.min() - target), 1e-12 * np.abs(target)
                )
                self.assertLessEqual(
                    np.abs(idp_per_h.max() - target), 1e-12 * np.abs(target)
                )

        # elevation gradient in y and constant buoyancy
        # analytical solution: idpdx / hu = buoy * d(elev)/dx
        lin = np.linspace(0.5, 1.5, y.size - 1)
        T.hn.values[-1, ...] = h_bck + lin[:, np.newaxis]
        T.hn.interp(U.hn)
        T.hn.interp(V.hn)
        pygetm._pygetm.thickness2vertical_coordinates(T.mask, T.H, T.hn, T.zc, T.zf)

        for const_buoy in (0.01, 0.0, -0.01):
            with self.subTest(const_buoy=const_buoy):
                buoy.all_values = const_buoy
                ip(buoy)
                idp_per_h = ip.idpdy.ma / V.hn.ma
                d_elev_dy = (lin[1] - lin[0]) / V.dy.values[0, 0]
                target = const_buoy * d_elev_dy
                self.assertLessEqual(
                    np.abs(idp_per_h.min() - target), 1e-12 * np.abs(target)
                )
                self.assertLessEqual(
                    np.abs(idp_per_h.max() - target), 1e-12 * np.abs(target)
                )


if __name__ == "__main__":
    unittest.main()
