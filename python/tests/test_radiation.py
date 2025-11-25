import unittest
import logging

import numpy as np

import pygetm

TOLERANCE = 1e-14


class TestRadiation(unittest.TestCase):
    def create_grid(self) -> tuple[pygetm.core.Grid, logging.Logger]:
        vc = pygetm.vertical_coordinates.Sigma(30, ddu=2.0)
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            H=50.0,
            interfaces=True,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        T = domain.create_grids(nz=vc.nz, halox=2, haloy=2)
        vc.initialize(T, logger=domain.root_logger.getChild("vertical_coordinates"))
        vc.update(None)
        pygetm._pygetm.thickness2vertical_coordinates(T.mask, T.H, T.hn, T.zc, T.zf)
        return T, domain.root_logger

    def test_downwelling(self):
        swr_sf_value = 100.0

        def _check():
            self.assertLess(
                np.abs(rad.rad.ma - expected).max() / swr_sf_value, TOLERANCE
            )
            land3d = np.broadcast_to(T._land, rad.rad.all_values.shape)
            self.assertTrue(
                (rad.rad.all_values[land3d] == pygetm.constants.FILL_VALUE).all()
            )
            expected_abs = expected[1:, ...] - expected[:-1, ...]
            expected_abs[0, ...] += expected[0, ...]
            self.assertLess(
                np.abs(rad.swr_abs.ma - expected_abs).max() / swr_sf_value, TOLERANCE
            )

        T, logger = self.create_grid()
        rad = pygetm.radiation.TwoBand()
        rad.initialize(T, logger.getChild("radiation"))
        rad.par0.saved = True
        rad.par.saved = True
        swr_sf = T.array(fill=np.nan)
        swr_sf.fill(swr_sf_value)
        self.assertRaises(AssertionError, rad, swr_sf)

        rad.kc1.fill(0.1)
        rad.kc2.fill(100.0)
        rad.A.fill(1.0)
        rad(swr_sf)
        expected = swr_sf.ma * np.exp(0.1 * T.zf.ma)
        _check()
        self.assertEqual(np.abs(rad.par.ma).max(), 0.0)
        self.assertEqual(np.abs(rad.par0.ma).max(), 0.0)

        rad.kc1.fill(100.0)
        rad.kc2.fill(0.2)
        rad.A.fill(0.0)
        rad(swr_sf)
        expected = swr_sf.ma * np.exp(0.2 * T.zf.ma)
        expected_par = swr_sf.ma * np.exp(0.2 * T.zc.ma)
        _check()
        self.assertEqual(np.abs(rad.par0.ma).max(), swr_sf_value)
        self.assertLess(
            np.abs(rad.par.ma - expected_par).max() / swr_sf_value, TOLERANCE
        )

        rad.kc1.fill(0.5)
        rad.kc2.fill(0.1)
        rad.A.fill(0.7)
        rad(swr_sf)
        expected = swr_sf.ma * (
            0.7 * np.exp(0.5 * T.zf.ma) + 0.3 * np.exp(0.1 * T.zf.ma)
        )
        expected_par = swr_sf.ma * 0.3 * np.exp(0.1 * T.zc.ma)
        _check()
        self.assertLess(
            np.abs(rad.par.ma - expected_par).max() / swr_sf_value, TOLERANCE
        )
        self.assertLess(
            np.abs(rad.par0.ma - 0.3 * swr_sf_value).max() / swr_sf_value, TOLERANCE
        )

        kc2_add = T.array(z=pygetm.CENTERS, fill=np.nan)
        kc2_add.fill(0.05)
        rad(swr_sf, kc2_add)
        expected = swr_sf.ma * (
            0.7 * np.exp(0.5 * T.zf.ma) + 0.3 * np.exp(0.15 * T.zf.ma)
        )
        expected_par = swr_sf.ma * 0.3 * np.exp(0.15 * T.zc.ma)
        _check()
        self.assertLess(
            np.abs(rad.par.ma - expected_par).max() / swr_sf_value, TOLERANCE
        )
        self.assertLess(
            np.abs(rad.par0.ma - 0.3 * swr_sf_value).max() / swr_sf_value, TOLERANCE
        )

    def test_bottom_reflection(self):
        swr_sf_value = 100.0
        T, logger = self.create_grid()
        rad = pygetm.radiation.TwoBand(reflect_at_bottom=True)
        rad.initialize(T, logger.getChild("radiation"))
        rad.par0.saved = True
        rad.par.saved = True
        swr_sf = T.array(fill=np.nan)
        swr_sf.fill(swr_sf_value)
        rad.kc1.fill(0.5)
        rad.kc2.fill(0.1)
        rad.A.fill(0.7)
        rad.bottom_albedo.fill(0.2)
        rad(swr_sf)
        down1 = swr_sf.ma * 0.7 * np.exp(0.5 * T.zf.ma)
        down2 = swr_sf.ma * 0.3 * np.exp(0.1 * T.zf.ma)
        expected_down = down1 + down2
        up1 = 0.2 * down1[0, ...] * np.exp(-0.5 * (T.zf.ma - T.zf.ma[0, ...]))
        up2 = 0.2 * down2[0, ...] * np.exp(-0.1 * (T.zf.ma - T.zf.ma[0, ...]))
        expected_up = up1 + up2
        expected_abs = (expected_down[1:, ...] - expected_down[:-1, ...]) + (
            expected_up[:-1, ...] - expected_up[1:, ...]
        )
        expected_abs[0, ...] += expected_down[0, ...] - expected_up[0, ...]
        self.assertLess(
            np.abs(rad.rad.ma - expected_down).max() / swr_sf_value, TOLERANCE
        )
        self.assertLess(
            np.abs(rad.rad_up.ma - expected_up).max() / swr_sf_value, TOLERANCE
        )
        self.assertLess(
            np.abs(rad.swr_abs.ma - expected_abs).max() / swr_sf_value, TOLERANCE
        )
        land3d = np.broadcast_to(T._land, rad.rad.all_values.shape)
        self.assertTrue(
            (rad.rad.all_values[land3d] == pygetm.constants.FILL_VALUE).all()
        )

        expected_par = swr_sf.ma * 0.3 * np.exp(0.1 * T.zc.ma)
        self.assertLess(
            np.abs(rad.par.ma - expected_par).max() / swr_sf_value, TOLERANCE
        )
        self.assertLess(
            np.abs(rad.par0.ma - down2[-1, ...]).max() / swr_sf_value, TOLERANCE
        )


if __name__ == "__main__":
    unittest.main()
