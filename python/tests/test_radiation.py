import unittest

import numpy as np

import pygetm

TOLERANCE = 1e-14


class TestRadiation(unittest.TestCase):
    def create_domain(self) -> pygetm.domain.Domain:
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            30,
            H=50.0,
            ddu=2.0,
            interfaces=True,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        domain.initialize(pygetm.BAROCLINIC)
        domain.update_depth(True)
        return domain

    def test_downwelling(self):
        swr_sf_value = 100.0

        def _check():
            self.assertLess(
                np.abs(rad.rad.ma - expected).max() / swr_sf_value, TOLERANCE
            )
            expected_abs = expected[1:, ...] - expected[:-1, ...]
            expected_abs[0, ...] += expected[0, ...]
            self.assertLess(
                np.abs(rad.swr_abs.ma - expected_abs).max() / swr_sf_value, TOLERANCE
            )

        domain = self.create_domain()
        rad = pygetm.radiation.TwoBand()
        rad.initialize(domain.T)
        swr_sf = domain.T.array(fill=np.nan)
        swr_sf.fill(swr_sf_value)
        self.assertRaises(AssertionError, rad, swr_sf)

        rad.kc1.fill(0.1)
        rad.kc2.fill(100.0)
        rad.A.fill(1.0)
        rad(swr_sf)
        expected = swr_sf.ma * np.exp(0.1 * domain.T.zf.ma)
        _check()

        rad.kc1.fill(100.0)
        rad.kc2.fill(0.2)
        rad.A.fill(0.0)
        rad(swr_sf)
        expected = swr_sf.ma * np.exp(0.2 * domain.T.zf.ma)
        _check()

        rad.kc1.fill(0.5)
        rad.kc2.fill(0.1)
        rad.A.fill(0.7)
        rad(swr_sf)
        expected = swr_sf.ma * (
            0.7 * np.exp(0.5 * domain.T.zf.ma) + 0.3 * np.exp(0.1 * domain.T.zf.ma)
        )
        _check()

        kc2_add = domain.T.array(z=pygetm.CENTERS, fill=np.nan)
        kc2_add.fill(0.05)
        rad(swr_sf, kc2_add)
        expected = swr_sf.ma * (
            0.7 * np.exp(0.5 * domain.T.zf.ma) + 0.3 * np.exp(0.15 * domain.T.zf.ma)
        )
        _check()

    def test_bottom_reflection(self):
        swr_sf_value = 100.0
        domain = self.create_domain()
        rad = pygetm.radiation.TwoBand(reflect_at_bottom=True)
        rad.initialize(domain.T)
        swr_sf = domain.T.array(fill=np.nan)
        swr_sf.fill(swr_sf_value)
        rad.kc1.fill(0.5)
        rad.kc2.fill(0.1)
        rad.A.fill(0.7)
        rad.bottom_albedo.fill(0.2)
        rad(swr_sf)
        down1 = swr_sf.ma * 0.7 * np.exp(0.5 * domain.T.zf.ma)
        down2 = swr_sf.ma * 0.3 * np.exp(0.1 * domain.T.zf.ma)
        expected_down = down1 + down2
        up1 = (
            0.2
            * down1[0, ...]
            * np.exp(-0.5 * (domain.T.zf.ma - domain.T.zf.ma[0, ...]))
        )
        up2 = (
            0.2
            * down2[0, ...]
            * np.exp(-0.1 * (domain.T.zf.ma - domain.T.zf.ma[0, ...]))
        )
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


if __name__ == "__main__":
    unittest.main()
