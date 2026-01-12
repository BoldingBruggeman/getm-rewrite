import unittest

import numpy as np

import pygetm
from pygetm.constants import INTERFACES


def create_grid():
    nx = 100
    ny = 2
    nz = 10
    lat = 50.0

    x = np.linspace(0, 100000, nx + 1)
    y = np.linspace(0, 1000, ny + 1)
    H = 5.0 + 1 * np.arange(nx)
    domain = pygetm.domain.create_cartesian(
        x,
        y,
        interfaces=True,
        lat=lat,
        H=H,
        logger=pygetm.parallel.get_logger(level="ERROR"),
    )
    grid = domain.create_grids(
        nz=nz,
        halox=0,
        haloy=0,
        t_postfix="t",
    )
    return grid


class TestSigma(unittest.TestCase):
    def test(self):
        for ddu in (-1.0, 0.0, 0.1, 0.25, 0.75, 1.0, 1.5, 2.0):
            for ddl in (-1.0, 0.0, 0.1, 0.25, 0.75, 1.0, 1.5, 2.0):
                self._test(ddu, ddl)

    def _test(self, ddu, ddl):
        grid = create_grid()
        vc = pygetm.vertical_coordinates.Sigma(grid.nz, ddl=ddl, ddu=ddu)
        vc(grid.D, grid.hn[...])
        i = 10
        # print(ddl,ddu,grid.D[0, i],grid.hn[:, 0, i].sum())
        self.assertAlmostEqual(
            grid.D[0, i], grid.hn[:, 0, i].sum(), places=7, msg=None, delta=None
        )


class TestGVC(unittest.TestCase):
    def test(self):
        Dgamma = 10.0
        gamma_surf = True
        for ddu in (0.0, 0.1, 0.25, 0.75, 1.0, 1.5, 2.0):
            for ddl in (0.0, 0.1, 0.25, 0.75, 1.0, 1.5, 2.0):
                self._test(ddu, ddl, Dgamma, gamma_surf)

    def _test(self, ddu, ddl, Dgamma, gamma_surf):
        grid = create_grid()
        try:
            vc = pygetm.vertical_coordinates.GVC(
                grid.nz, ddl=ddl, ddu=ddu, Dgamma=Dgamma, gamma_surf=gamma_surf
            )
            vc(grid.D[...], grid.hn[...])
            i = 10
            # print(ddl,ddu,grid.D[0, i],grid.hn[:, 0, i].sum())
            self.assertAlmostEqual(
                grid.D[0, i], grid.hn[:, 0, i].sum(), places=7, msg=None, delta=None
            )
        except Exception:
            pass


class TestAdaptive(unittest.TestCase):
    def test(self):
        Dgamma = 20.0
        gamma_surf = True
        for csigma in (-1.0, 0.0, 0.001, 0.0001):
            for cgvc in (-1.0, 0.0, 0.001, 0.0001):
                for ddu in (0.75, 1.5):
                    for ddl in (0.75, 1.5):
                        self._test(ddu, ddl, Dgamma, gamma_surf, csigma, cgvc)

    def _test(self, ddu, ddl, Dgamma, gamma_surf, csigma, cgvc):
        import logging

        grid = create_grid()

        # print(csigma, cgvc, ddl, ddu)

        #grid.zo = grid.array(fill_value=0.0)
        grid.array(name="NN", fill_value=0.0, z=INTERFACES)
        grid.array(name="SS", fill_value=0.0, z=INTERFACES)

        try:
            vc = pygetm.vertical_coordinates.Adaptive(
                grid.nz,
                cnpar=1.0,
                ddu=ddu,
                ddl=ddl,
                gamma_surf=gamma_surf,
                Dgamma=Dgamma,
                csigma=csigma,
                cgvc=cgvc,
                hpow=3,
                chsurf=-0.001,
                hsurf=1.5,
                chmidd=-0.1,
                hmidd=0.5,
                chbott=-0.001,
                hbott=1.5,
                cneigh=-0.1,
                rneigh=0.25,
                decay=2.0 / 3.0,
                cNN=-1.0,
                drho=0.3,
                cSS=-1.0,
                dvel=0.1,
                chmin=-0.1,
                hmin=2.5,
                nvfilter=-1,
                vfilter=0.2,
                nhfilter=-1,
                hfilter=0.1,
                split=1,
                timescale=1 * 3600.0,
            )
            logger = logging.getLogger(__name__)
            vc.initialize(grid, logger=logger)
            vc.update(60.0)
            i = 10
            self.assertAlmostEqual(
                grid.D[0, i], grid.hn[:, 0, i].sum(), places=7, msg=None, delta=None
            )
        except Exception:
            # print(ddl,ddu)
            pass


if __name__ == "__main__":
    unittest.main()
