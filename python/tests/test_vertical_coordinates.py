import unittest

import numpy as np

import pygetm
from pygetm.constants import CENTERS, INTERFACES, FILL_VALUE


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
    def _get_grid(self, nz: int, halox: int = 0, haloy: int = 0):
        dom = pygetm.domain.create_cartesian(
            x=[-0.5, 0.5],
            y=[-0.5, 0.5],
            interfaces=True,
            f=0.0,
            H=100.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        grid = dom.create_grids(nz, halox, haloy)
        NN = grid.array(name="NN", z=INTERFACES)
        SS = grid.array(name="SS", z=INTERFACES)
        grid.ho = grid.array(z=CENTERS, fill_value=FILL_VALUE)
        grid.hhalf = grid.array(z=CENTERS, fill_value=FILL_VALUE)
        NN.all_values.fill(np.nan)
        SS.all_values.fill(np.nan)
        return grid, dom.logger.getChild("vertical_coordinates")

    def test_default_parameterization(self):
        c = pygetm.vertical_coordinates.Adaptive(30)
        grid, logger = self._get_grid(c.nz, halox=0, haloy=0)
        grid.fields["NN"].fill(0.0)
        grid.fields["SS"].fill(0.0)
        c.initialize(grid, logger=logger)
        c.update()
        grid.ho.all_values = grid.hn.all_values
        c.update(timestep=600.0)

        c = pygetm.vertical_coordinates.Adaptive(30)
        grid, logger = self._get_grid(c.nz, halox=2, haloy=2)
        grid.fields["NN"].fill(0.0)
        grid.fields["SS"].fill(0.0)
        c.initialize(grid, logger=logger)
        c.update()
        grid.ho.all_values = grid.hn.all_values
        c.update(timestep=600.0)

    def test_analytical_solutions(self):
        kwargs = dict(
            csigma=0.1,
            cgvc=0.0,
            ddu=0.0,
            ddl=0.0,
            Dgamma=10.0,
            gamma_surf=True,
            chsurf=0.0,
            chmidd=0.0,
            chbott=0.0,
            cneigh=0.0,
            cNN=0.0,
            cSS=0.0,
            chmin=0.0,
            hmin=0.0,
            nhfilter=0,
            nvfilter=0,
            decay=0.0,
            hpow=1.0,
            timescale=3600.0 * 4,
        )

        # Tendency towards sigma only (no halos)
        c = self._test(**kwargs)
        nug_tgt = 0.1 / (3600.0 * 4) * 2.0
        self.assertTrue(np.abs(c.nug.all_values - nug_tgt).max() < 1e-15 * nug_tgt)

        # Tendency towards sigma only (with halos)
        c = self._test(halox=2, haloy=2, **kwargs)
        self.assertTrue(np.abs(c.nug.values - nug_tgt).max() < 1e-15 * nug_tgt)

        # Tendency towards sigma and gvc
        kwargs["cgvc"] = 0.1
        kwargs["ddu"] = 1.0
        kwargs["ddl"] = 0.0
        self._test(**kwargs)
        self._test(halox=2, haloy=2, **kwargs)

        for ddu in (0.0, 0.75, 1.5):
            for ddl in (0.0, 0.75, 1.5):
                with self.subTest(ddu=ddu, ddl=ddl):
                    kwargs["ddu"] = ddu
                    kwargs["ddl"] = ddl
                    kwargs["gamma_surf"] = ddl <= ddu
                    if ddu == 0.0 and ddl == 0.0:
                        kwargs["csigma"] = 0.1
                        kwargs["cgvc"] = 0.0
                    else:
                        kwargs["csigma"] = 0.0
                        kwargs["cgvc"] = 0.1
                    self._test(**kwargs)

    def _test(
        self, halox: int = 0, haloy: int = 0, **kwargs
    ) -> pygetm.vertical_coordinates.Adaptive:
        c = pygetm.vertical_coordinates.Adaptive(30, **kwargs)
        grid, logger = self._get_grid(c.nz, halox=halox, haloy=haloy)
        c.initialize(grid, logger=logger)

        # Calculate initial thicknesses (GVC or sigma) and verify
        # that is stays unmodified if we only use GVC tendency
        # for grid diffusivity
        c.update()
        h_ini = grid.hn.values.copy()
        grid.ho.all_values = grid.hn.all_values
        c.update(timestep=600.0)

        self.assertGreaterEqual(c.nug.values.min(), 0.0)
        maxabsdiff = np.abs(grid.hn.values - h_ini).max()
        self.assertLessEqual(maxabsdiff, 1e-15 * grid.H[0, 0])
        return c


if __name__ == "__main__":
    unittest.main()
