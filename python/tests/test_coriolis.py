import unittest

import numpy as np
import pygetm


class Test(unittest.TestCase):
    def test(self):
        for lat in (-50, -20.0, 0.0, 20.0, 50.0):
            for u in (-5, -0.0, 5.0):
                for v in (-5, -0.0, 5.0):
                    self._test(u, v, lat)

    def _test(self, u, v, lat):
        H = 100.0
        nz = 30

        x = np.linspace(0, 100000, 100)
        y = np.linspace(0, 100000, 100)
        domain = pygetm.domain.create_cartesian(
            x,
            y,
            nz,
            interfaces=True,
            lat=lat,
            H=H,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        sim = pygetm.Simulation(domain, pygetm.BAROCLINIC)
        sim.momentum.U.all_values[...] = np.where(
            domain.U.mask.all_values != 0, u * H, 0.0
        )
        sim.momentum.V.all_values[...] = np.where(
            domain.V.mask.all_values != 0, v * H, 0.0
        )

        sim.momentum.coriolis(sim.momentum.U, sim.momentum.fU)
        sim.momentum.coriolis(sim.momentum.V, sim.momentum.fV)

        OMEGA = 2.0 * np.pi / 86164.0  # 86164 is number of seconds in sidereal day

        f = 2.0 * OMEGA * np.sin(np.pi * lat / 180.0)
        self.assertTrue((domain.cor == f).all())

        fu = sim.momentum.fU / domain.V.H
        fv = sim.momentum.fV / domain.U.H

        self.assertTrue((fu.values[:-1, 1:-1] == f * u).all())
        self.assertTrue((fu.values[:-1, 0] == 0.5 * f * u).all())
        self.assertTrue((fu.values[:-1, -1] == 0.5 * f * u).all())
        # self.assertTrue((fu.values[-1, :] == 0.0).all())

        self.assertTrue((fv.values[1:-1, :-1] == f * v).all())
        self.assertTrue((fv.values[0, :-1] == 0.5 * f * v).all())
        self.assertTrue((fv.values[-1, :-1] == 0.5 * f * v).all())
        # self.assertTrue((fv.values[:, -1] == 0.0).all())

        if u == v:
            self.assertTrue((fu.values == fv.values.T).all())

        h = H / nz
        sim.momentum.pk.all_values[...] = np.where(
            domain.U.mask.all_values != 0, u * h, 0.0
        )
        sim.momentum.qk.all_values[...] = np.where(
            domain.V.mask.all_values != 0, v * h, 0.0
        )

        self.assertTrue((domain.U.hn.ma == h).all())
        self.assertTrue((domain.V.hn.ma == h).all())

        sim.momentum.coriolis(sim.momentum.pk, sim.momentum.fpk)
        sim.momentum.coriolis(sim.momentum.qk, sim.momentum.fqk)

        fu = sim.momentum.fpk / domain.V.hn
        fv = sim.momentum.fqk / domain.U.hn

        self.assertTrue((fu.values[:, :-1, 1:-1] == f * u).all())
        self.assertTrue((fu.values[:, :-1, 0] == 0.5 * f * u).all())
        self.assertTrue((fu.values[:, :-1, -1] == 0.5 * f * u).all())
        # self.assertTrue((fu.values[:, -1, :] == 0.0).all())

        self.assertTrue((fv.values[:, 1:-1, :-1] == f * v).all())
        self.assertTrue((fv.values[:, 0, :-1] == 0.5 * f * v).all())
        self.assertTrue((fv.values[:, -1, :-1] == 0.5 * f * v).all())
        # self.assertTrue((fv.values[:, :, -1] == 0.0).all())


if __name__ == "__main__":
    unittest.main()
