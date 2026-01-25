import unittest
import gc

import numpy as np
import pygetm


class TestCoriolis(unittest.TestCase):
    def test(self):
        for lat in (-50, -20.0, 0.0, 20.0, 50.0):
            for u in (-5, -0.0, 5.0):
                for v in (-5, -0.0, 5.0):
                    self._test(u, v, lat)
                    gc.collect()

    def _test(self, u, v, lat):
        H = 100.0
        nz = 30

        x = np.linspace(0, 100000, 100)
        y = np.linspace(0, 100000, 100)
        domain = pygetm.domain.create_cartesian(
            x,
            y,
            interfaces=True,
            lat=lat,
            H=H,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        sim = pygetm.Simulation(
            domain, vertical_coordinates=pygetm.vertical_coordinates.Sigma(nz)
        )
        sim.momentum.U.all_values = np.where(sim.U.mask.all_values != 0, u * H, 0.0)
        sim.momentum.V.all_values = np.where(sim.V.mask.all_values != 0, v * H, 0.0)

        sim.momentum.coriolis(sim.momentum.U, sim.momentum.corV, True)
        sim.momentum.coriolis(sim.momentum.V, sim.momentum.corU, False)

        OMEGA = 2.0 * np.pi / 86164.0  # 86164 is number of seconds in sidereal day

        f = 2.0 * OMEGA * np.sin(np.pi * lat / 180.0)
        self.assertTrue((sim.T.cor.values == f).all())
        self.assertTrue((sim.U.cor.values == f).all())
        self.assertTrue((sim.V.cor.values == f).all())

        fu = sim.momentum.corV / sim.V.H
        fv = sim.momentum.corU / sim.U.H

        self.assertTrue((fu.values[:-1, 1:-1] == -f * u).all())
        self.assertTrue((fu.values[:-1, 0] == -0.5 * f * u).all())
        self.assertTrue((fu.values[:-1, -1] == -0.5 * f * u).all())
        # self.assertTrue((fu.values[-1, :] == 0.0).all())

        self.assertTrue((fv.values[1:-1, :-1] == f * v).all())
        self.assertTrue((fv.values[0, :-1] == 0.5 * f * v).all())
        self.assertTrue((fv.values[-1, :-1] == 0.5 * f * v).all())
        # self.assertTrue((fv.values[:, -1] == 0.0).all())

        if u == v:
            self.assertTrue((fu.values == -fv.values.T).all())

        h = H / nz
        sim.momentum.pk.all_values = np.where(sim.U.mask.all_values != 0, u * h, 0.0)
        sim.momentum.qk.all_values = np.where(sim.V.mask.all_values != 0, v * h, 0.0)

        self.assertTrue(
            (sim.U.hn.all_values == h).all(where=sim.U.mask.all_values != 0)
        )
        self.assertTrue(
            (sim.V.hn.all_values == h).all(where=sim.V.mask.all_values != 0)
        )

        sim.momentum.coriolis(sim.momentum.pk, sim.momentum.corqk, True)
        sim.momentum.coriolis(sim.momentum.qk, sim.momentum.corpk, False)

        fu = sim.momentum.corqk / sim.V.hn
        fv = sim.momentum.corpk / sim.U.hn

        self.assertTrue((fu.values[:, :-1, 1:-1] == -f * u).all())
        self.assertTrue((fu.values[:, :-1, 0] == 0.5 * -f * u).all())
        self.assertTrue((fu.values[:, :-1, -1] == 0.5 * -f * u).all())
        # self.assertTrue((fu.values[:, -1, :] == 0.0).all())

        self.assertTrue((fv.values[:, 1:-1, :-1] == f * v).all())
        self.assertTrue((fv.values[:, 0, :-1] == 0.5 * f * v).all())
        self.assertTrue((fv.values[:, -1, :-1] == 0.5 * f * v).all())
        # self.assertTrue((fv.values[:, :, -1] == 0.0).all())


if __name__ == "__main__":
    unittest.main()
