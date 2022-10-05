import unittest
import datetime

import numpy as np
import cftime

import pygetm


TOLERANCE = 1e-14
START, STOP = cftime.datetime(2000, 1, 1), cftime.datetime(2000, 1, 2)
DT = (STOP - START).total_seconds()

TIMESTEP = 30
SPLIT_FACTOR = 15


class Base(unittest.TestCase):
    def create_domain(self) -> pygetm.domain.Domain:
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            30,
            H=50.0,
            interfaces=True,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        domain.logger.parent.setLevel("ERROR")
        return domain

    def create_simulation(self, domain) -> pygetm.Simulation:
        sim = pygetm.Simulation(
            domain,
            pygetm.BAROCLINIC,
            airsea=pygetm.airsea.Fluxes(),
            radiation=pygetm.radiation.Radiation(),
        )
        tracer = sim.tracers.add("dum")
        sim.tracer_totals.append(pygetm.tracer.TracerTotal(tracer))
        tracer.fill(1.0)
        return sim


class TestRiversNoTracer(Base):
    def test_single_river(self):
        flow = 100.0

        domain = self.create_domain()
        river = domain.rivers.add_by_index("dummy", 25, 25)
        sim = self.create_simulation(domain)
        river.flow.set(flow)

        sim.start(START, TIMESTEP, SPLIT_FACTOR, report=datetime.timedelta(days=1))
        total_volume, total_tracers = sim.totals
        while sim.time < STOP:
            sim.advance()
        sim.finish()
        total_volume2, total_tracers2 = sim.totals

        self.assertLess(
            np.abs((total_volume2 - total_volume) / (DT * flow) - 1.0), TOLERANCE,
        )
        tt1, tot1, mean1 = total_tracers[-1]
        tt2, tot2, mean2 = total_tracers2[-1]
        self.assertLess(np.abs(tot2 / tot1 - 1.0), TOLERANCE)

    def test_collocated_rivers(self):
        flow1 = 100.0
        flow2 = 200.0
        flow = flow1 + flow2

        domain = self.create_domain()
        river1 = domain.rivers.add_by_index("dummy1", 25, 25)
        river2 = domain.rivers.add_by_index("dummy2", 25, 25)
        sim = self.create_simulation(domain)
        river1.flow.set(flow1)
        river2.flow.set(flow2)

        sim.start(START, TIMESTEP, SPLIT_FACTOR, report=datetime.timedelta(days=1))
        total_volume, total_tracers = sim.totals
        while sim.time < STOP:
            sim.advance()
        sim.finish()
        total_volume2, total_tracers2 = sim.totals

        self.assertLess(
            np.abs((total_volume2 - total_volume) / (DT * flow) - 1.0), TOLERANCE,
        )
        tt1, tot1, mean1 = total_tracers[-1]
        tt2, tot2, mean2 = total_tracers2[-1]
        self.assertLess(np.abs(tot2 / tot1 - 1.0), TOLERANCE)

    def test_multiple_rivers(self):
        n = 100
        flow = np.random.uniform(0.0, 500.0, n)

        domain = self.create_domain()
        i_all = np.random.randint(0, domain.nx, n)
        j_all = np.random.randint(0, domain.ny, n)

        rivers = []
        for iriver, (i, j) in enumerate(zip(i_all, j_all)):
            rivers.append(domain.rivers.add_by_index("dummy%i" % iriver, i, j))
        sim = self.create_simulation(domain)
        for river, f in zip(rivers, flow):
            river.flow.set(f)

        sim.start(START, TIMESTEP, SPLIT_FACTOR, report=datetime.timedelta(days=1))
        total_volume, total_tracers = sim.totals
        while sim.time < STOP:
            sim.advance()
        sim.finish()
        total_volume2, total_tracers2 = sim.totals

        self.assertLess(
            np.abs((total_volume2 - total_volume) / (DT * flow.sum()) - 1.0), TOLERANCE,
        )
        tt1, tot1, mean1 = total_tracers[-1]
        tt2, tot2, mean2 = total_tracers2[-1]
        self.assertLess(np.abs(tot2 / tot1 - 1.0), TOLERANCE)


if __name__ == "__main__":
    unittest.main()
