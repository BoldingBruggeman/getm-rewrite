import unittest
from typing import Tuple

import numpy as np
import cftime

import pygetm
import pygetm.fabm

TOLERANCE = 1e-13
START, STOP = cftime.datetime(2000, 1, 1), cftime.datetime(2000, 1, 3)
DT = (STOP - START).total_seconds()

TIMESTEP = 30
SPLIT_FACTOR = 30


class TestFABM(unittest.TestCase):
    def simulate(
        self, fabm_yaml, bioshade_feedback: bool = False
    ) -> Tuple[pygetm.domain.Domain, pygetm.Simulation]:
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
        sim = pygetm.Simulation(
            domain,
            pygetm.BAROCLINIC,
            airsea=pygetm.airsea.Fluxes(),
            radiation=pygetm.radiation.Radiation(),
            fabm=pygetm.fabm.FABM(fabm_yaml, bioshade_feedback=bioshade_feedback),
        )
        sim.start(START, TIMESTEP, SPLIT_FACTOR)
        if not bioshade_feedback:
            self.assertIsNone(sim.fabm.kc)
        while sim.time < STOP:
            sim.advance()
        sim.finish()
        return domain, sim

    def test_constant(self):
        initial = 1.0
        fabm_yaml = {
            "instances": {
                "tracer": {"model": "bb/passive", "initialization": {"c": initial}}
            }
        }
        domain, sim = self.simulate(fabm_yaml)
        self.assertLess(np.abs(sim["tracer_c"].ma / initial - 1).max(), TOLERANCE)

    def test_source(self):
        initial = 1.0
        source = 5.0 / 86400
        fabm_yaml = {
            "instances": {
                "tracer": {"model": "bb/passive", "initialization": {"c": initial}},
                "source": {
                    "model": "interior_constant",
                    "parameters": {"value": source},
                },
                "apply_source": {
                    "model": "interior_source",
                    "coupling": {"target": "tracer/c", "source": "source/data"},
                },
            }
        }
        domain, sim = self.simulate(fabm_yaml)
        expected = initial + source * DT
        self.assertLess(np.abs(sim["tracer_c"].ma / expected - 1).max(), TOLERANCE)

    def test_surface_flux(self):
        initial = 1.0
        flux = 5.0 / 86400
        fabm_yaml = {
            "instances": {
                "tracer": {"model": "bb/passive", "initialization": {"c": initial}},
                "flux": {"model": "surface_constant", "parameters": {"value": flux}},
                "apply_flux": {
                    "model": "external_surface_flux",
                    "coupling": {"target": "tracer/c", "flux": "flux/data"},
                },
            }
        }
        domain, sim = self.simulate(fabm_yaml)
        expected_total = (
            domain.T.area.ma * ((domain.T.hn.ma * initial).sum(axis=0) + flux * DT)
        ).sum()
        total = (domain.T.area.ma * domain.T.hn.ma * sim["tracer_c"].ma).sum()
        self.assertLess(np.abs(total / expected_total - 1).max(), TOLERANCE)

    def test_bottom_flux(self):
        initial = 1.0
        flux = 5.0 / 86400
        fabm_yaml = {
            "instances": {
                "tracer": {"model": "bb/passive", "initialization": {"c": initial}},
                "flux": {"model": "bottom_constant", "parameters": {"value": flux}},
                "apply_flux": {
                    "model": "external_bottom_flux",
                    "coupling": {"target": "tracer/c", "flux": "flux/data"},
                },
            }
        }
        domain, sim = self.simulate(fabm_yaml)
        expected_total = (
            domain.T.area.ma * ((domain.T.hn.ma * initial).sum(axis=0) + flux * DT)
        ).sum()
        total = (domain.T.area.ma * domain.T.hn.ma * sim["tracer_c"].ma).sum()
        self.assertLess(np.abs(total / expected_total - 1).max(), TOLERANCE)

    def test_bioshade_feedback(self):
        initial = 2.0
        specific_light_attenuation = 0.05
        fabm_yaml = {
            "instances": {
                "tracer": {
                    "model": "bb/passive",
                    "initialization": {"c": initial},
                    "parameters": {
                        "specific_light_attenuation": specific_light_attenuation
                    },
                },
            }
        }
        domain, sim = self.simulate(fabm_yaml, bioshade_feedback=True)
        expected = initial * specific_light_attenuation
        self.assertLess(np.abs(sim.fabm.kc.ma / expected - 1).max(), TOLERANCE)


if __name__ == "__main__":
    unittest.main()
