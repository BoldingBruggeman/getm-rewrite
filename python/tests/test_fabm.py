import unittest

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
    def setup(
        self,
        fabm_yaml,
        bioshade_feedback: bool = False,
        repair: bool = True,
        mask: bool = False,
    ) -> tuple[pygetm.domain.Domain, pygetm.Simulation]:
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            H=50.0,
            interfaces=True,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        if mask:
            domain.mask = np.random.rand(domain.ny, domain.nx) < 0.5
        sim = pygetm.Simulation(
            domain,
            airsea=pygetm.airsea.Fluxes(),
            radiation=pygetm.radiation.Radiation(),
            fabm=pygetm.fabm.FABM(
                fabm_yaml, bioshade_feedback=bioshade_feedback, repair=repair
            ),
            vertical_coordinates=pygetm.vertical_coordinates.Sigma(30, ddu=2.0),
        )
        sim.start(START, TIMESTEP, SPLIT_FACTOR)
        if not bioshade_feedback:
            self.assertIsNone(sim.fabm.kc)
        return domain, sim

    def simulate(
        self, *args, **kwargs
    ) -> tuple[pygetm.domain.Domain, pygetm.Simulation]:
        domain, sim = self.setup(*args, **kwargs)
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
            sim.T.area.ma * ((sim.T.hn.ma * initial).sum(axis=0) + flux * DT)
        ).sum()
        total = (sim.T.area.ma * sim.T.hn.ma * sim["tracer_c"].ma).sum()
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
            sim.T.area.ma * ((sim.T.hn.ma * initial).sum(axis=0) + flux * DT)
        ).sum()
        total = (sim.T.area.ma * sim.T.hn.ma * sim["tracer_c"].ma).sum()
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

    def test_check_state(self):
        fabm_yaml = {
            "instances": {
                "tracer": {"model": "bb/passive", "initialization": {"c": 1.0}}
            }
        }

        domain, sim = self.setup(fabm_yaml, repair=False)
        sim["tracer_c"].set(-1.0)
        with self.assertRaisesRegex(Exception, "FABM state contains invalid values."):
            sim.fabm.update_sources(0.0)

        domain, sim = self.setup(fabm_yaml, repair=True)
        mask = sim["tracer_c"].all_mask
        sim["tracer_c"].set(-1.0)
        sim.fabm.update_sources(0.0)
        self.assertTrue((sim["tracer_c"].ma == 0.0).all())
        self.assertTrue(
            (sim["tracer_c"].all_values[mask] == sim["tracer_c"].fill_value).all()
        )

        # Verify that check_state leaves masked values unchanged when repair=True
        # For that purpose, we set the masked values to a known fill value.
        domain, sim = self.setup(fabm_yaml, repair=True, mask=True)
        mask = sim["tracer_c"].all_mask
        FILL_CHECK = -2.0
        sim["tracer_c"].all_values.fill(FILL_CHECK)
        sim["tracer_c"].all_values[~mask] = -1.0
        sim.fabm.update_sources(0.0)
        self.assertTrue((sim["tracer_c"].ma == 0.0).all())
        self.assertTrue((sim["tracer_c"].all_values[mask] == FILL_CHECK).all())


if __name__ == "__main__":
    unittest.main()
