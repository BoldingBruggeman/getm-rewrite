import unittest

import numpy as np

import pygetm


class Test2DTransport(unittest.TestCase):
    def check_range(self, name, values, rtol=1e-12, atol=1e-12, target_value=None):
        values = np.asarray(values)
        absrange = values.max() - values.min()
        mean = values.mean()
        crit = rtol * abs(mean) + atol
        self.assertLess(absrange, crit, "%s is not homogeneous" % name)
        if target_value is not None:
            self.assertLess(
                abs(mean - target_value),
                rtol * abs(target_value) + atol,
                "failed to match target value %s" % target_value,
            )

    def test_all(self):
        taus = [0.0, 0.01, -0.01]
        for apply_bottom_friction in [True, False]:
            for tau_x in taus:
                for tau_y in taus:
                    with self.subTest(tau_x=tau_x, tau_y=tau_y):
                        self._test(
                            tau_x=tau_x,
                            tau_y=tau_y,
                            periodic_x=tau_x != 0.0,
                            periodic_y=tau_y != 0.0,
                            apply_bottom_friction=apply_bottom_friction,
                        )

    def _test(
        self,
        periodic_x=False,
        periodic_y=False,
        tau_x=0.0,
        tau_y=0.0,
        timestep=10.0,
        ntime=360,
        apply_bottom_friction=False,
    ):
        # Set up rectangular domain (all points unmasked)
        domain = pygetm.domain.create_cartesian(
            500.0 * np.arange(100),
            500.0 * np.arange(30),
            1,
            f=0,
            H=50,
            periodic_x=periodic_x,
            periodic_y=periodic_y,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        sim = pygetm.Simulation(
            domain,
            runtype=pygetm.BAROTROPIC_2D,
            apply_bottom_friction=apply_bottom_friction,
        )

        # Idealized surface forcing
        tausx = domain.U.array(fill=tau_x)
        tausy = domain.V.array(fill=tau_y)
        sp = domain.T.array(fill=0.0)

        for _ in range(ntime):
            sim.update_surface_pressure_gradient(domain.T.z, sp)
            sim.momentum.advance_depth_integrated(
                timestep, tausx, tausy, sim.dpdx, sim.dpdy
            )
            sim.advance_surface_elevation(
                timestep, sim.momentum.U, sim.momentum.V, sim.fwf
            )
            sim.domain.update_depth()

        # Expected transport at time t: t * acceleration * water depth D
        # = t * (force/area = tau) / (mass/area = rho*D) * D = t * tau / rho
        U_ref = ntime * timestep * tau_x / pygetm.RHO0
        V_ref = ntime * timestep * tau_y / pygetm.RHO0
        U_tgt = None if apply_bottom_friction else U_ref
        V_tgt = None if apply_bottom_friction else V_ref
        self.check_range(
            "U", sim.momentum.U, target_value=U_tgt,
        )
        self.check_range(
            "V", sim.momentum.V, target_value=V_tgt,
        )
        self.check_range("z", domain.T.z, target_value=0.0)


if __name__ == "__main__":
    unittest.main()
