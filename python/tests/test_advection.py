import logging
import unittest
from functools import wraps

import numpy as np

import pygetm.core
import pygetm.domain
import pygetm.vertical_coordinates
import pygetm.parallel
import pygetm.operators

handler = logging.StreamHandler()
handler.setLevel(level=logging.ERROR)
logging.basicConfig(handlers=(handler,))

rng = np.random.default_rng()


def create_grid(**kwargs) -> pygetm.core.Grid:
    EXTENT = 50000
    domain = pygetm.domain.create_cartesian(
        np.linspace(0, EXTENT, 2),
        np.linspace(0, EXTENT, 2),
        f=0,
        H=50.0,
        logger=pygetm.parallel.get_logger(level="ERROR"),
        interfaces=True,
    )
    kwargs.setdefault("nz", 25)
    vc = pygetm.vertical_coordinates.Sigma(**kwargs)
    T = domain.create_grids(vc.nz, halox=0, haloy=0, velocity_grids=0)
    vc.initialize(T, logger=domain.root_logger.getChild("vertical_coordinates"))
    vc.update(0.0)
    T.ho = T.array(z=T.hn.z)
    T.ho.all_values[:, :, :] = T.hn.all_values
    return T


def for_each_grid(test_func):
    @wraps(test_func)
    def wrapper(self: unittest.TestCase, *args, **kwargs):
        for grid_args in ({}, {"ddu": 1}, {"ddl": 1}, {"ddu": 2}, {"ddl": 2}):
            with self.subTest(grid=grid_args):
                test_func(self, create_grid(**grid_args), *args, **kwargs)

    return wrapper


class TestAdvection(unittest.TestCase):
    @for_each_grid
    def test_vertical_advection_equal_thicknesses(self, grid):
        op = pygetm.operators.VerticalAdvection(grid)
        w = grid.array(z=pygetm.INTERFACES)
        w_var = grid.array(z=pygetm.INTERFACES)
        w.fill(0.0)
        TIMESTEP = 3600.0
        tracer = grid.array(z=pygetm.CENTERS)

        # Just below CFL
        w_var.fill(-0.99 * grid.ho.all_values[-1, 0, 0] / TIMESTEP)
        tracer.fill(0.0)
        tracer.all_values[-1, :, :] = 1.0
        tracer_int = (tracer.all_values * grid.ho.all_values).sum()
        op(w, w_var, TIMESTEP, tracer)
        tracer_int_after = (tracer.all_values * grid.ho.all_values).sum()
        self.assertAlmostEqual(tracer_int, tracer_int_after, places=13)

        # Above CFL
        NITER_MAX = 200
        minh = grid.ho.all_values.min(axis=0)
        for niter in range(2, NITER_MAX + 1):
            with self.subTest(niter=niter):
                w_var.fill(-niter * 0.9999 * minh / TIMESTEP)
                tracer.fill(0.0)
                tracer.all_values[-1, :, :] = 1.0
                tracer_int = (tracer.all_values * grid.ho.all_values).sum()
                op(w, w_var, TIMESTEP, tracer)
                tracer_int_after = (tracer.all_values * grid.ho.all_values).sum()
                self.assertAlmostEqual(tracer_int, tracer_int_after, places=13)
                tgt_tracer = tracer.all_values.copy()

                tracer.fill(0.0)
                tracer.all_values[-1, :, :] = 1.0
                w_var /= niter
                for _ in range(niter):
                    op(w, w_var, TIMESTEP, tracer)

                diff = np.abs(tracer.all_values[:, 0, 0] - tgt_tracer[:, 0, 0])
                self.assertAlmostEqual(diff.max(), 0.0, places=13)


if __name__ == "__main__":
    unittest.main()
