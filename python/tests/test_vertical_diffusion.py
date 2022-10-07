import unittest
import logging
from functools import wraps
from typing import Optional

import numpy as np

import pygetm


handler = logging.StreamHandler()
handler.setLevel(level=logging.ERROR)
logging.basicConfig(handlers=(handler,))

rng = np.random.default_rng()


def for_each_grid(test_func):
    @wraps(test_func)
    def wrapper(self: unittest.TestCase, *args, **kwargs):
        for grid_args in ({}, {"ddu": 1}, {"ddl": 1}):
            with self.subTest(grid=grid_args):
                EXTENT = 50000
                domain = pygetm.domain.create_cartesian(
                    np.linspace(0, EXTENT, 50),
                    np.linspace(0, EXTENT, 52),
                    25,
                    f=0,
                    H=50.0,
                    **grid_args
                )
                # randomly mask half of the domain
                domain.mask[...] = rng.random(domain.mask.shape) > 0.5
                self.sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D)
                assert (domain.T.ho.all_values == domain.T.hn.all_values).all()
                test_func(self, domain.T, *args, **kwargs)

    return wrapper


def for_each_cnpar(test_func):
    @wraps(test_func)
    def wrapper(self: unittest.TestCase, *args, **kwargs):
        for cnpar in [0.5, 0.75, 1.0]:
            with self.subTest(cnpar=cnpar):
                test_func(self, cnpar, *args, **kwargs)

    return wrapper


def repeat(test_func):
    @wraps(test_func)
    def wrapper(self: unittest.TestCase, *args, **kwargs):
        for i in range(10):
            with self.subTest(repeat=i):
                test_func(self, *args, **kwargs)

    return wrapper


class TestVerticalDiffusion(unittest.TestCase):
    DT = 600.0
    NSTEP = 100

    def diffuse(
        self,
        tracer_in: pygetm.core.Array,
        nuh: pygetm.core.Array,
        cnpar: float,
        tolerance: float = 1e-13,
        sources: Optional[pygetm.core.Array] = None,
        change_ho: bool = False,
    ):
        self.assertTrue(
            np.isfinite(tracer_in.ma).all(),
            "tracer contains non-finite values before diffusion",
        )

        tracer = tracer_in.grid.array(z=tracer_in.z)
        tracer.all_values[...] = tracer_in.all_values
        vdif = pygetm.operators.VerticalDiffusion(tracer.grid, cnpar=cnpar)

        # Set diffusivity at all masked points to NaN
        nuh.all_values[:, nuh.grid.mask.all_values != 1] = np.nan

        # Set diffusivity at the very surface and bottom to NaN,
        # so we can later check that this value has not been propagated (used)
        nuh.all_values[0, ...] = np.nan
        nuh.all_values[-1, ...] = np.nan

        expected_integral = (tracer_in.ma * tracer_in.grid.hn).sum(axis=0)
        for _ in range(self.NSTEP):
            if change_ho:
                tracer.grid.ho.all_values[...] = tracer.grid.hn.all_values
                hn_pert = np.random.lognormal(
                    np.log(tracer.grid.hn.ma.mean(axis=(1, 2))), 0.1
                )
                tracer.grid.hn.fill(hn_pert[:, np.newaxis, np.newaxis])
            vdif(nuh, self.DT, tracer, ea4=sources, use_ho=change_ho)

        self.assertTrue(
            np.isfinite(tracer.ma)[...].all(),
            "tracer contains non-finite values after diffusion",
        )

        col_min = tracer.ma.min(axis=(1, 2))
        col_max = tracer.ma.max(axis=(1, 2))

        col_range = col_max - col_min
        self.assertEqual(
            col_range.max(), 0.0, "horizontal variability in tracer after diffusion"
        )
        if sources is None and not change_ho:
            ini_min, ini_max = tracer_in.ma.min(), tracer_in.ma.max()
            global_min, global_max = col_min.min(), col_max.max()
            eps = tolerance * max(abs(global_min), abs(global_max))
            self.assertLessEqual(
                global_max,
                ini_max + eps,
                "final global maximum value exceeds initial maximum",
            )
            self.assertGreaterEqual(
                global_min,
                ini_min - eps,
                "final global minimum value below initial minimum",
            )

        if sources is not None:
            expected_integral += sources.values.sum(axis=0) * self.NSTEP
        integral = (tracer.ma * tracer.grid.hn).sum(axis=0)
        delta = np.abs(integral - expected_integral).max()
        reldelta = delta / np.abs(expected_integral).max()
        self.assertLessEqual(
            reldelta, tolerance, "depth integral differs from expected value"
        )

    @for_each_grid
    @for_each_cnpar
    def test_mixing_from_bottom(self, cnpar: float, grid: pygetm.domain.Grid):
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        tracer = grid.array(z=pygetm.CENTERS, fill=0.0)
        tracer.values[0, ...] = 1.0
        self.diffuse(tracer, nuh, cnpar)

    @for_each_grid
    @for_each_cnpar
    def test_mixing_from_surface(self, cnpar: float, grid: pygetm.domain.Grid):
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        tracer = grid.array(z=pygetm.CENTERS, fill=0.0)
        tracer.values[-1, ...] = 1.0
        self.diffuse(tracer, nuh, cnpar)

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_mixing_of_random_state(self, cnpar: float, grid: pygetm.domain.Grid):
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        tracer = grid.array(z=pygetm.CENTERS)
        tracer[...] = rng.uniform(0.0, 1.0, (tracer.shape[0], 1, 1))
        self.diffuse(tracer, nuh, cnpar)

    @for_each_grid
    @for_each_cnpar
    def test_source(self, cnpar: float, grid: pygetm.domain.Grid):
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        tracer = grid.array(z=pygetm.CENTERS, fill=0.0)
        # note that sources should be time- and layer-integrated!
        sources = grid.array(fill=1.0 / self.DT, z=pygetm.CENTERS)
        self.diffuse(tracer, nuh, cnpar, sources=sources * self.DT * grid.hn)

    @for_each_grid
    def test_relative_source(self, grid):
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        vdif = pygetm.operators.VerticalDiffusion(grid, cnpar=1.0)

        tracer = grid.array(z=pygetm.CENTERS, fill=1.0)

        # Now try without spatial gradient and relative [linear] source term only
        tolerance = 1e-13
        r = 0.1  # relative rate of increase
        # note that sources should be time- and layer-integrated!
        rel_sources = grid.array(fill=r / self.DT, z=pygetm.CENTERS) * self.DT * grid.hn
        for _ in range(self.NSTEP):
            vdif(nuh, self.DT, tracer, ea2=rel_sources)
        expected = 1.0 / (1.0 - r) ** self.NSTEP
        delta = tracer.ma - expected
        rel_delta = delta / expected
        rel_delta_min = rel_delta.min(axis=(1, 2))
        rel_delta_max = rel_delta.max(axis=(1, 2))
        self.assertFalse(
            (rel_delta_min - rel_delta_max).any(),
            "relative error in tracer varies horizontally after using vertical diffusion solver to integrate relative sources: %s vs %s"
            % (rel_delta_min, rel_delta_max),
        )
        rel_error = np.abs(rel_delta).max()
        self.assertLessEqual(
            rel_error,
            tolerance,
            "maximum relative error %s in tracer after using vertical diffusion solver to integrate relative sources exceeds tolerance %s"
            % (rel_error, tolerance),
        )

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_random_layer_height_change(self, cnpar: float, grid: pygetm.domain.Grid):
        tracer = grid.array(z=pygetm.CENTERS, fill_value=1.0)
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=0.01)
        self.diffuse(tracer, nuh, cnpar, change_ho=True)

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_random_diffusivity(self, cnpar: float, grid: pygetm.domain.Grid):
        tracer = grid.array(z=pygetm.CENTERS, fill_value=np.nan)
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=np.nan)
        nuh.fill(10.0 ** rng.uniform(-6.0, 0.0, (nuh.shape[0], 1, 1)))
        tracer.fill(35.0)
        self.diffuse(tracer, nuh, cnpar, tolerance=1e-11)

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_random_diffusivity_random_tracer(
        self, cnpar: float, grid: pygetm.domain.Grid
    ):
        tracer = grid.array(z=pygetm.CENTERS, fill_value=np.nan)
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=np.nan)
        nuh.fill(10.0 ** rng.uniform(-6.0, 0.0, (nuh.shape[0], 1, 1)))
        tracer.fill(rng.uniform(0.0, 1.0, (tracer.shape[0], 1, 1)))
        self.diffuse(tracer, nuh, cnpar, tolerance=1e-11)

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_random_diffusivity_random_tracer_with_sources(
        self, cnpar: float, grid: pygetm.domain.Grid
    ):
        tracer = grid.array(z=pygetm.CENTERS, fill_value=np.nan)
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=np.nan)
        nuh.fill(10.0 ** rng.uniform(-6.0, 0.0, (nuh.shape[0], 1, 1)))
        tracer.fill(rng.uniform(0.0, 1.0, (tracer.shape[0], 1, 1)))
        sources = grid.array(fill=1.0, z=pygetm.CENTERS)
        sources.fill(rng.uniform(0.0, 1.0, (sources.shape[0], 1, 1)))
        self.diffuse(
            tracer, nuh, cnpar, tolerance=1e-11, sources=sources * self.DT * grid.hn
        )

    @repeat
    @for_each_grid
    @for_each_cnpar
    def test_random_diffusivity_random_tracer_with_sources_changing_h(
        self, cnpar: float, grid: pygetm.domain.Grid
    ):
        tracer = grid.array(z=pygetm.CENTERS, fill_value=np.nan)
        nuh = grid.array(z=pygetm.INTERFACES, fill_value=np.nan)
        nuh.fill(10.0 ** rng.uniform(-6.0, 0.0, (nuh.shape[0], 1, 1)))
        tracer.fill(rng.uniform(0.0, 1.0, (tracer.shape[0], 1, 1)))
        sources = grid.array(fill=1.0, z=pygetm.CENTERS)
        sources.fill(rng.uniform(0.0, 1.0, (sources.shape[0], 1, 1)))
        self.diffuse(
            tracer,
            nuh,
            cnpar,
            tolerance=1e-11,
            sources=sources * self.DT * grid.hn,
            change_ho=True,
        )


if __name__ == "__main__":
    unittest.main()

