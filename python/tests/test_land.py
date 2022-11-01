import sys
import os.path
import unittest

import cftime
import numpy as np

import pygetm
import pygetm.util.compare_nc

sys.path.append(os.path.join(os.path.dirname(__file__), "../examples"))

import north_sea


class TestLandMask(unittest.TestCase):
    def setUp(self) -> None:
        setups_dir = "../../../getm-setups"
        if "GETM_SETUPS_DIR" in os.environ:
            setups_dir = os.environ["GETM_SETUPS_DIR"]
        self.setup_dir = os.path.join(setups_dir, "NorthSea")
        self.domain = north_sea.create_domain(
            self.setup_dir, logger=pygetm.parallel.get_logger(level="ERROR")
        )

    def test_nan(self):
        start = cftime.datetime(2006, 1, 2)
        stop = cftime.datetime(2006, 1, 3)

        sim = north_sea.create_simulation(
            self.domain, pygetm.BAROCLINIC, self.setup_dir
        )

        # Set land points (mask==0) of every variable to NaN
        U, V = self.domain.U, self.domain.V
        for array in self.domain.fields.values():
            grid = array.grid
            readonly = not array.all_values.flags.writeable
            if readonly:
                array.all_values.flags.writeable = True
            if (
                array.all_values.dtype == float
                and not array.on_boundary
                and array.ndim > 0
            ):
                array.all_values[..., grid._land] = np.nan
            if readonly:
                if grid in (U, V) and array in (grid.dx, grid.dy):
                    # dx and dy at land-water interface need to be finite
                    # for u,v advection and vertical velocity calculation
                    # The simulation is already protected against that as dx, dy
                    # were readonly. Quietly restore finite values at these locations
                    edges = grid._water_contact & grid._land
                    array.all_values[edges] = array.fill_value
                array.all_values.flags.writeable = False

        sim.start(start, timestep=60.0, split_factor=30, report=60)
        while sim.time < stop:
            sim.advance()
        sim.check_finite()
        sim.finish()

    def test_masked(self):
        skip = "u10", "v10", "t2m", "tcc", "tp", "sp", "zen", "w"
        grid_skip = (
            "lon",
            "lat",
            "dx",
            "dy",
            "area",
            "idx",
            "idy",
            "iarea",
            "rotation",
            "cor",
            "zc",
            "zf",
        )
        stop = cftime.datetime(2006, 1, 3)
        sim = north_sea.create_simulation(
            self.domain, pygetm.BAROCLINIC, self.setup_dir
        )
        north_sea.run(sim, stop=stop)

        for array in self.domain.fields.values():
            skip_this = array.name in skip
            for s in grid_skip:
                if array is getattr(array.grid, s):
                    skip_this = True
            if (
                array.on_boundary
                or array.ndim == 0
                or array.attrs.get("_mask_output", False)
                or skip_this
            ):
                continue
            with self.subTest(name=array.name):
                land_values = array.all_values[..., array.grid.mask.all_values == 0]
                self.assertTrue(np.isfinite(land_values).all())
                self.assertTrue((land_values == array.fill_value).all())


if __name__ == "__main__":
    unittest.main()

