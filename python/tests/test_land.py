import sys
import os
import unittest
from pathlib import Path

import cftime
import numpy as np

import pygetm
from pygetm.constants import CellType

sys.path.append(str(Path(__file__).parent / "../examples"))

import north_sea


RUNTYPES = (
    pygetm.RunType.BAROTROPIC_2D,
    pygetm.RunType.BAROTROPIC_3D,
    pygetm.RunType.BAROCLINIC,
)


def _invalidate(sim: pygetm.Simulation, startup: bool):
    """Set unused points (e.g. land) of every variable to NaN"""
    for array in sim._fields.values():
        if array.all_values.dtype != float:
            continue
        all_values = array.all_values.view()
        all_values.flags.writeable = True
        mask = array.all_mask
        if (
            startup
            and array.all_values.flags.writeable
            and not array.on_boundary
            and array.ndim > 0
        ):
            mask = mask | array.grid.get_mask((CellType.ACTIVE, CellType.BOUNDARY))
        all_values[mask] = np.nan


class TestLandMask(unittest.TestCase):
    def setUp(self) -> None:
        setups_dir = os.environ.get("GETM_SETUPS_DIR", "../../../getm-setups")
        self.setup_dir = Path(setups_dir) / "NorthSea"
        self.domain = north_sea.create_domain(
            self.setup_dir, logger=pygetm.parallel.get_logger(level="ERROR")
        )

    def test_nan(self):
        start = cftime.datetime(2006, 1, 2)
        stop = cftime.datetime(2006, 1, 2, 2)

        for runtype in RUNTYPES:
            with self.subTest(runtype=runtype.name):
                sim = north_sea.create_simulation(self.domain, runtype, self.setup_dir)

                _invalidate(sim, startup=True)
                sim.start(start, timestep=60.0, split_factor=30, report=60)
                _invalidate(sim, startup=False)
                while sim.time < stop:
                    sim.advance(check_finite=True)
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

        for runtype in RUNTYPES:
            with self.subTest(runtype=runtype.name):
                sim = north_sea.create_simulation(self.domain, runtype, self.setup_dir)
                north_sea.run(sim, stop=stop)

                for array in sim._fields.values():
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
                        land_values = array.all_values[array.all_mask]
                        self.assertTrue(np.isfinite(land_values).all())
                        self.assertTrue((land_values == array.fill_value).all())


if __name__ == "__main__":
    unittest.main()
