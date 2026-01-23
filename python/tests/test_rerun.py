import sys
import os
import unittest
import datetime
from pathlib import Path

import cftime

import pygetm

sys.path.append(str(Path(__file__).parent / "../examples"))

import north_sea


class TestRerun(unittest.TestCase):
    def setUp(self) -> None:
        setups_dir = os.environ.get("GETM_SETUPS_DIR", "../../../getm-setups")
        self.setup_dir = Path(setups_dir) / "NorthSea"
        self.domain = north_sea.create_domain(
            self.setup_dir,
            use_boundaries=True,
            use_rivers=True,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

    def test_2d(self):
        stop = cftime.datetime(2006, 1, 3)
        sim = north_sea.create_simulation(
            self.domain, pygetm.RunType.BAROTROPIC_2D, self.setup_dir
        )
        output = sim.output_manager.add_netcdf_file(
            "north_sea_2d.nc", interval=datetime.timedelta(hours=1), sync_interval=None
        )
        output.request("zt", "u1", "v1")

        north_sea.run(sim, stop=stop)
        self.assertFalse(sim.output_manager._active_files)
        self.assertFalse(sim.output_manager._startable_files)
        self.assertFalse(sim.output_manager._stoppable_files)

        output = sim.output_manager.add_netcdf_file(
            "north_sea_2d_.nc", interval=datetime.timedelta(hours=1), sync_interval=None
        )
        output.request("zt", "u1", "v1")

        north_sea.run(sim, stop=stop)


if __name__ == "__main__":
    unittest.main()
