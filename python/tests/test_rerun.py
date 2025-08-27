import sys
import os.path
import unittest
import datetime

import cftime

import pygetm
import pygetm.util.compare_nc

sys.path.append(os.path.join(os.path.dirname(__file__), "../examples"))

import north_sea


class TestRerun(unittest.TestCase):
    def setUp(self) -> None:
        setups_dir = "../../../getm-setups"
        if "GETM_SETUPS_DIR" in os.environ:
            setups_dir = os.environ["GETM_SETUPS_DIR"]
        self.setup_dir = os.path.join(setups_dir, "NorthSea")
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
