import sys
import os
import unittest
from pathlib import Path

import cftime

import pygetm
import pygetm.util.compare_nc

sys.path.append(str(Path(__file__).parent / "../examples"))

import north_sea


class TestRotation(unittest.TestCase):
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
        outputs = ("u1", "v1", "zt", "U", "V", "u_bdy", "v_bdy", "zt_bdy")
        stop = cftime.datetime(2006, 1, 3)
        name_map = {
            "u1": "v1",
            "v1": "u1",
            "U": "V",
            "V": "U",
            "latu": "latv",
            "latv": "latu",
            "lonu": "lonv",
            "lonv": "lonu",
        }
        flip = ("v1", "V")

        domain = self.domain
        domain_rot = domain.rotate()
        self.assertTrue((domain.lon.T[::-1, :] == domain_rot.lon)[1:-1, 1:-1].all())
        self.assertTrue((domain.lat.T[::-1, :] == domain_rot.lat)[1:-1, 1:-1].all())
        self.assertTrue((domain.dx.T[::-1, :] == domain_rot.dy)[1:-1, 1:-1].all())
        self.assertTrue((domain.dy.T[::-1, :] == domain_rot.dx)[1:-1, 1:-1].all())
        self.assertTrue((domain.area.T[::-1, :] == domain_rot.area)[1:-1, 1:-1].all())
        sim = north_sea.create_simulation(
            domain, pygetm.RunType.BAROTROPIC_2D, self.setup_dir
        )

        output = sim.output_manager.add_netcdf_file("result_ref.nc", interval=-1)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        sim = north_sea.create_simulation(
            domain_rot, pygetm.RunType.BAROTROPIC_2D, self.setup_dir
        )
        sim.momentum._ufirst = not sim.momentum._ufirst
        sim.momentum._u3dfirst = not sim.momentum._u3dfirst
        sim.momentum.uadv.ufirst = not sim.momentum.uadv.ufirst
        sim.momentum.vadv.ufirst = not sim.momentum.vadv.ufirst
        sim.tracers._advection.ufirst = not sim.tracers._advection.ufirst
        output = sim.output_manager.add_netcdf_file("result_rot.nc", interval=-1)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        match = pygetm.util.compare_nc.compare(
            "result_rot.nc",
            "result_ref.nc",
            rotate=True,
            name_map=name_map,
            flip=flip,
            tolerance=1e-14,
        )
        self.assertTrue(match)

    def test_3d(self):
        outputs = "uk", "vk", "temp", "salt", "zt", "nuh"
        stop = cftime.datetime(2006, 1, 2, 0, 30)
        name_map = {
            "uk": "vk",
            "vk": "uk",
            "latu": "latv",
            "latv": "latu",
            "lonu": "lonv",
            "lonv": "lonu",
            "zcu": "zcv",
            "zcv": "zcu",
        }
        flip = ("vk",)

        domain = self.domain
        domain_rot = domain.rotate()
        sim = north_sea.create_simulation(
            domain, pygetm.RunType.BAROCLINIC, self.setup_dir
        )
        output = sim.output_manager.add_netcdf_file("result_ref.nc", interval=30)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        sim = north_sea.create_simulation(
            domain_rot, pygetm.RunType.BAROCLINIC, self.setup_dir
        )
        sim.momentum._ufirst = not sim.momentum._ufirst
        sim.momentum._u3dfirst = not sim.momentum._u3dfirst
        sim.momentum.uadv.ufirst = not sim.momentum.uadv.ufirst
        sim.momentum.vadv.ufirst = not sim.momentum.vadv.ufirst
        sim.tracers._advection.ufirst = not sim.tracers._advection.ufirst
        output = sim.output_manager.add_netcdf_file("result_rot.nc", interval=30)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        match = pygetm.util.compare_nc.compare(
            "result_rot.nc",
            "result_ref.nc",
            rotate=True,
            name_map=name_map,
            flip=flip,
            tolerance=1e-12,
        )
        self.assertTrue(match)


if __name__ == "__main__":
    unittest.main()
