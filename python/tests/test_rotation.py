import sys
import os.path
import unittest

import pygetm
import pygetm.util.compare_nc

sys.path.append(os.path.join(os.path.dirname(__file__), "../examples"))

import north_sea


class TestRotation(unittest.TestCase):
    def setUp(self) -> None:
        self.setup_dir = "../../../getm-setups/NorthSea"
        self.domain = north_sea.create_domain(
            self.setup_dir,
            use_boundaries=False,
            use_rivers=False,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

    def test_2d(self):
        outputs = "u1", "v1", "zt", "U", "V"
        stop = "2006-01-03 00:00:00"
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
        sim = north_sea.create_simulation(domain, pygetm.BAROTROPIC_2D, self.setup_dir)
        output = sim.output_manager.add_netcdf_file("result_ref.nc", interval=-1)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        sim = north_sea.create_simulation(
            domain_rot, pygetm.BAROTROPIC_2D, self.setup_dir
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
        stop = "2006-01-02 00:30:00"
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
        sim = north_sea.create_simulation(domain, pygetm.BAROCLINIC, self.setup_dir)
        output = sim.output_manager.add_netcdf_file("result_ref.nc", interval=30)
        output.request(*outputs)
        north_sea.run(sim, stop=stop)

        sim = north_sea.create_simulation(domain_rot, pygetm.BAROCLINIC, self.setup_dir)
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

