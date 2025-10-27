#!/usr/bin/env python

import argparse
import datetime
from pathlib import Path

import pygetm
import pygetm.legacy

parser = argparse.ArgumentParser()
parser.add_argument(
    "setup_dir",
    type=Path,
    default=".",
    help="Path to configuration files (sylt directory from https://sourceforge.net/p/getm/getm-setups)",
)
args = parser.parse_args()

domain = pygetm.legacy.domain_from_topo(
    args.setup_dir / "topo.nc", z0=0.01, f=pygetm.domain.coriolis(50.0)
)
pygetm.legacy.load_bdyinfo(domain, args.setup_dir / "bdyinfo.dat")

sim = pygetm.Simulation(
    domain,
    runtype=pygetm.BAROTROPIC_3D,
    vertical_coordinates=pygetm.vertical_coordinates.Sigma(10),
    airsea=pygetm.airsea.Fluxes(),
    gotm=args.setup_dir / "gotmturb.nml",
    Dcrit=0.2,
    Dmin=0.05,
)

sim.logger.info("Reading 2D boundary data from file")
bdy_2d_path = args.setup_dir / "sylt_bdy.nc"
sim.open_boundaries.z.set(pygetm.input.from_nc(bdy_2d_path, "elev"))

sim.logger.info("Setting up output")
output = sim.output_manager.add_netcdf_file("sylt_2d.nc", interval=100)
output.request("zt", "u1", "v1", grid=sim.T)
output = sim.output_manager.add_netcdf_file("sylt_3d.nc", interval=500)
output.request("uk", "vk", "tke", "num", "nuh", "eps", grid=sim.T)

sim.start(
    datetime.datetime(2000, 1, 1),
    timestep=4.4714,
    split_factor=10,
    report=60,
    report_totals=600,
    profile="sylt",
)
for _ in range(30000):
    sim.advance()
sim.finish()
