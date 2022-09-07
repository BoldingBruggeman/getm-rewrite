#!/usr/bin/env python

import argparse
import datetime
import os.path

import pygetm
import pygetm.legacy

parser = argparse.ArgumentParser()
parser.add_argument("setup_dir", help="Path to configuration files", default=".")
args = parser.parse_args()

domain = pygetm.legacy.domain_from_topo(
    os.path.join(args.setup_dir, "Topo/NS6nm.v01.nc"),
    nlev=30,
    vertical_coordinate_method=pygetm.VerticalCoordinates.GVC,
    Dgamma=40.0,
    ddu=0.75,
    ddl=0.5,
    Dcrit=0.2,
    Dmin=0.05,
    z0=0.001,
)
pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, "bdyinfo.dat"))
pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, "riverinfo.dat"))

sim = pygetm.Simulation(
    domain,
    runtype=pygetm.BAROCLINIC,
    advection_scheme=pygetm.AdvectionScheme.HSIMT,
    gotm=os.path.join(args.setup_dir, "gotmturb.nml"),
    airsea=pygetm.airsea.FluxesFromMeteo(
        humidity_measure=pygetm.HumidityMeasure.SPECIFIC_HUMIDITY,
        calculate_evaporation=True,
    ),
    internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
)

sim.logger.info("Reading 2D boundary data from file")
bdy_2d_path = os.path.join(args.setup_dir, "Forcing/2D/bdy.2d.2006.nc")
domain.open_boundaries.z.set(pygetm.input.from_nc(bdy_2d_path, "elev"))
domain.open_boundaries.u.set(pygetm.input.from_nc(bdy_2d_path, "u"))
domain.open_boundaries.v.set(pygetm.input.from_nc(bdy_2d_path, "v"))

bdy_3d_path = os.path.join(args.setup_dir, "Forcing/3D/bound_3D.CFSR.2006.nc")
sim.temp.open_boundaries.type = pygetm.SPONGE
sim.temp.open_boundaries.values.set(pygetm.input.from_nc(bdy_3d_path, "temp"))
sim.salt.open_boundaries.type = pygetm.SPONGE
sim.salt.open_boundaries.values.set(pygetm.input.from_nc(bdy_3d_path, "salt"))

river_path = os.path.join(args.setup_dir, "Forcing/River/rivers.nc")
for name, river in domain.rivers.items():
    river.flow.set(pygetm.input.from_nc(river_path, name))
    river["salt"].set(pygetm.input.from_nc(river_path, "%s_salt" % name))

sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)
sim.temp.set(11.6)
sim.salt.set(35.2)
sim.density.convert_ts(sim.salt, sim.temp)

sim.logger.info("Setting up NS original meteorological forcing")
met_path = os.path.join(args.setup_dir, "Forcing/Meteo/CFSR.daymean.2006.nc")
sim.airsea.tcc.set(pygetm.input.from_nc(met_path, "tcc"))
sim.airsea.t2m.set(pygetm.input.from_nc(met_path, "t2"))
sim.airsea.qa.set(pygetm.input.from_nc(met_path, "sh"))
sim.airsea.sp.set(pygetm.input.from_nc(met_path, "slp"))
sim.airsea.u10.set(pygetm.input.from_nc(met_path, "u10"))
sim.airsea.v10.set(pygetm.input.from_nc(met_path, "v10"))
sim.airsea.tp.set(pygetm.input.from_nc(met_path, "precip"))

sim.logger.info("Setting up output")
output = sim.output_manager.add_netcdf_file(
    "meteo.nc", interval=datetime.timedelta(hours=1)
)
output.request("u10", "v10", "sp", "t2m", "qa", "tcc")
output = sim.output_manager.add_netcdf_file(
    "north_sea_2d.nc", interval=datetime.timedelta(hours=1)
)
output.request("zt", "Dt", "u1", "v1", "tausxu", "tausyv")
output = sim.output_manager.add_netcdf_file(
    "north_sea_3d.nc", interval=datetime.timedelta(hours=6)
)
output.request("uk", "vk", "ww", "SS", "num")
output.request("temp", "salt", "rho", "NN", "rad", "sst", "hnt", "nuh")

sim.start(datetime.datetime(2006, 1, 2), timestep=60.0, split_factor=30, report=60)
while sim.time < datetime.datetime(2007, 1, 1):
    sim.advance()
sim.finish()
