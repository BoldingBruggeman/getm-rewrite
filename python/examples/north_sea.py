#!/usr/bin/env python

import argparse
import datetime
import os.path
import pathlib

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

parser = argparse.ArgumentParser()
parser.add_argument(
    "setup_dir", type=pathlib.Path, help="Path to configuration files", default="."
)
parser.add_argument(
    "--start",
    help="Simulation start time - yyyy-mm-dd hh:mi:ss",
    default="2006-01-02 00:00:00",
)
parser.add_argument(
    "--stop",
    help="Simulation stop time - yyyy-mm-dd hh:mi:ss",
    default="2007-01-01 00:00:00",
)
parser.add_argument(
    "--meteo_dir", type=pathlib.Path, help="Path to ERA5 meteo forcing files"
)
# parser.add_argument('--input_dir', type=pathlib.Path, help='Path to input files', default='input' )
parser.add_argument(
    "--tpxo9_dir", type=pathlib.Path, help="Path to TPXO9 configuration files"
)
parser.add_argument(
    "--tiling", type=argparse.FileType("r"), help="Path to tiling pickle file"
)
parser.add_argument(
    "--initial",
    action="store_true",
    help="Read initial salinity and temerature conditions from file",
)
# parser.add_argument('--no_meteo', action='store_true', help='No meteo forcing')
parser.add_argument(
    "--no_boundaries",
    action="store_false",
    dest="boundaries",
    help="No open boundaries",
)
parser.add_argument(
    "--no_rivers", action="store_false", dest="rivers", help="No river input"
)
parser.add_argument(
    "--no_output",
    action="store_false",
    dest="output",
    help="Do not save any results to NetCDF",
)
parser.add_argument(
    "--debug_output", action="store_true", help="Do not save any results to NetCDF"
)
parser.add_argument("--profile", help="File to save profiling report to")
parser.add_argument(
    "--runtype",
    type=int,
    choices=(pygetm.BAROTROPIC_2D, pygetm.BAROTROPIC_3D, pygetm.BAROCLINIC),
    help="Run type",
    default=pygetm.BAROCLINIC,
)
parser.add_argument("--save_restart", help="File to save restart to")
parser.add_argument("--load_restart", help="File to load restart from")
args = parser.parse_args()

simstart = datetime.datetime.strptime(args.start, "%Y-%m-%d %H:%M:%S")
simstop = datetime.datetime.strptime(args.stop, "%Y-%m-%d %H:%M:%S")
# if args.input_dir is None: args.input_dir = os.path.join(args.setup_dir,'input')
tiling = args.tiling if args.tiling is not None else None
# if args.meteo_dir is None: args.no_meteo = True

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

if args.boundaries:
    kwargs = {"type_2d": -4} if args.tpxo9_dir else {}
    pygetm.legacy.load_bdyinfo(
        domain, os.path.join(args.setup_dir, "bdyinfo.dat"), **kwargs
    )
if args.rivers:
    pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, "riverinfo.dat"))

# Select between original or ERA5 meteo forcing
if args.meteo_dir:
    airsea = pygetm.airsea.FluxesFromMeteo(
        humidity_measure=pygetm.airsea.HumidityMeasure.DEW_POINT_TEMPERATURE,
        calculate_evaporation=True,
    )
else:
    airsea = pygetm.airsea.FluxesFromMeteo(
        humidity_measure=pygetm.airsea.HumidityMeasure.SPECIFIC_HUMIDITY,
        calculate_evaporation=True,
    )

sim = pygetm.Simulation(
    domain,
    runtype=args.runtype,
    advection_scheme=pygetm.AdvectionScheme.HSIMT,
    gotm=os.path.join(args.setup_dir, "gotmturb.nml"),
    airsea=airsea,
    internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
    #        fabm='../../extern/fabm/testcases/fabm-jrc-med_ergom.yaml',
)

sim.logger.info("Setting up TPXO tidal boundary forcing")
if domain.open_boundaries:
    if args.tpxo9_dir is None:
        sim.logger.info("Reading 2D boundary data from file")
        domain.open_boundaries.z.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/2D/bdy.2d.2006.nc"), "elev"
            )
        )
        domain.open_boundaries.u.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/2D/bdy.2d.2006.nc"), "u"
            )
        )
        domain.open_boundaries.v.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/2D/bdy.2d.2006.nc"), "v"
            )
        )
    else:
        sim.logger.info("set handling of TPXO boundary input")
        bdy_lon = domain.open_boundaries.lon
        bdy_lat = domain.open_boundaries.lat
        domain.open_boundaries.z.set(
            pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=args.tpxo9_dir)
        )
        domain.open_boundaries.u.set(
            pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="u", root=args.tpxo9_dir)
        )
        domain.open_boundaries.v.set(
            pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable="v", root=args.tpxo9_dir)
        )

    if sim.runtype == pygetm.BAROCLINIC:
        sim.temp.open_boundaries.type = pygetm.SPONGE
        sim.temp.open_boundaries.values.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/3D/bound_3D.CFSR.2006.nc"), "temp"
            )
        )
        sim.salt.open_boundaries.type = pygetm.SPONGE
        sim.salt.open_boundaries.values.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/3D/bound_3D.CFSR.2006.nc"), "salt"
            )
        )

if domain.rivers:
    for name, river in domain.rivers.items():
        river.flow.set(
            pygetm.input.from_nc(
                os.path.join(args.setup_dir, "Forcing/River/rivers.nc"), name
            )
        )
        if sim.runtype == pygetm.BAROCLINIC:
            river["salt"].set(
                pygetm.input.from_nc(
                    os.path.join(args.setup_dir, "Forcing/River/rivers.nc"),
                    "%s_salt" % name,
                )
            )

if sim.runtype < pygetm.BAROCLINIC:
    sim.sst = sim.airsea.t2m
    if sim.runtype > pygetm.BAROTROPIC_2D:
        sim.turbulence.num[...] = 1e-2
if sim.runtype == pygetm.BAROCLINIC:
    sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)
    sim.temp.set(11.6)
    sim.salt.set(35.2)

if args.meteo_dir:
    sim.logger.info("Setting up ERA5 meteorological forcing")
    ERA_path = os.path.join(args.meteo_dir, "era5_2006.nc")
    sim.airsea.u10.set(pygetm.input.from_nc(ERA_path, "u10"))
    sim.airsea.v10.set(pygetm.input.from_nc(ERA_path, "v10"))
    sim.airsea.t2m.set(pygetm.input.from_nc(ERA_path, "t2m") - 273.15)
    sim.airsea.d2m.set(pygetm.input.from_nc(ERA_path, "d2m") - 273.15)
    sim.airsea.sp.set(pygetm.input.from_nc(ERA_path, "sp"))
    sim.airsea.tcc.set(pygetm.input.from_nc(ERA_path, "tcc"))
    sim.airsea.tp.set(pygetm.input.from_nc(ERA_path, "tp") / 3600.0)
else:
    sim.logger.info("Setting up NS original meteorological forcing")
    met_path = os.path.join(args.setup_dir, "Forcing/Meteo/CFSR.daymean.2006.nc")
    sim.airsea.tcc.set(pygetm.input.from_nc(met_path, "tcc"))
    sim.airsea.t2m.set(pygetm.input.from_nc(met_path, "t2"))
    sim.airsea.qa.set(pygetm.input.from_nc(met_path, "sh"))
    sim.airsea.sp.set(pygetm.input.from_nc(met_path, "slp"))
    sim.airsea.u10.set(pygetm.input.from_nc(met_path, "u10"))
    sim.airsea.v10.set(pygetm.input.from_nc(met_path, "v10"))
    sim.airsea.tp.set(pygetm.input.from_nc(met_path, "precip"))

if sim.fabm:
    sim.logger.info("Setting up FABM dependencies that GETM does not provide")
    if sim.runtype == pygetm.BAROTROPIC_3D:
        sim.fabm.get_dependency("temperature").set(5.0)
        sim.fabm.get_dependency("practical_salinity").set(35.0)

if args.output:
    sim.logger.info("Setting up output")
    output = sim.output_manager.add_netcdf_file(
        "meteo.nc", interval=datetime.timedelta(hours=1), sync_interval=None
    )
    output.request("u10", "v10", "sp", "t2m", "qa", "tcc", "e", "tp", "pe")
    # output.request('qe', 'qh', 'ql', 'swr', 'albedo', 'zen')
    output = sim.output_manager.add_netcdf_file(
        "north_sea_2d.nc", interval=datetime.timedelta(hours=1), sync_interval=None
    )
    output.request("zt", "Dt", "u1", "v1", "tausxu", "tausyv")
    if args.debug_output:
        output.request("maskt", "masku", "maskv")
        output.request("U", "V")
        output.request("Du", "Dv", "dpdx", "dpdy", "z0bu", "z0bv", "z0bt")  # , 'u_taus'
        output.request("ru", "rru", "rv", "rrv")
    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file(
            "north_sea_3d.nc", interval=datetime.timedelta(hours=6), sync_interval=None
        )
        output.request("uk", "vk", "ww", "SS", "num")
        if args.debug_output:
            output.request("fpk", "fqk", "advpk", "advqk")  # 'diffpk', 'diffqk')
    if sim.runtype == pygetm.BAROCLINIC:
        output.request("temp", "salt", "rho", "NN", "rad", "sst", "hnt", "nuh")
        if args.debug_output:
            output.request("idpdx", "idpdy")
        if sim.fabm:
            output.request("par", "med_ergom_o2", "med_ergom_OFL", "med_ergom_dd")

if args.save_restart:
    sim.output_manager.add_restart(args.save_restart)

if args.load_restart:
    simstart = sim.load_restart(args.load_restart)

sim.start(simstart, timestep=60.0, split_factor=30, report=60, profile=args.profile)
while sim.time < simstop:
    sim.advance()
sim.finish()
