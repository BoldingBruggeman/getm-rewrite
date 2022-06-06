#!/usr/bin/env python

import argparse
import datetime
import os.path
import pathlib

import numpy as np

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

parser = argparse.ArgumentParser()
parser.add_argument('setup_dir', type=pathlib.Path, help='Path to configuration files', default='.')
parser.add_argument('--start', help='Simulation start time - yyyy-mm-dd hh:mi:ss', default='2004-01-01 00:00:00')
parser.add_argument('--stop', help='Simulation stop time - yyyy-mm-dd hh:mi:ss', default='2004-01-11 00:00:00')
parser.add_argument('--nx', type=int, help='Number of tracer points in x direction', default=100)
parser.add_argument('--ny', type=int, help='Number of tracer points in y direction', default=30)
parser.add_argument('--profile', help='File to save profiling report to')
parser.add_argument('--no_output', action='store_false', dest='output', help='Do not save any results to NetCDF')
args = parser.parse_args()

simstart = datetime.datetime.strptime(args.start, '%Y-%m-%d %H:%M:%S')
simstop = datetime.datetime.strptime(args.stop, '%Y-%m-%d %H:%M:%S')

domain = pygetm.domain.create_spherical(lon=np.linspace(0., 1.25, args.nx + 1), lat=np.linspace(45., 45.25, args.ny + 1), interfaces=True, nz=30,
    z0=0.01, f=pygetm.domain.coriolis(45.),
    vertical_coordinate_method=pygetm.VerticalCoordinates.GVC, Dgamma=40., ddu=1.0, ddl=1.0, Dcrit=0.1, Dmin=0.02
)

ry, rx = np.meshgrid(np.linspace(-0.5, 0.5, args.nx), np.linspace(-0.5, 0.5, args.ny))
H = 100 * np.exp(-0.5 * 10 * (rx**2 + ry**2))
domain.set_bathymetry(H)

sim = pygetm.Simulation(domain,
        runtype=pygetm.BAROCLINIC,
        advection_scheme=pygetm.AdvectionScheme.HSIMT,
        gotm=os.path.join(args.setup_dir, 'gotmturb.nml'),
        airsea=pygetm.airsea.Fluxes(),
        internal_pressure_method=pygetm.InternalPressure.SHCHEPETKIN_MCWILLIAMS,
)

sim.airsea.taux[...] = 0.1
sim.airsea.tauy[...] = 0.
sim.airsea.sp[...] = 0.
sim.airsea.shf[...] = 0.
sim.airsea.swr[...] = 0.
sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_I)

sim.salt.set(0.)

dat = np.loadtxt(os.path.join(args.setup_dir, 'tprof.dat'), skiprows=1)
z, t = dat[::-1, 0], dat[::-1, 1]
sim.temp.values[...] = np.interp(domain.T.zc.values.ravel(), z, t).reshape(sim.temp.shape)

if args.output:
    sim.logger.info('Setting up output')
    basename = os.path.splitext(os.path.basename(__file__))[0]
    output = sim.output_manager.add_netcdf_file('%s_2d.nc' % basename, interval=datetime.timedelta(hours=1), sync_interval=None)
    output.request(('zt', 'Dt', 'u1', 'v1', 'tausxu', 'tausyv', ))
    output = sim.output_manager.add_netcdf_file('%s_3d.nc' % basename, interval=datetime.timedelta(hours=1), sync_interval=None)
    output.request(('uk', 'vk', 'ww', 'SS', 'num',))
    output.request(('temp', 'salt', 'rho', 'NN', 'sst', 'hnt', 'nuh',))

sim.start(simstart, timestep=10., split_factor=60, report=360, profile=args.profile)
while sim.time < simstop:
    sim.advance()
sim.finish()
