#!/usr/bin/env python

import argparse
import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

parser = argparse.ArgumentParser()
parser.add_argument('setup_dir', help='Path to configuration files', default='.')
args = parser.parse_args()

domain = pygetm.legacy.domain_from_topo(os.path.join(args.setup_dir, 'Topo/NS6nm.v01.nc'), nlev=30, z0_const=0.001)

pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, 'bdyinfo.dat'))
pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, 'riverinfo.dat'))

sim = pygetm.Simulation(domain,
        runtype=pygetm.BAROCLINIC,
        advection_scheme=pygetm.AdvectionScheme.HSIMT,
        gotm=os.path.join(args.setup_dir, 'gotmturb.nml'),
        airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.SPECIFIC_HUMIDITY),
        internal_pressure_method=pygetm.InternalPressure.BLUMBERG_MELLOR,
)

sim.logger.info('Reading 2D boundary data from file')
domain.open_boundaries.z.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'elev'))
domain.open_boundaries.u.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'u'))
domain.open_boundaries.v.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'v'))

sim.temp.open_boundaries.type = pygetm.SPONGE
sim.temp.open_boundaries.values.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/3D/bound_3D.CFSR.2006.nc'), 'temp'))
sim.salt.open_boundaries.type = pygetm.SPONGE
sim.salt.open_boundaries.values.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/3D/bound_3D.CFSR.2006.nc'), 'salt'))

for name, river in domain.rivers.items():
    river.flow.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/River/rivers.nc'), name))
    river['salt'].set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/River/rivers.nc'), '%s_salt' % name))

sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)
sim.temp.set(11.6)
sim.salt.set(35.2)

sim.logger.info('Setting up NS original meteorological forcing')
met_path = os.path.join(args.setup_dir, 'Forcing/Meteo/CFSR.daymean.2006.nc')
sim.airsea.tcc.set(pygetm.input.from_nc(met_path, 'tcc'))
sim.airsea.t2m.set(pygetm.input.from_nc(met_path, 't2'))
sim.airsea.qa.set(pygetm.input.from_nc(met_path, 'sh'))
sim.airsea.sp.set(pygetm.input.from_nc(met_path, 'slp'))
sim.airsea.u10.set(pygetm.input.from_nc(met_path, 'u10'))
sim.airsea.v10.set(pygetm.input.from_nc(met_path, 'v10'))

sim.logger.info('Setting up output')
output = sim.output_manager.add_netcdf_file('meteo.nc', interval=datetime.timedelta(hours=1), sync_interval=None)
output.request(('u10', 'v10', 'sp', 't2m', 'qa', 'tcc'))
output = sim.output_manager.add_netcdf_file('north_sea_2d.nc', interval=datetime.timedelta(hours=1), sync_interval=None)
output.request(('zt', 'Dt', 'u1', 'v1', 'tausxu', 'tausyv', ))
output = sim.output_manager.add_netcdf_file('north_sea_3d.nc', interval=datetime.timedelta(hours=6), sync_interval=None)
output.request(('uk', 'vk', 'ww', 'SS', 'num',))
output.request(('temp', 'salt', 'rho', 'NN', 'rad', 'sst', 'hnt', 'nuh',))

sim.start(datetime.datetime(2006, 1, 2), timestep=60., split_factor=30)
while sim.time < datetime.datetime(2007, 1, 1):
    sim.advance()
sim.finish()