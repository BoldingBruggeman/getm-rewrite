import datetime
import os.path
import argparse

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

parser = argparse.ArgumentParser()
parser.add_argument('setup_dir', nargs='?', help='Path to GETM Norht Sea setup')
parser.add_argument('--no-output', action='store_false', dest='output', help='Path to GETM Norht Sea setup')
args = parser.parse_args()

if args.setup_dir is None:
    getm_setups_dirs = ('/data/kb/getm-setups.SF', '../../../getm-setups')
    getm_setups_dir = next(filter(os.path.isdir, getm_setups_dirs))
    args.setup_dir = os.path.join(getm_setups_dir, 'NorthSea')

domain = pygetm.legacy.domain_from_topo(os.path.join(args.setup_dir, 'Topo/NS6nm.v01.nc'), nlev=30, z0_const=0.001)
pygetm.legacy.load_bdyinfo(domain, os.path.join(args.setup_dir, 'bdyinfo.dat'))
pygetm.legacy.load_riverinfo(domain, os.path.join(args.setup_dir, 'riverinfo.dat'))

sim = pygetm.Simulation(domain, runtype=pygetm.BAROCLINIC, advection_scheme=pygetm.AdvectionScheme.HSIMT,
    gotm=os.path.join(args.setup_dir, 'gotmturb.nml'),
    airsea=pygetm.airsea.FluxesFromMeteo(humidity_measure=pygetm.airsea.HumidityMeasure.SPECIFIC_HUMIDITY),
#    fabm='../../extern/fabm/testcases/fabm-jrc-med_ergom.yaml',
)
sim.radiation.set_jerlov_type(pygetm.radiation.JERLOV_II)

domain.plot().savefig('domain-%02i.png' % domain.tiling.rank)
#sim.input_manager.debug_nc_reads()

if args.output:
    debug_output = False
    sim.logger.info('Setting up output')
    output = sim.output_manager.add_netcdf_file('meteo.nc', interval=60, sync_interval=200000)
    output.request(('u10', 'v10', 'sp', 't2m', 'qa', 'tcc'))
    #output.request(('qe', 'qh', 'ql', 'swr', 'albedo', 'zen'))
    output = sim.output_manager.add_netcdf_file('north_sea_2d.nc', interval=60, sync_interval=200000)
    output.request(('U', 'V'), mask=True)
    output.request(('zt', 'Dt', 'tausxu', 'tausyv', ))
    if debug_output:
        output.request(('maskt', 'masku', 'maskv', ))   #, 'u_taus'
        output.request(('Du', 'Dv', 'dpdx', 'dpdy', 'z0bu', 'z0bv', 'z0bt'))   #, 'u_taus'
        output.request(('ru', 'rru', 'rv', 'rrv'))
    if sim.runtype > pygetm.BAROTROPIC_2D:
        output = sim.output_manager.add_netcdf_file('north_sea_3d.nc', interval=360, sync_interval=200000)
        output.request(('uk', 'vk', 'ww', 'SS', 'num',))
        if debug_output:
            output.request(('fpk', 'fqk', 'advpk', 'advqk',))
    if sim.runtype == pygetm.BAROCLINIC:
        output.request(('temp', 'salt', 'rho', 'NN', 'rad', 'sst', 'hnt', 'nuh',))
        if debug_output:
            output.request(( ))
        if sim.fabm_model:
            output.request(('par', 'med_ergom_o2', 'med_ergom_OFL', 'med_ergom_dd'))

sim.logger.info('Setting up ERA meteorological forcing')
met_path = os.path.join(args.setup_dir, 'Forcing/Meteo/CFSR.daymean.2006.nc')
sim.airsea.tcc.set(pygetm.input.from_nc(met_path, 'tcc'))
sim.airsea.t2m.set(pygetm.input.from_nc(met_path, 't2'))
sim.airsea.qa.set(pygetm.input.from_nc(met_path, 'sh'))
sim.airsea.sp.set(pygetm.input.from_nc(met_path, 'slp'))
sim.airsea.u10.set(pygetm.input.from_nc(met_path, 'u10'))
sim.airsea.v10.set(pygetm.input.from_nc(met_path, 'v10'))
if sim.runtype < pygetm.BAROCLINIC:
    sim.sst = sim.airsea.t2m
    sim.turbulence.num[...]=1e-2
if sim.runtype == pygetm.BAROCLINIC:
    sim.temp.set(10.)

sim.logger.info('Setting up TPXO tidal boundary forcing')
bdy_lon = domain.T.lon.all_values[domain.bdy_j, domain.bdy_i]
bdy_lat = domain.T.lat.all_values[domain.bdy_j, domain.bdy_i]
if domain.open_boundaries:
    sim.zbdy.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'elev'), on_grid=True)
    sim.bdyu.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'u'), on_grid=True)
    sim.bdyv.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/2D/bdy.2d.2006.nc'), 'v'), on_grid=True)
    if sim.runtype == pygetm.BAROCLINIC:
        sim.temp.boundaries.type = pygetm.SPONGE
        sim.temp.boundaries.values.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/3D/bound_3D.CFSR.2006.nc'), 'temp'), on_grid=True)
        sim.salt.boundaries.type = pygetm.SPONGE
        sim.salt.boundaries.values.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/3D/bound_3D.CFSR.2006.nc'), 'salt'), on_grid=True)

for name, river in domain.rivers.items():
    river.flow.set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/River/rivers.nc'), name))
    river['salt'].set(pygetm.input.from_nc(os.path.join(args.setup_dir, 'Forcing/River/rivers.nc'), '%s_salt' % name))

if sim.fabm_model:
    sim.logger.info('Setting up FABM dependencies that GETM does not provide')
    sim.get_fabm_dependency('bottom_stress').set(0)
    if sim.runtype == pygetm.BAROTROPIC_3D:
        sim.get_fabm_dependency('temperature').set(5.)
        sim.get_fabm_dependency('practical_salinity').set(35.)

sim.start(datetime.datetime(2006, 1, 2), timestep=60., split_factor=30) #, profile='north_sea')
while sim.time < datetime.datetime(2006, 1, 3):
    sim.advance()

sim.finish()
