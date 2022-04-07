import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

getm_setups_dirs = ('/data/kb/getm-setups.SF', '../../../getm-setups')
getm_setups_dir = next(filter(os.path.isdir, getm_setups_dirs))

igotm_data_dirs = ('/server/data', '../../../igotm/data')
igotm_data_dir = next(filter(os.path.isdir, igotm_data_dirs))

domain = pygetm.legacy.domain_from_topo(os.path.join(getm_setups_dir, 'NorthSea/Topo/NS6nm.v01.nc'), nlev=30, z0_const=0.001)
pygetm.legacy.load_bdyinfo(domain, os.path.join(getm_setups_dir, 'NorthSea/bdyinfo.dat'))
pygetm.legacy.load_riverinfo(domain, os.path.join(getm_setups_dir, 'NorthSea/riverinfo.dat'))
sim = pygetm.Simulation(domain, runtype=pygetm.BAROCLINIC, advection_scheme=pygetm.HSIMT,
    gotm=os.path.join(getm_setups_dir, 'NorthSea/gotmturb.nml'),
#    fabm='../../extern/fabm/testcases/fabm-jrc-med_ergom.yaml',
)

#sim.input_manager.debug_nc_reads()

if True:
    sim.logger.info('Setting up output')
    output = sim.output_manager.add_netcdf_file('north_sea.nc', interval=60, sync_interval=200000)
    output.request(('u10', 'v10', 'sp', 't2m', 'd2m', 'tcc'))
    #output.request(('qe', 'qh', 'ql', 'swr', 'albedo', 'zen'))
    output.request(('U', 'V'), mask=True)
    output.request(('zt', 'Dt', 'Du', 'Dv', 'masku', 'maskv'))
    output.request(('dpdx', 'dpdy', 'tausxu', 'tausyv', 'z0bu', 'z0bv', 'z0bt'))   #, 'u_taus'
    output.request(('ru', 'rru', 'rv', 'rrv'))
    if sim.runtype > pygetm.BAROTROPIC_2D:
        output.request(('uk', 'vk', 'ww', 'SS', 'fpk', 'fqk', 'advpk', 'advqk', 'nuh',))
    if False and sim.do_fabm:
        output.request(('med_ergom_o2', 'med_ergom_OFL', 'med_ergom_dd'))
    if sim.runtype == pygetm.BAROCLINIC:
        output.request(('temp', 'salt', 'rho', 'NN', 'sst', 'hnt', 'rad', 'par'))

sim.logger.info('Setting up ERA meteorological forcing')
met_path = os.path.join(getm_setups_dir, 'NorthSea/Forcing/Meteo/CFSR.daymean.2006.nc')
sim.airsea.tcc.set(pygetm.input.from_nc(met_path, 'tcc'))
sim.airsea.t2m.set(pygetm.input.from_nc(met_path, 't2'))
sim.airsea.d2m.set(pygetm.input.from_nc(met_path, 'sh'))
sim.airsea.sp.set(pygetm.input.from_nc(met_path, 'slp'))
sim.airsea.u10.set(pygetm.input.from_nc(met_path, 'u10'))
sim.airsea.v10.set(pygetm.input.from_nc(met_path, 'v10'))
if sim.runtype < pygetm.BAROCLINIC:
    sim.sst = sim.airsea.t2m
    sim.turbulence.num[...]=1e-2

sim.logger.info('Setting up TPXO tidal boundary forcing')
tpxo_dir = os.path.join(igotm_data_dir, 'TPXO9')
bdy_lon = domain.T.lon.all_values[domain.bdy_j, domain.bdy_i]
bdy_lat = domain.T.lat.all_values[domain.bdy_j, domain.bdy_i]
if domain.open_boundaries:
    sim.zbdy.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/2D/bdy.2d.2006.nc'), 'elev'), on_grid=True)
    sim.bdyu.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/2D/bdy.2d.2006.nc'), 'u'), on_grid=True)
    sim.bdyv.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/2D/bdy.2d.2006.nc'), 'v'), on_grid=True)
    if sim.runtype == pygetm.BAROCLINIC:
        sim.temp.boundaries.type = pygetm.SPONGE
        sim.temp.boundaries.values.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/3D/bound_3D.CFSR.2006.nc'), 'temp'), on_grid=True)
        sim.salt.boundaries.type = pygetm.SPONGE
        sim.salt.boundaries.values.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/3D/bound_3D.CFSR.2006.nc'), 'salt'), on_grid=True)

for name, river in domain.rivers.items():
    river.flow.set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/River/rivers.nc'), name))
    river['salt'].set(pygetm.input.from_nc(os.path.join(getm_setups_dir, 'NorthSea/Forcing/River/rivers.nc'), '%s_salt' % name))

if sim.fabm_model:
    sim.logger.info('Setting up FABM dependencies that GETM does not provide')
    sim.get_fabm_dependency('bottom_stress').set(0)
    if sim.runtype == pygetm.BAROTROPIC_3D:
        sim.get_fabm_dependency('temperature').set(5.)
        sim.get_fabm_dependency('practical_salinity').set(35.)

#sim.start(datetime.datetime(2006, 2, 1), timestep=60., split_factor=30, profile='north_sea') #, profile='north_sea')
sim.start(datetime.datetime(2006, 2, 1), timestep=60., split_factor=30) #, profile='north_sea')
while sim.time < datetime.datetime(2006, 2, 2):
    sim.advance()

sim.finish()
