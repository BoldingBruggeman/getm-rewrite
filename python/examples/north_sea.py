import datetime
import os.path

import pygetm
import pygetm.legacy
import pygetm.input.tpxo

getm_setups_dir = '../../../getm-setups'
igotm_data_dirs = ('/server/data', '../../../igotm/data')

igotm_data_dir = next(filter(os.path.isdir, igotm_data_dirs))
era_path = os.path.join(igotm_data_dir, 'ERA-interim/2016.nc')

domain = pygetm.legacy.domain_from_topo(os.path.join(getm_setups_dir, 'NorthSea/Topo/NS6nm.v01.nc'), nlev=30, z0_const=0.001)
pygetm.legacy.load_bdyinfo(domain, os.path.join(getm_setups_dir, 'NorthSea/bdyinfo.dat'))
sim = pygetm.Simulation(domain, runtype=4, advection_scheme=1, fabm='../../extern/fabm/testcases/fabm-jrc-med_ergom.yaml')

#sim.input_manager.debug_nc_reads()

sim.logger.info('Setting up output')
output = sim.output_manager.add_netcdf_file('north_sea.nc', interval=60)
output.request(('u10', 'v10', 'sp', 't2m', 'd2m', 'tcc'))
output.request(('U', 'V'))
output.request('zt')
output.request(('dpdx', 'dpdy', 'tausxu', 'tausyv'))
output.request(('med_ergom_o2', 'med_ergom_OFL'))

sim.logger.info('Setting up ERA meteorological forcing')
era_kwargs = {'preprocess': lambda ds: ds.isel(time=slice(4, -4))}
sim.airsea.tcc.set(pygetm.input.from_nc(era_path, 'tcc', **era_kwargs))
sim.airsea.t2m.set(pygetm.input.from_nc(era_path, 't2m', **era_kwargs) - 273.15)
sim.airsea.d2m.set(pygetm.input.from_nc(era_path, 'd2m', **era_kwargs) - 273.15)
sim.airsea.sp.set(pygetm.input.from_nc(era_path, 'sp', **era_kwargs))
sim.airsea.u10.set(pygetm.input.from_nc(era_path, 'u10', **era_kwargs))
sim.airsea.v10.set(pygetm.input.from_nc(era_path, 'v10', **era_kwargs))

sim.logger.info('Setting up TPXO tidal boundary forcing')
tpxo_dir = os.path.join(igotm_data_dir, 'TPXO9')
bdy_lon = domain.T.lon[domain.bdy_j, domain.bdy_i]
bdy_lat = domain.T.lat[domain.bdy_j, domain.bdy_i]
sim.zbdy.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, root=tpxo_dir), on_grid=True)
sim.bdyu.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='u', root=tpxo_dir), on_grid=True)
sim.bdyv.set(pygetm.input.tpxo.get(bdy_lon, bdy_lat, variable='v', root=tpxo_dir), on_grid=True)

sim.logger.info('Setting up FABM dependencies that GETM does not provide')
sim.get_fabm_dependency('downwelling_photosynthetic_radiative_flux').set(0)
sim.get_fabm_dependency('surface_downwelling_photosynthetic_radiative_flux').set(0)
sim.get_fabm_dependency('bottom_stress').set(0)

sim.start(datetime.datetime(2016, 1, 1), timestep=60., split_factor=30)

while sim.time < datetime.datetime(2016, 1, 5):
    sim.advance()

sim.finish()
