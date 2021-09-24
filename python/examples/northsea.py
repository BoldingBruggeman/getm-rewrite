import datetime
import os.path

import numpy

import pygetm
import pygetm.domain
import pygetm.input
import pygetm.legacy
import pygetm.input.tpxo

getm_setups_dir = '../../../getm-setups'

domain = pygetm.legacy.domain_from_topo(os.path.join(getm_setups_dir, 'NorthSea/Topo/NS6nm.v01.nc'), nlev=1, z0_const=0.001)
pygetm.legacy.load_bdyinfo(domain, os.path.join(getm_setups_dir, 'NorthSea/bdyinfo.dat'))
sim = pygetm.Simulation(domain, runtype=1, advection_scheme=1)

bdy_lon = domain.T.lon[domain.bdy_j, domain.bdy_i]
bdy_lat = domain.T.lat[domain.bdy_j, domain.bdy_i]
tpxo = pygetm.input.tpxo.Dataset(bdy_lon, bdy_lat, verbose=True, root='../../../igotm/data/TPXO9')

# No surface forcing
tausx = domain.U.array(fill=0.0)
tausy = domain.V.array(fill=0.0)
sp = domain.T.array(fill=0.0)

timestep = 60.
timedelta = datetime.timedelta(seconds=timestep)
date = datetime.datetime(2000, 1, 1)
for itime in range(60 *24 * 30):
    print(date, domain.T.z.ma.min(), domain.T.z.ma.max())
    date += timedelta

    sim.zbdy[:] = 0.001 * tpxo.update(date)
    sim.update_sealevel_boundaries(timestep)

    sim.update_surface_pressure_gradient(domain.T.z, sp)
    sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
    sim.U.update_halos()
    sim.V.update_halos()
    sim.update_sealevel(timestep, sim.U, sim.V)
    sim.update_depth()
