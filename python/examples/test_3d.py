import datetime
import numpy
import pygetm

# Set up rectangular domain with outer points masked
domain = pygetm.domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 50, lat=0, H=50)
sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D, advection_scheme=1)

nc = sim.output_manager.add_netcdf_file('test_3d.nc')
nc.request(['ww', 'uk', 'vk', 'zt'])

# Idealized surface forcing
tausx = domain.U.array(fill=0.01)
tausy = domain.V.array(fill=0.)
sp = domain.T.array(fill=0.)

idpdx = domain.U.array(fill=0., z=pygetm.CENTERS)
idpdy = domain.V.array(fill=0., z=pygetm.CENTERS)
viscosity = domain.T.array(fill=0., z=pygetm.INTERFACES)

# Time
timestep = 10.
ntime = int(3600. // timestep)

plotting_interval = 5
times = timestep * numpy.arange(ntime)
mode_split = 10
domain.T.zio.all_values[...] = 0
domain.T.zin.all_values[...] = 0
sim.start(datetime.datetime(2000, 1, 1), timestep, mode_split)

for istep, time in enumerate(times):
    sim.update_surface_pressure_gradient(domain.T.z, sp)
    sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
    sim.update_sealevel(timestep, sim.U, sim.V)
    sim.update_depth()

    if istep % mode_split == 0:
        print('sim.uvw_momentum_3d')
        sim.Ui.all_values[...] /= mode_split
        sim.Vi.all_values[...] /= mode_split
        sim.start_3d()
        domain.do_vertical()
        sim.update_surface_pressure_gradient(domain.T.zio, sp)
        print('Ui', sim.Ui.values.min(), sim.Ui.values.max())
        print('Vi', sim.Vi.values.min(), sim.Vi.values.max())
        sim.uvw_momentum_3d(timestep * mode_split, tausx, tausy, sim.dpdx, sim.dpdy, idpdx, idpdy, viscosity)
        print('pk', sim.pk.values.min(), sim.pk.values.max())
        print('qk', sim.qk.values.min(), sim.qk.values.max())
        print('ww', sim.ww.values.min(), sim.ww.values.max())
        sim.Ui.all_values[...] = 0
        sim.Vi.all_values[...] = 0
        sim.output_manager.save()
sim.output_manager.close()