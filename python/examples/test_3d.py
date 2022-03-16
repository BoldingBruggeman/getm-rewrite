import numpy
import pygetm

# Set up rectangular domain with outer points masked
domain = pygetm.domain.create_cartesian(500.*numpy.arange(100), 500.*numpy.arange(30), 50, lat=0, H=50)
sim = pygetm.Simulation(domain, runtype=pygetm.BAROTROPIC_3D, advection_scheme=1)

# Idealized surface forcing
tausx = domain.U.array(fill=0.01)
tausy = domain.V.array(fill=0.)
sp = domain.T.array(fill=0.)

idpdx = domain.U.array(fill=0., is_3d=True)
idpdy = domain.V.array(fill=0., is_3d=True)
viscosity = domain.W.array(fill=0., is_3d=True)

# Time
timestep = 10.
ntime = int(3600. // timestep)

plotting_interval = 5
times = timestep * numpy.arange(ntime)
mode_split = 10
for istep, time in enumerate(times):
    sim.update_surface_pressure_gradient(domain.T.z, sp)
    sim.uv_momentum_2d(timestep, tausx, tausy, sim.dpdx, sim.dpdy)
    sim.update_sealevel(timestep, sim.U, sim.V)
    sim.update_depth()

    if istep % mode_split == 0:
        print('sim.uvw_momentum_3d')
        sim.uvw_momentum_3d(timestep, tausx, tausy, sim.dpdx, sim.dpdy, idpdx, idpdy, viscosity)
