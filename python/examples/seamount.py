import datetime

import numpy as np
import cftime
import pygetm

# Parameters for horizontal grid and bottom topography
nx = 66
ny = 66
nz = 20
dx = 8000.0
dy = 8000.0
HMAX = 4500.0
A = 0.6
L = 50000.0

# horizontal grid
x = np.linspace(-0.5 * nx, 0.5 * nx, nx) * dx
y = np.linspace(-0.5 * ny, 0.5 * ny, ny) * dy

# bottom topography
r = -(x[np.newaxis, :] ** 2 + y[:, np.newaxis] ** 2) / L**2
H = HMAX * (1.0 - A * np.exp(r))

domain = pygetm.domain.create_cartesian(
    x, y, f=0.0, H=H, z0=0.01, periodic_x=True, periodic_y=True
)

sim = pygetm.Simulation(
    domain,
    airsea=pygetm.airsea.Fluxes(taux=0.0),
    gotm="../../../getm-setups/seamount/gotmturb.nml",
    #internal_pressure=pygetm.internal_pressure.BlumbergMellor(),
    vertical_coordinates=pygetm.vertical_coordinates.Sigma(nz),
    Dcrit=0.1,
    Dmin=0.02,
)

sim.salt.fill(35.0)
sim.temp.fill(5.0 + 15.0 * np.exp(sim.T.zc / 1000.0))

sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)

output = sim.output_manager.add_netcdf_file(
    "seamount.nc", interval=datetime.timedelta(hours=1)
)
output.request("zt", "temp", "Ht", "idpdx", "idpdy", "SxB", "SyB")

stop = cftime.datetime(2001, 1, 4)
sim.start(cftime.datetime(2001, 1, 1), 12.0, 30, report=datetime.timedelta(hours=1))
while sim.time < stop:
    sim.advance()
sim.finish()
