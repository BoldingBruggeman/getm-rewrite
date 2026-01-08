import numpy as np
import logging

from pygetm import domain
from pygetm import vertical_coordinates
from pygetm.constants import CENTERS, INTERFACES

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["lines.linewidth"] = 0.05

nx = 100
ny = 1
nz = 25
nz = 50
nz = 10

x = np.linspace(0.0, 100000, nx + 1)
y = np.linspace(0.0, 100, ny + 1)
# x = np.linspace(0.0, 100000, nx)
# y = np.linspace(0.0, 100, ny)
H = 5.0 + 1 * np.arange(nx)
if True:
    A = 0.6
    n0 = nx // 2
    L = 10000
    r = (x - x[n0]) ** 2 / 10000**2
    M = 40 * (-A * np.exp(-r))
    H += M[1:]


domain = domain.create_cartesian(x=x, y=y, f=0.0, H=H, interfaces=True)
# domain = domain.create_cartesian(x=x, y=y, f=0.0, H=H)
grid = domain.create_grids(
    nz=nz,
    halox=0,
    haloy=0,
    t_postfix="t",
    velocity_grids=1,
)
grid.ho = grid.array(fill_value=0.0, z=CENTERS)
grid.array(name="NN", fill_value=0.0, z=INTERFACES)
grid.array(name="SS", fill_value=0.0, z=INTERFACES)

# Adaptive Coordinates
timestep = 20 * 60.0
cnpar = 1.0
ddl = 0.00001
ddu = 0.00001
ddu = 1.0
ddl = 1.0
Dgamma = 30.0
gamma_surf = True

vc = vertical_coordinates.Adaptive(
    nz,
    cnpar=cnpar,
    ddu=ddu,
    ddl=ddl,
    gamma_surf=gamma_surf,
    Dgamma=Dgamma,
    csigma=0.001,
    cgvc=0.001,
    hpow=3,
    chsurf=-0.001,
    hsurf=1.5,
    chmidd=-0.1,
    hmidd=0.5,
    chbott=-0.001,
    hbott=1.5,
    cneigh=-0.1,
    rneigh=0.25,
    decay=2.0 / 3.0,
    cNN=-1.0,
    drho=0.3,
    cSS=-1.0,
    dvel=0.1,
    chmin=-0.1,
    hmin=2.5,
    nvfilter=-1,
    vfilter=0.2,
    nhfilter=-1,
    hfilter=0.1,
    split=1,
    timescale=1.0 * 3600.0,
)
logger = logging.getLogger(__name__)
vc.initialize(grid, logger=logger)

print(f"nx={nx}, ny={ny}, nz={nz}")
print(f"depth: {H.min()} to {H.max()}")
print(f"ddl={ddl}, ddu={ddu}, Dgamma={Dgamma}")
print("note - specified on super grid")

vc.update()
grid.ho.all_values = grid.hn.all_values

fig, axs = plt.subplots(3, 2)
first = True
xp = 50
surf = []
midd = []
bott = []
for i in range(10000):
    #vc(grid.D[...], grid.hn[...])
    vc.update(timestep)
    grid.ho.all_values = grid.hn.all_values

    if i % 100 == 0:
        fig.suptitle(f"nug, dga and hn (nz={nz}): {i*timestep/3600:.1f} hours")
        f0 = axs[0][0].plot(
            vc.nug[1:, 0, ::10], -np.cumsum(grid.hn[:, 0, ::10], axis=0), marker=None
        )
        f1 = axs[0][1].pcolormesh(
            x[1:], -np.cumsum(grid.hn[:, 0, :], axis=0), vc.nug[1:, 0, :]
        )
        if first:
            plt.colorbar(f1, ax=axs[0][1])

        f2 = axs[1][0].plot(
            vc.dga_t[:, 0, ::10], -np.cumsum(grid.hn[:, 0, ::10], axis=0), marker=None
        )
        f3 = axs[1][1].pcolormesh(vc.dga_t[:, 0, :])
        if first:
            plt.colorbar(f3, ax=axs[1][1])

        f4 = axs[2][0].plot(
            grid.hn[:, 0, ::10], -np.cumsum(grid.hn[:, 0, ::10], axis=0), marker=None
        )
        # f5 = axs[2][1].pcolormesh(x[1:], -np.cumsum(grid.hn[:,0,:], axis=0), grid.hn[:,0,:])
        # if first:
        #    plt.colorbar(f5, ax=axs[2][1])
        f5 = axs[2][1].plot(
            x[1::2], -np.cumsum(grid.hn[:, 0, ::2], axis=0).T, linewidth=0.1
        )
        first = False
        plt.pause(0.1)
        # plt.show()

        surf.append(grid.hn[nz-1, 0, xp])
        midd.append(grid.hn[nz//2, 0, xp])
        bott.append(grid.hn[0, 0, xp])

if False:
    print("nug ", vc.nug[1:, 0, 25])
    print("ga  ", vc.ga[:, 0, 5])
    print("dga ", vc.dga_t[:, 0, 5])
    print("hn  ", grid.hn[:, 0, 21])
    print(grid.hn[:, 0, 21].sum())

    print(vc.dga_t[:, 0, 25::50])

plt.show()

mpl.rcParams["lines.linewidth"] = 0.5
plt.title(f"layer height (x={xp})")
plt.xlabel(f"time [hours]")
plt.ylabel(f"hn [m]")
x=100*timestep*np.arange(len(surf))/3600
plt.plot(x, surf, label='Surface')
plt.plot(x, midd, label='Middle')
plt.plot(x, bott, label='Bottom')
plt.legend(loc="upper right")
plt.show()
