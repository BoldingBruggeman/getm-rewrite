import argparse
import timeit
import cProfile

import numpy
import pygetm
import pygetm.parallel
from mpi4py import MPI

import matplotlib.pyplot

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--debug', action='store_true', help='diagnose individual subdomains')
parser.add_argument('-i', '--interval', type=int, help='number of timesteps between plots', default=10)
parser.add_argument('--nmax', type=int, help='total number of timesteps', default=None)
parser.add_argument('--nrow', type=int, help='number of rows in subdomain division', default=2)
parser.add_argument('--ncol', type=int, help='number of columns in subdomain division', default=2)
parser.add_argument('--noplot', action='store_true', help='skip plotting (useful for performance testing)')
parser.add_argument('--profile', action='store_true', help='use profiler to time function calls')
args = parser.parse_args()

halo = 2
Lx, Ly = 100., 100.
nx, ny, nlev = 612, 600, 1

rank = MPI.COMM_WORLD.Get_rank()

global_domain, f_glob = None, None
if rank == 0:
    # Set up global 2D coordinate arrays, bathymetry and initial tracer field
    global_domain = pygetm.Domain.create_cartesian(numpy.linspace(-Lx/2, Lx/2, nx), numpy.linspace(-Ly/2, Ly/2, ny), nlev, H=0.)
    global_domain.T.H[1:-1, 1:-1] = 1.
    global_domain.T.mask[1:-1, 1:-1] = 1
    global_domain.initialize()
    f_glob, _ = global_domain.array()
    f_glob[int(0.2 * ny):int(0.4 * ny), int(0.2 * nx):int(0.4 * nx)] = 5.

tiling = pygetm.parallel.Tiling(args.nrow, args.ncol)

# Set up local subdomain 
subdomain = pygetm.Domain.partition(tiling, nx, ny, nlev, global_domain)
halo = subdomain.halo

# Knut's mask hack, implemented by setting up global mask and scattering that
mask_glob = None
subdomain.T.mask_[...] = 0  # this ensure that outer halos [outside global domain] are masked too - should not be needed
if rank == 0:
    mask_glob, _ = global_domain.array(dtype=int)
    mask_glob[2:-2, 2:-2] = 1
tiling.wrap(subdomain.T.mask_, halo=halo).scatter(mask_glob)

# Set up velocities
period = 600
omega = 2 * numpy.pi / period
cfl = 1.
umax = omega * Lx / 2
dt_cfl = cfl * min(Lx / nx, Ly / ny) / umax
no_of_revolutions = 5
Nmax = no_of_revolutions * round(2 * numpy.pi / omega / dt_cfl)
tmax = no_of_revolutions * 2 * numpy.pi / omega
timestep = tmax / Nmax

u_ = -omega * subdomain.T.y_
v_ = omega * subdomain.T.x_
u_[subdomain.U.mask_ == 0] = 0.
v_[subdomain.V.mask_ == 0] = 0.

if args.nmax:
    Nmax = args.nmax

# Set up tracer field for subdomain, wrap it for halo updates, and MPI-scatter it from root node
f, f_ = subdomain.array(fill=0.)
distf = tiling.wrap(f_, halo=halo)
distf.scatter(f_glob)

# Gather and plot global velocities
u_glob = tiling.wrap(u_, halo=halo).gather()
v_glob = tiling.wrap(v_, halo=halo).gather()
if u_glob is not None and not args.noplot:
    fig = matplotlib.pyplot.figure()
    ax = fig.gca()
    ax.quiver(u_glob[::10, ::10], v_glob[::10, ::10], angles='xy')
    fig.savefig('vel.png')

if args.debug:
    # Plot local velocities
    fig = matplotlib.pyplot.figure()
    ax = fig.gca()
    ax.quiver(u_[::10, ::10], v_[::10, ::10], angles='xy')
    fig.savefig('vel_%i.png' % rank)

    # Set up figure for plotting tracer per subdomain
    fig_sub = matplotlib.pyplot.figure()
    ax_sub = fig_sub.gca()
    pc_sub = ax_sub.pcolormesh(f_)
    cb_sub = fig_sub.colorbar(pc_sub)

# Set up figure for plotting global tracer field
if f_glob is not None and not args.noplot:
    fig = matplotlib.pyplot.figure()
    ax = fig.gca()
    pc = ax.pcolormesh(f_glob)
    cb = fig.colorbar(pc)

def main():
    adv = pygetm.Advection(subdomain, scheme=1)
    ifig = 0
    start = timeit.default_timer()
    for i in range(Nmax):
        if i % args.interval == 0:
            if rank == 0:
                print('time step %i of %i' % (i, Nmax), flush=True)

            if args.debug:
                # Print tracer max along boundaries, inside and outsde halo 
                print(i, rank, 'inside', f[0, :].max(), f[-1, :].max(), f[:, 0].max(), f[:, -1].max(), flush=True)
                print(i, rank, 'outside', rank, f_[halo-1, :].max(), f_[-halo, :].max(), f_[:, halo-1].max(), f_[:, -halo].max(), flush=True)

                # Plot local tracer field
                pc_sub.set_array(f_.ravel())
                fig_sub.savefig('subadv_%i_%04i.png' % (rank, ifig))

            # Gather and plot global tracer field
            if not args.noplot:
                distf.gather(out=f_glob)
                if f_glob is not None:
                    pc.set_array(f_glob.ravel())
                    fig.savefig('adv_%04i.png' % ifig)

            ifig += 1

        # Advect
        adv.calculate(u_, v_, timestep, f_)

        # Update halos
        distf.update_halos()

    print('%.4f s per iteration' % ((timeit.default_timer() - start) / Nmax,))

if args.profile:
    cProfile.run('main()', sort='tottime')
else:
    main()