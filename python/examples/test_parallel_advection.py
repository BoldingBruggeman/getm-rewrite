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
nlev = 1

domain_shape = 600, 604  # nj, ni

rank = MPI.COMM_WORLD.Get_rank()

if rank == 0:
    # Set up global 2D coordinate arrays, bathymetry and initial tracer field
    c1_glob = numpy.broadcast_to(numpy.linspace(-Lx/2, Lx/2, domain_shape[1])[None, :], domain_shape)
    c2_glob = numpy.broadcast_to(numpy.linspace(-Ly/2, Ly/2, domain_shape[0])[:, None], domain_shape)
    H_glob = numpy.zeros(domain_shape)
    H_glob[1:-1, 1:-1] = 1.
    f_glob = numpy.zeros(domain_shape)
    f_glob[int(0.2 * domain_shape[-2]):int(0.4 * domain_shape[-2]), int(0.2 * domain_shape[-1]):int(0.4 * domain_shape[-1])] = 5.
else:
    c1_glob, c2_glob, H_glob, mask_glob, f_glob = None, None, None, None, None

tiling = pygetm.parallel.Tiling(args.nrow, args.ncol)

# Set up local subdomain 
domain = pygetm.Domain(1, nlev, 1, domain_shape[0] // tiling.nrow, 1, domain_shape[1] // tiling.ncol)
halo = domain.halo

# Scatter global bathymetry to subdomains
tiling.wrap(domain.T.H_, halo=halo).scatter(H_glob)

# Scatter global 2D coordinates to subdomains and extract 1D c1, c2
c1_ = numpy.zeros_like(domain.T.H_)
c2_ = numpy.zeros_like(domain.T.H_)
tiling.wrap(c1_, halo=halo).scatter(c1_glob)
tiling.wrap(c2_, halo=halo).scatter(c2_glob)
domain.T.c1_[:] = c1_[halo, :]
domain.T.c2_[:] = c2_[:, halo]

# Infer local mask from bathymetry
domain.T.mask_[...] = domain.T.H_ == 1.

#print(rank, domain.T.c1_[0], domain.T.c1_[-1], domain.T.c2_[0], domain.T.c2_[-1])

domain.initialize()

# Knut's mask hack, implemnted by setting up global mask and scattering that
domain.T.mask_[...] = 0  # this ensure that outer halos [outside global domain] are masked too - should not be needed
if rank == 0:
    mask_glob = numpy.zeros(H_glob.shape, dtype=int)
    mask_glob[2:-2, 2:-2] = 1
tiling.wrap(domain.T.mask_, halo=halo).scatter(mask_glob)

# Set up velocities
period = 600
omega = 2 * numpy.pi / period
cfl = 1.
umax = omega * Lx / 2
dt_cfl = cfl * min(Lx / (domain_shape[1] + 1), Ly / (domain_shape[0] + 1)) / umax
no_of_revolutions = 5
Nmax = no_of_revolutions * round(2 * numpy.pi / omega / dt_cfl)
tmax = no_of_revolutions * 2 * numpy.pi / omega
timestep = tmax / Nmax

# Update U and V mask in halos
tiling.wrap(domain.U.mask_, halo=halo).update_halos()
tiling.wrap(domain.V.mask_, halo=halo).update_halos()
tiling.wrap(domain.U.H_, halo=halo).update_halos()
tiling.wrap(domain.V.H_, halo=halo).update_halos()
tiling.wrap(domain.U.D_, halo=halo).update_halos()
tiling.wrap(domain.V.D_, halo=halo).update_halos()

# These extra halo updates should not be needed - just for testing
# tiling.wrap(u_, halo=halo).update_halos()
# tiling.wrap(v_, halo=halo).update_halos()
# tiling.wrap(domain.T.mask_, halo=halo).update_halos()
# tiling.wrap(domain.U.dx_, halo=halo).update_halos()
# tiling.wrap(domain.V.dx_, halo=halo).update_halos()
# tiling.wrap(domain.T.dx_, halo=halo).update_halos()
# tiling.wrap(domain.U.dy_, halo=halo).update_halos()
# tiling.wrap(domain.V.dy_, halo=halo).update_halos()
# tiling.wrap(domain.T.dy_, halo=halo).update_halos()

f, f_ = domain.array(fill=0.)
u_ = -omega * c2_
v_ = omega * c1_
u_[domain.U.mask_ == 0] = 0.
v_[domain.V.mask_ == 0] = 0.

# Wrap tracer for halo updates
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
    global Nmax
    adv = pygetm.Advection(domain, scheme=6)
    if args.nmax:
        Nmax = args.nmax
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
                f_glob = distf.gather()
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