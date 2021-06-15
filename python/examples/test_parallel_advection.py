import argparse
import timeit
import cProfile

import numpy
import pygetm
import pygetm.parallel

import matplotlib.pyplot

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--debug', action='store_true', help='diagnose individual subdomains')
parser.add_argument('-i', '--interval', type=int, help='number of timesteps between plots', default=10)
parser.add_argument('--nmax', type=int, help='total number of timesteps', default=None)
parser.add_argument('--nrow', type=int, help='number of rows in subdomain division', default=2)
parser.add_argument('--ncol', type=int, help='number of columns in subdomain division', default=2)
parser.add_argument('--noplot', action='store_false', dest='plot', help='skip plotting (useful for performance testing)')
parser.add_argument('--profile', action='store_true', help='use profiler to time function calls')
args = parser.parse_args()

halo = 2
Lx, Ly = 100., 100.
nx, ny, nlev = 612, 600, 1

tiling = pygetm.parallel.Tiling(args.nrow, args.ncol)
rank = tiling.rank

global_domain, f_glob = None, None
if rank == 0:
    # Set up global domain and initial tracer field
    global_domain = pygetm.domain.Domain.create_cartesian(numpy.linspace(-Lx/2, Lx/2, nx), numpy.linspace(-Ly/2, Ly/2, ny), nlev, H=1, f=0.)
    global_domain.initialize(runtype=1)
    f_glob = global_domain.T.array()
    f_glob[int(0.2 * ny):int(0.4 * ny), int(0.2 * nx):int(0.4 * nx)] = 5.

# Set up local subdomain (this also calls subdomain.initialize)
subdomain = pygetm.domain.Domain.partition(tiling, nx, ny, nlev, global_domain, runtype=1)
halo = subdomain.halo

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

u = -omega * subdomain.U.y
v = omega * subdomain.V.x
u[(2 * subdomain.U.x / Lx)**2 + (2 * subdomain.U.y / Ly)**2 > 1] = 0.
v[(2 * subdomain.V.x / Lx)**2 + (2 * subdomain.V.y / Ly)**2 > 1] = 0.

if args.nmax:
    Nmax = args.nmax

# Set up tracer field for subdomain, wrap it for halo updates and MPI-scatter it from root node
f = subdomain.T.array(fill=0.)
f.scatter(f_glob)

# Gather and plot global velocities
u_glob = u.gather(domain=global_domain)
v_glob = v.gather(domain=global_domain)
if u_glob is not None and args.plot:
    fig, ax = matplotlib.pyplot.subplots()
    u_destag, v_destag = u_glob.interp(global_domain.T), v_glob.interp(global_domain.T)
    ax.quiver(u_destag[::10, ::10], v_destag[::10, ::10], angles='xy')
    fig.savefig('vel.png')

if args.debug:
    # Plot local velocities
    fig, ax = matplotlib.pyplot.subplots()
    u_destag, v_destag = u.interp(subdomain.T), v.interp(subdomain.T)
    ax.quiver(subdomain.T.x[::10, ::10], subdomain.T.y[::10, ::10], u_destag[::10, ::10], v_destag[::10, ::10], angles='xy')
    fig.savefig('vel_%i.png' % rank)

    # Set up figure for plotting tracer per subdomain
    fig_sub, ax_sub = matplotlib.pyplot.subplots()
    pc_sub = ax_sub.pcolormesh(subdomain.T.xi, subdomain.T.yi, f)
    cb_sub = fig_sub.colorbar(pc_sub)

# Set up figure for plotting global tracer field
if f_glob is not None and args.plot:
    fig, ax = matplotlib.pyplot.subplots()
    pc = ax.pcolormesh(global_domain.T.xi, global_domain.T.yi, f_glob)
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
                # Print tracer max along boundaries, inside and outside halo 
                print(i, rank, 'inside', f[0, :].max(), f[-1, :].max(), f[:, 0].max(), f[:, -1].max(), flush=True)
                print(i, rank, 'outside', rank, f.all_values[halo-1, :].max(), f.all_values[-halo, :].max(), f.all_values[:, halo-1].max(), f.all_values[:, -halo].max(), flush=True)

                # Plot local tracer field
                pc_sub.set_array(f[...].ravel())
                fig_sub.savefig('subadv_%i_%04i.png' % (rank, ifig))

            # Gather and plot global tracer field
            if args.plot:
                f.gather(out=f_glob)
                if f_glob is not None:
                    pc.set_array(f_glob[...].ravel())
                    fig.savefig('adv_%04i.png' % ifig)

            ifig += 1

        # Advect
        adv.calculate(u, v, timestep, f)

        # Update halos
        f.update_halos()

    print('%.4f s per iteration' % ((timeit.default_timer() - start) / Nmax,))

if args.profile:
    cProfile.run('main()', sort='tottime')
else:
    main()