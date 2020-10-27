import numpy
import pygetm
import pygetm.parallel
from mpi4py import MPI

import matplotlib.pyplot

halo = 2
Lx, Ly = 100., 100.
nlev = 1

domain_shape = 400, 600  # nj, ni

rank = MPI.COMM_WORLD.Get_rank()

if rank == 0:
    # Set up global 2D coordinate arrays and bathymetry
    c1_glob = numpy.broadcast_to(numpy.linspace(-Lx/2, Lx/2, domain_shape[1])[None, :], domain_shape)
    c2_glob = numpy.broadcast_to(numpy.linspace(-Ly/2, Ly/2, domain_shape[0])[:, None], domain_shape)
    H_glob = numpy.zeros(domain_shape)
    H_glob[1:-1, 1:-1] = 1.
else:
    c1_glob, c2_glob, H_glob, mask_glob = None, None, None, None

tiling = pygetm.parallel.Tiling(2, 2)

# Set up local subdomain 
domain = pygetm.Domain(1, nlev, 1, domain_shape[0] // tiling.nrow, 1, domain_shape[1] // tiling.ncol)
halo = domain.halo

# Scatter global bathymetry to subdomains
tiling.wrap(domain.T.H_, halo=halo).scatter(H_glob)

# Scatter global coordinates to subdomains and extract 1D c1, c2
c1_ = numpy.zeros_like(domain.T.H_)
c2_ = numpy.zeros_like(domain.T.H_)
tiling.wrap(c1_, halo=halo).scatter(c1_glob)
tiling.wrap(c2_, halo=halo).scatter(c2_glob)
domain.T.c1_[:] = c1_[halo, :]
domain.T.c2_[:] = c2_[:, halo]

# Infer local mask from bathymetry
domain.T.mask_[...] = domain.T.H_ == 1.

#print(rank, domain.T.c1_.min(), domain.T.c1_.max(), domain.T.c2_.min(), domain.T.c2_.max())

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

u, u_ = domain.array(fill=0.)
v, v_ = domain.array(fill=0.)
var, var_ = domain.array(fill=1.)
u_[:, :] = -omega * c2_
v_[:, :] = omega * c1_
u_[domain.U.mask_ == 0] = 0.
v_[domain.V.mask_ == 0] = 0.
if rank == 0:
    var[int(tiling.nrow * 0.2 * domain.shape[1]):int(tiling.nrow * 0.4 * domain.shape[1]), int(tiling.ncol * 0.2 * domain.shape[2]):int(tiling.ncol * 0.4 * domain.shape[2])] = 5.

# These extra halo updates should not be needed - just for testing
tiling.wrap(u_, halo=halo).update_halos()
tiling.wrap(v_, halo=halo).update_halos()
tiling.wrap(domain.U.mask_, halo=halo).update_halos()
tiling.wrap(domain.V.mask_, halo=halo).update_halos()
tiling.wrap(domain.T.mask_, halo=halo).update_halos()

# Gather and plot global velocities
u_glob = tiling.wrap(u_, halo=halo).gather()
v_glob = tiling.wrap(v_, halo=halo).gather()
if u_glob is not None:
    fig = matplotlib.pyplot.figure()
    ax = fig.gca()
    ax.quiver(u_glob[::10, ::10], v_glob[::10, ::10], angles='xy')
    fig.savefig('vel.png')
fig = matplotlib.pyplot.figure()
ax = fig.gca()
ax.quiver(u_[::10, ::10], v_[::10, ::10], angles='xy')
fig.savefig('vel_%i.png' % rank)

# Wrap tracer for halo updates
distvar = tiling.wrap(var_, halo=halo)

# Set up figure for plotting tracer per subdomain
fig_sub = matplotlib.pyplot.figure()
ax_sub = fig_sub.gca()
pc_sub = ax_sub.pcolormesh(var_)
cb_sub = fig_sub.colorbar(pc_sub)

# Set up figure for plotting global tracer field
var_glob = distvar.gather()
if var_glob is not None:
    fig = matplotlib.pyplot.figure()
    ax = fig.gca()
    pc = ax.pcolormesh(var_glob)
    cb = fig.colorbar(pc)

adv = pygetm.Advection(domain, scheme=1)
ifig = 0
for i in range(Nmax):
    if i % 10 == 0:
        # Print tracer max along boundaries, inside and outsde halo 
        print(i, rank, 'inside', var[0, :].max(), var[-1, :].max(), var[:, 0].max(), var[:, -1].max(), flush=True)
        print(i, rank, 'outside', rank, var_[halo-1, :].max(), var_[-halo, :].max(), var_[:, halo-1].max(), var_[:, -halo].max(), flush=True)

        # Plot local tracer field
        pc_sub.set_array(var_.ravel())
        fig_sub.savefig('subadv_%i_%04i.png' % (rank, ifig))

        # Gather and plot global tracer field
        var_glob = distvar.gather()
        if var_glob is not None:
            pc.set_array(var_glob.ravel())
            fig.savefig('adv_%04i.png' % ifig)

        ifig += 1

    # Advect
    adv.calculate(u_, v_, timestep, var_)

    # Update halos
    distvar.update_halos()


