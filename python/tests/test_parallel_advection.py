import argparse
import timeit
import cProfile
import unittest
from typing import Optional
import sys
import hashlib
import gc

import numpy as np
import pygetm
import pygetm.parallel
import pygetm.util.compare_nc


def calculate_md5(path: str) -> str:
    hash = hashlib.md5()
    with open(path, "rb") as f:
        hash.update(f.read())
    return hash.hexdigest()


class TestParallelAdvection(unittest.TestCase):
    args = None

    def test(self):
        kwargs = {}
        if self.args is not None:
            kwargs.update(
                nmax=self.args.nmax,
                output=self.args.output,
                debug=self.args.debug,
                plot=self.args.plot,
                n=self.args.n,
                profile=self.args.profile,
            )
        combos = [(1, 1)]
        if pygetm.parallel.MPI.COMM_WORLD.size >= 2:
            combos += [(1, 2), (2, 1)]
        if pygetm.parallel.MPI.COMM_WORLD.size >= 4:
            combos += [(2, 2), (1, 4), (4, 1)]
        for scheme in [pygetm.AdvectionScheme.DEFAULT]:
            ref_output = None
            for nrow, ncol in combos:
                with self.subTest(scheme=scheme.name, nrow=nrow, ncol=ncol):
                    tiling = pygetm.parallel.Tiling(nrow, ncol, ncpus=nrow * ncol)
                    kwargs["output"] = "par_adv_%ix%i.nc" % (nrow, ncol)
                    self._test(tiling, scheme=scheme, **kwargs)
                    gc.collect()
                    if tiling.rank == 0:
                        if ref_output is None:
                            ref_output = kwargs["output"]
                        else:
                            self.assertTrue(
                                pygetm.util.compare_nc.compare(
                                    ref_output, kwargs["output"]
                                )
                            )

    def _test(
        self,
        tiling: pygetm.parallel.Tiling,
        nmax: Optional[int] = None,
        output: Optional[str] = None,
        debug: bool = False,
        plot: bool = False,
        n: int = 1,
        profile: bool = False,
        scheme: pygetm.AdvectionScheme = pygetm.AdvectionScheme.DEFAULT,
    ):
        Lx, Ly = 100.0, 100.0
        nx, ny = 612, 600

        rank = tiling.rank

        domain = pygetm.domain.create_cartesian(
            np.linspace(-Lx / 2, Lx / 2, nx),
            np.linspace(-Ly / 2, Ly / 2, ny),
            H=1,
            f=0.0,
            comm=tiling.comm,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        T = domain.create_grids(1, halox=2, haloy=2, velocity_grids=1, tiling=tiling)

        if T is None:
            # unused rank
            return

        halox = T.halox
        haloy = T.haloy
        outman = pygetm.output.OutputManager(T.fields, rank=rank)

        # Set up velocities
        period = 600
        omega = 2 * np.pi / period
        cfl = 1.0
        umax = omega * Lx / 2
        dt_cfl = cfl * min(Lx / nx, Ly / ny) / umax
        Nmax = n * round(2 * np.pi / omega / dt_cfl)
        tmax = n * 2 * np.pi / omega
        timestep = tmax / Nmax

        # Calculate u and v
        # Note that this only sets values in the interior of the domain.
        # Therefore, a halo exchange is needed to ensure u and v are also valid
        # in the innermost halo strip (needed for advection scheme)
        u = -omega * T.ugrid.y
        v = omega * T.vgrid.x
        u[(2 * T.ugrid.x / Lx) ** 2 + (2 * T.ugrid.y / Ly) ** 2 >= 1] = 0.0
        v[(2 * T.vgrid.x / Lx) ** 2 + (2 * T.vgrid.y / Ly) ** 2 >= 1] = 0.0
        u.update_halos()
        v.update_halos()

        if nmax is not None:
            Nmax = nmax

        # Set up tracer field for subdomain, wrap it for halo updates
        # and MPI-scatter it from root node
        f = T.array(fill=0.0, name="tracer")
        values_glob = None
        if rank == 0:
            # Root rank: calculate rectangular global initial tracer field
            values_glob = np.zeros((ny, nx))
            imin, imax = int(0.2 * nx), int(0.4 * nx)
            jmin, jmax = int(0.2 * ny), int(0.4 * ny)
            values_glob[jmin:jmax, imin:imax] = 5.0
        f.scatter(values_glob)

        # Gather and plot global velocities
        if plot:
            u_glob = u.interp(T).gather()
            v_glob = v.interp(T).gather()
            if u_glob is not None:
                import matplotlib.pyplot

                fig, ax = matplotlib.pyplot.subplots()
                ax.quiver(u_glob[::10, ::10], v_glob[::10, ::10], angles="xy")
                fig.savefig("vel.png")

        if debug:
            # Plot local velocities
            import matplotlib.pyplot

            fig, ax = matplotlib.pyplot.subplots()
            u_destag, v_destag = u.interp(T), v.interp(T)
            ax.quiver(
                T.x[::10, ::10],
                T.y[::10, ::10],
                u_destag[::10, ::10],
                v_destag[::10, ::10],
                angles="xy",
            )
            fig.savefig("vel_%i.png" % rank)

            # Set up figure for plotting tracer per subdomain
            fig_sub, ax_sub = matplotlib.pyplot.subplots()
            pc_sub = ax_sub.pcolormesh(T.xgrid.x, T.xgrid.y, f)
            cb_sub = fig_sub.colorbar(pc_sub)

        # Set up figure for plotting global tracer field
        if values_glob is not None and plot:
            fig, ax = matplotlib.pyplot.subplots()
            pc = ax.pcolormesh(T.xgrid.x, T.xgrid.y, values_glob)
            cb = fig.colorbar(pc)

        if output:
            ncf = outman.add_netcdf_file(output, interval=10)
            ncf.request("tracer")
        outman.start()

        advect = pygetm.operators.Advection(T, scheme=scheme)

        if profile:
            prof = cProfile.Profile()
            prof.enable()

        ifig = 0
        start = timeit.default_timer()
        for i in range(Nmax):
            if i % int(0.1 * Nmax) == 0:
                if rank == 0:
                    domain.root_logger.info("time step %i of %i" % (i, Nmax))

                if debug:
                    # Print tracer max along boundaries, inside and outside halo
                    print(
                        i,
                        rank,
                        "inside",
                        f[0, :].max(),
                        f[-1, :].max(),
                        f[:, 0].max(),
                        f[:, -1].max(),
                        flush=True,
                    )
                    print(
                        i,
                        rank,
                        "outside",
                        rank,
                        f.all_values[haloy - 1, :].max(),
                        f.all_values[-haloy, :].max(),
                        f.all_values[:, halox - 1].max(),
                        f.all_values[:, -halox].max(),
                        flush=True,
                    )

                    # Plot local tracer field
                    pc_sub.set_array(f[...].ravel())
                    fig_sub.savefig("subadv_%i_%04i.png" % (rank, ifig))

                # Gather and plot global tracer field
                if plot:
                    f.gather(out=values_glob)
                    if values_glob is not None:
                        pc.set_array(values_glob[...].ravel())
                        fig.savefig("adv_%04i.png" % ifig)

                ifig += 1

            # Advect
            advect(u, v, timestep, f)

            outman.save(i * timestep, i)

        duration = timeit.default_timer() - start
        domain.root_logger.info("Time spent in loop: %.4f s" % duration)
        domain.root_logger.info("%.4f ms per iteration" % (1000 * duration / Nmax,))

        if profile:
            prof.disable()
            prof.print_stats("tottime")

        outman.close(Nmax * timestep)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--debug", action="store_true", help="diagnose individual subdomains"
    )
    parser.add_argument(
        "--nmax", type=int, help="total number of timesteps", default=None
    )
    parser.add_argument(
        "--nrow", type=int, help="number of rows in subdomain division", default=1
    )
    parser.add_argument(
        "--ncol", type=int, help="number of columns in subdomain division", default=1
    )
    parser.add_argument("-n", type=int, help="number of revolutions", default=1)
    parser.add_argument(
        "--plot",
        action="store_true",
        help="create stills that show the tracer at every time step",
    )
    parser.add_argument(
        "--profile", action="store_true", help="use profiler to time function calls"
    )
    parser.add_argument("-o", "--output", help="NetCDF file to save result to")
    TestParallelAdvection.args, remaining = parser.parse_known_args()

    unittest.main(argv=sys.argv[:1] + remaining)
