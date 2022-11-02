import argparse
import timeit
import cProfile
import unittest
from typing import Optional
import sys
import hashlib

import numpy as np
import pygetm
import pygetm.parallel


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
        combos = [(1, 1)] + {2: [(1, 2), (2, 1)]}.get(
            pygetm.parallel.MPI.COMM_WORLD.size, []
        )
        for scheme in [pygetm.AdvectionScheme.DEFAULT]:
            ref_output = None
            for nrow, ncol in combos:
                with self.subTest(scheme=scheme.name, nrow=nrow, ncol=ncol):
                    tiling = pygetm.parallel.Tiling(nrow, ncol, ncpus=nrow * ncol)
                    kwargs["output"] = "par_adv_%ix%i.nc" % (nrow, ncol)
                    self._test(tiling, scheme=scheme, **kwargs)
                    if tiling.rank == 0:
                        md5 = calculate_md5(kwargs["output"])
                        if ref_output is None:
                            ref_output = md5
                        else:
                            self.assertEqual(ref_output, md5)

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
        nx, ny, nlev = 612, 600, 1

        rank = tiling.rank

        subdomain = pygetm.domain.create_cartesian(
            np.linspace(-Lx / 2, Lx / 2, nx),
            np.linspace(-Ly / 2, Ly / 2, ny),
            nlev,
            H=1,
            f=0.0,
            tiling=tiling,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

        if subdomain is None:
            # unused rank
            return

        halo = subdomain.halo
        subdomain.initialize(runtype=pygetm.BAROTROPIC_2D)
        outman = pygetm.output.OutputManager(subdomain.fields, rank=rank)

        f_glob = None
        if subdomain.glob:
            # Root rank: set up global domain and initial tracer field
            f_glob = subdomain.glob.T.array(fill=0.0)
            f_glob[int(0.2 * ny) : int(0.4 * ny), int(0.2 * nx) : int(0.4 * nx)] = 5.0

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
        u = -omega * subdomain.U.y
        v = omega * subdomain.V.x
        u[(2 * subdomain.U.x / Lx) ** 2 + (2 * subdomain.U.y / Ly) ** 2 >= 1] = 0.0
        v[(2 * subdomain.V.x / Lx) ** 2 + (2 * subdomain.V.y / Ly) ** 2 >= 1] = 0.0
        u.update_halos()
        v.update_halos()

        if nmax is not None:
            Nmax = nmax

        # Set up tracer field for subdomain, wrap it for halo updates
        # and MPI-scatter it from root node
        f = subdomain.T.array(fill=0.0, name="tracer")
        f.scatter(f_glob)

        # Gather and plot global velocities
        u_glob = u.gather()
        v_glob = v.gather()
        if u_glob is not None and plot:
            import matplotlib.pyplot

            fig, ax = matplotlib.pyplot.subplots()
            u_destag, v_destag = (
                u_glob.interp(u_glob.grid.domain.T),
                v_glob.interp(u_glob.grid.domain.T),
            )
            ax.quiver(u_destag[::10, ::10], v_destag[::10, ::10], angles="xy")
            fig.savefig("vel.png")

        if debug:
            # Plot local velocities
            import matplotlib.pyplot

            fig, ax = matplotlib.pyplot.subplots()
            u_destag, v_destag = u.interp(subdomain.T), v.interp(subdomain.T)
            ax.quiver(
                subdomain.T.x[::10, ::10],
                subdomain.T.y[::10, ::10],
                u_destag[::10, ::10],
                v_destag[::10, ::10],
                angles="xy",
            )
            fig.savefig("vel_%i.png" % rank)

            # Set up figure for plotting tracer per subdomain
            fig_sub, ax_sub = matplotlib.pyplot.subplots()
            pc_sub = ax_sub.pcolormesh(subdomain.T.xi, subdomain.T.yi, f)
            cb_sub = fig_sub.colorbar(pc_sub)

        # Set up figure for plotting global tracer field
        if f_glob is not None and plot:
            fig, ax = matplotlib.pyplot.subplots()
            pc = ax.pcolormesh(f_glob.grid.domain.T.xi, f_glob.grid.domain.T.yi, f_glob)
            cb = fig.colorbar(pc)

        if output:
            ncf = outman.add_netcdf_file(output, interval=10)
            ncf.request("tracer")
        outman.start()

        advect = pygetm.operators.Advection(subdomain.T, scheme=scheme)

        if profile:
            prof = cProfile.Profile()
            prof.enable()

        ifig = 0
        start = timeit.default_timer()
        for i in range(Nmax):
            if i % int(0.1 * Nmax) == 0:
                if rank == 0:
                    subdomain.logger.info("time step %i of %i" % (i, Nmax))

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
                        f.all_values[halo - 1, :].max(),
                        f.all_values[-halo, :].max(),
                        f.all_values[:, halo - 1].max(),
                        f.all_values[:, -halo].max(),
                        flush=True,
                    )

                    # Plot local tracer field
                    pc_sub.set_array(f[...].ravel())
                    fig_sub.savefig("subadv_%i_%04i.png" % (rank, ifig))

                # Gather and plot global tracer field
                if plot:
                    f.gather(out=f_glob)
                    if f_glob is not None:
                        pc.set_array(f_glob[...].ravel())
                        fig.savefig("adv_%04i.png" % ifig)

                ifig += 1

            # Advect
            advect(u, v, timestep, f)

            outman.save(i * timestep, i)

        duration = timeit.default_timer() - start
        subdomain.logger.info("Time spent in loop: %.4f s" % duration)
        subdomain.logger.info("%.4f ms per iteration" % (1000 * duration / Nmax,))

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
    TestAdvection.args, remaining = parser.parse_known_args()

    unittest.main(argv=sys.argv[:1] + remaining)

