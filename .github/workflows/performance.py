import logging
import timeit
import functools
from collections.abc import Callable
import argparse

import numpy as np

import pygetm.core
import pygetm.domain
import pygetm.vertical_coordinates
import pygetm.parallel
import pygetm.operators

handler = logging.StreamHandler()
handler.setLevel(level=logging.ERROR)
logging.basicConfig(handlers=(handler,))

rng = np.random.default_rng()


def create_grid(nx: int, ny: int, nz: int, **kwargs) -> pygetm.core.Grid:
    EXTENT = 50000
    domain = pygetm.domain.create_cartesian(
        np.linspace(0, EXTENT, nx + 1),
        np.linspace(0, EXTENT, ny + 1),
        f=0,
        H=50.0,
        logger=pygetm.parallel.get_logger(level="ERROR"),
        interfaces=True,
    )
    vc = pygetm.vertical_coordinates.Sigma(nz=nz, **kwargs)
    T = domain.create_grids(vc.nz, halox=0, haloy=0, velocity_grids=0)
    vc.initialize(T, logger=domain.root_logger.getChild("vertical_coordinates"))
    vc.update(0.0)
    T.ho = T.array(z=T.hn.z)
    T.ho.all_values[:, :, :] = T.hn.all_values
    return T


test_grids = []
for n in [10, 20, 50]:
    for nz in [10, 25, 50]:
        test_grids.append(create_grid(nx=n, ny=n, nz=nz))


def test_sinking(grid: pygetm.core.Grid):
    # Tracer moves downwards due to sinking (w_var<0)
    # but water remains stationary (w=0)
    op = pygetm.operators.VerticalAdvection(grid)
    w = grid.array(z=pygetm.INTERFACES)
    w_var = grid.array(z=pygetm.INTERFACES)
    w.fill(0.0)
    TIMESTEP = 3600.0
    tracer = grid.array(z=pygetm.CENTERS)

    # Just below CFL
    w_var.fill(-0.99 * grid.ho.all_values[-1, 0, 0] / TIMESTEP)
    tracer.fill(0.0)
    tracer.all_values[-1, :, :] = 1.0
    return timeit.timeit(lambda: op(w, w_var, TIMESTEP, tracer), number=10000)


tests: list[tuple[str, Callable[[], float]]] = []
for grid in test_grids:
    tests.append(
        (f"test_sinking{grid.ho.shape}", functools.partial(test_sinking, grid))
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run advection tests.")
    parser.add_argument(
        "outfile", type=str, help="Output file to write test results to", nargs="?"
    )
    parser.add_argument(
        "--compare", type=str, help="Previous results file to compare against"
    )
    args = parser.parse_args()

    if args.compare:
        with open(args.compare, "r") as f:
            previous_results = {
                line.split(":")[0]: float(line.split(":")[1]) for line in f
            }

    result = {}
    for name, test in tests:
        elapsed = test()
        result[name] = elapsed
        postfix = ""
        if args.compare and name in previous_results:
            change = 100 * (elapsed - previous_results[name]) / previous_results[name]
            postfix = (
                f" (previous: {previous_results[name]:.3f} s) Change: {change:.1f} %"
            )
        print(f"{name}: {elapsed:.3f} s{postfix}")

    if args.outfile:
        with open(args.outfile, "w") as f:
            for name, test in tests:
                elapsed = result[name]
                f.write(f"{name}: {elapsed}\n")
