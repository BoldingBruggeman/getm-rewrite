import logging
import timeit
import functools
from collections.abc import Callable
import argparse
import glob

import numpy as np

import pygetm.core
import pygetm.domain
import pygetm.vertical_coordinates
import pygetm.parallel
import pygetm.operators

handler = logging.StreamHandler()
handler.setLevel(level=logging.ERROR)
logging.basicConfig(handlers=(handler,))

print(f"Importing pygetm from {pygetm.__file__}")


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


def parse_results(fn: str) -> dict[str, float]:
    with open(fn, "r") as f:
        return {line.split(":")[0]: float(line.split(":")[1]) for line in f}


def write_results(outfile: str, results: dict[str, float]) -> None:
    with open(outfile, "w") as f:
        for name, elapsed in results.items():
            f.write(f"{name}: {elapsed}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run advection tests.")
    subparsers = parser.add_subparsers(required=True, dest="command")
    parser_run = subparsers.add_parser("run", help="Run the performance tests")
    parser_run.add_argument(
        "out", type=str, help="Output file to write test results to", nargs="?"
    )
    parser_compare = subparsers.add_parser(
        "compare", help="Compare performance test results"
    )
    parser_compare.add_argument(
        "ref", type=str, help="Reference results file to compare against"
    )
    parser_compare.add_argument(
        "new", type=str, help="New results file to compare against"
    )
    parser_average = subparsers.add_parser(
        "average", help="Average performance test results"
    )
    parser_average.add_argument(
        "files", type=str, nargs="+", help="Performance test result files to average"
    )
    parser_average.add_argument(
        "--out", type=str, help="Output file to write test results to"
    )
    args = parser.parse_args()

    if args.command == "run":
        result = {}
        for name, test in tests:
            elapsed = test()
            result[name] = elapsed
            print(f"{name}: {elapsed:.3f}", flush=True)

        if args.out:
            write_results(args.out, result)

    elif args.command == "compare":
        previous_results = parse_results(args.ref)
        new_results = parse_results(args.new)
        for name in previous_results:
            if name in new_results:
                change = (
                    100
                    * (new_results[name] - previous_results[name])
                    / previous_results[name]
                )
                print(
                    f"{name}: previous = {previous_results[name]:.3f} s, new = {new_results[name]:.3f} s, change = {change:.1f} %"
                )
    elif args.command == "average":
        all_results = []
        for fn in args.files:
            all_results += [parse_results(f) for f in glob.glob(fn)]
        n = len(all_results)
        averaged_results = {}
        for name in all_results[0]:
            if all(name in result for result in all_results):
                averaged_results[name] = sum(r[name] for r in all_results) / n
        for name, avg in averaged_results.items():
            print(f"{name}: {avg:.3f}")
        if args.out:
            write_results(args.out, averaged_results)
