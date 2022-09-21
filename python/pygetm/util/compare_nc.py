import sys
from typing import Container, Mapping, Optional
import numpy as np
import netCDF4


def compare(
    path1: str,
    path2: str,
    rotate: bool = False,
    check_constant: bool = True,
    name_map: Mapping[str, str] = {},
    flip: Container[str] = (),
    tolerance: float = 0.0,
) -> bool:
    mismatches = []

    def load(slc=(Ellipsis,)):
        values1 = ncvar1[slc]
        values2 = ncvar2[slc]
        if rotate and ncvar1.ndim >= 2:
            values2 = np.swapaxes(values2, -1, -2)[..., ::-1, :]
            if ncvar1.dimensions[-2:] == ("yv", "xv") and ncvar2.dimensions[-2:] == (
                "yu",
                "xu",
            ):
                values1 = values1[..., :-1, :]
                values2 = values2[..., 1:, :]
            if name in flip:
                values2 = -values2
        return np.ma.masked_invalid(values1), np.ma.masked_invalid(values2)

    with netCDF4.Dataset(path1) as nc1, netCDF4.Dataset(path2) as nc2:
        for name in frozenset(nc1.variables) & frozenset(nc2.variables):
            ncvar1 = nc1.variables[name]
            ncvar2 = nc2.variables[name_map.get(name, name)]
            if "time" in ncvar1.dimensions:
                # Quick check for very last time point; if values match,
                # values at all preceding time points are assumed to match too.
                # Then skip the full values check (expensive!)
                values1, values2 = load((-1, Ellipsis))
                if np.ma.filled(values1 == values2, True).all():
                    continue
            elif not check_constant:
                continue

            # Values at last time did not match. Now check all time points
            values1, values2 = load()
            mismatch = np.ma.filled(values1 != values2, False)
            maxval = max(np.ma.abs(values1).max(), np.ma.abs(values2).max())
            idim_start = 1 if "time" in ncvar1.dimensions else 0
            mismatch_count = mismatch.sum(axis=tuple(range(idim_start, mismatch.ndim)))
            absdiff = np.abs(np.ma.filled(values1 - values2, 0.0))
            mismatch_max = absdiff.max(axis=tuple(range(idim_start, mismatch.ndim)))
            if mismatch_count.sum() > 0:
                ifirst_mismatch = mismatch_count.nonzero()[0][0]
                mismatches.append(
                    (name, mismatch_count, ifirst_mismatch, mismatch_max / maxval)
                )

    success = True
    for name, mismatch_count, ifirst, mismatch_max in sorted(
        mismatches, key=lambda x: x[2]
    ):
        msg = (
            "%s: mismatches at %i of %i times, first at time=%i with %i mismatches (max = %.3e, final max = %.3e)"
            % (
                name,
                (mismatch_count > 0).sum(),
                mismatch_count.size,
                ifirst,
                mismatch_count.flat[ifirst],
                mismatch_max.flat[ifirst],
                mismatch_max.flat[-1],
            )
        )
        if mismatch_max.flat[-1] < tolerance:
            msg += " - WITHIN TOLERANCE (%.3e)" % (tolerance,)
        else:
            success = False
        print(msg)
    return success


def compare_at_time(
    path1: str, path2: str, itime1: int, itime2: int, verbose: bool = False
) -> bool:
    success = True
    with netCDF4.Dataset(path1) as nc1, netCDF4.Dataset(path2) as nc2:
        for name in frozenset(nc1.variables) & frozenset(nc2.variables):
            ncvar1, ncvar2 = nc1.variables[name], nc2.variables[name]
            if "time" not in ncvar1.dimensions:
                continue
            values1 = np.ma.masked_invalid(ncvar1[itime1, ...])
            values2 = np.ma.masked_invalid(ncvar2[itime2, ...])
            mismatch = np.ma.filled(values1 != values2, False)
            if mismatch.any():
                success = False
                iworst = np.ma.filled(np.abs(values1 - values2), 0.0).ravel().argmax()
                print(
                    "%s: %i mismatches (out of %i)"
                    " - maximum absolute difference %s vs %s"
                    % (
                        name,
                        mismatch.sum(),
                        mismatch.size,
                        values1.flat[iworst],
                        values2.flat[iworst],
                    )
                )
            elif verbose:
                print("%s: identical" % name)
    return success


def compare_command():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("file1")
    parser.add_argument("file2")
    parser.add_argument("name_map", nargs="*")  # , action="append", default=[])
    parser.add_argument("--itime1", type=int)
    parser.add_argument("--itime2", type=int)
    parser.add_argument("--rotate", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_intermixed_args()
    if (args.itime1 is None and args.itime2 is not None) or (
        args.itime1 is not None and args.itime2 is None
    ):
        print(
            "ERROR: --itime1 and --itime2 must either be absent, or both be provided."
        )
        sys.exit(2)
    name_map = {}
    flip = set()
    for s in args.name_map:
        try:
            v1, v2 = s.split("=")
        except ValueError:
            raise Exception("Invalid assignment: %s" % s)
        if ".flip" in v2:
            flip.add(v1)
        name_map[v1] = v2.split(".", 1)[0]
    if args.itime1 is None:
        success = compare(
            args.file1, args.file2, args.rotate, name_map=name_map, flip=flip
        )
    else:
        success = compare_at_time(
            args.file1, args.file2, args.itime1, args.itime2, args.verbose
        )
    if not success:
        sys.exit(1)
