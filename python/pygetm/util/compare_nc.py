import sys
import numpy as np
import netCDF4


def compare(path1: str, path2: str) -> bool:
    mismatches = []
    with netCDF4.Dataset(path1) as nc1, netCDF4.Dataset(path2) as nc2:
        for name in frozenset(nc1.variables) & frozenset(nc2.variables):
            ncvar1, ncvar2 = nc1.variables[name], nc2.variables[name]
            if "time" not in ncvar1.dimensions:
                continue

            # Quick check for very last time point; if values match,
            # values at all preceding time points are assumed to match too.
            # Then skip the full values check (expensive!)
            values1 = np.ma.masked_invalid(ncvar1[-1, ...])
            values2 = np.ma.masked_invalid(ncvar2[-1, ...])
            if np.ma.filled(values1 == values2, True).all():
                continue

            # Values at last time did not match. Now check all time points
            values1 = np.ma.masked_invalid(ncvar1[...])
            values2 = np.ma.masked_invalid(ncvar2[...])
            mismatch = np.ma.filled(values1 != values2, False)
            mismatch_count = mismatch.sum(axis=tuple(range(1, mismatch.ndim)))
            if mismatch_count.sum() > 0:
                ifirst_mismatch = mismatch_count.nonzero()[0][0]
                mismatches.append((name, mismatch_count, ifirst_mismatch))

    for name, mismatch_count, ifirst in sorted(mismatches, key=lambda x: x[2]):
        print(
            "%s: mismatches at %i of %i times, first at time=%i with %i mismatches"
            % (
                name,
                (mismatch_count > 0).sum(),
                mismatch_count.size,
                ifirst,
                mismatch_count[ifirst],
            )
        )
    return len(mismatches) == 0


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
    parser.add_argument("--itime1", type=int)
    parser.add_argument("--itime2", type=int)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    if (args.itime1 is None and args.itime2 is not None) or (
        args.itime1 is not None and args.itime2 is None
    ):
        print(
            "ERROR: --itime1 and --itime2 must either be absent, or both be provided."
        )
        sys.exit(2)
    if args.itime1 is None:
        success = compare(args.file1, args.file2)
    else:
        success = compare_at_time(
            args.file1, args.file2, args.itime1, args.itime2, args.verbose
        )
    if not success:
        sys.exit(1)
