import sys
import numpy
import netCDF4

def compare(path1: str, path2: str) -> bool:
    mismatches = []
    with netCDF4.Dataset(path1) as nc1, netCDF4.Dataset(path2) as nc2:
        for name in frozenset(nc1.variables) & frozenset(nc2.variables):
            ncvar1, ncvar2 = nc1.variables[name], nc2.variables[name]
            if 'time' not in ncvar1.dimensions:
                continue
            if numpy.ma.filled(numpy.ma.masked_invalid(ncvar1[-1, ...]) == numpy.ma.masked_invalid(ncvar2[-1, ...]), True).all():
                continue
            values1 = numpy.ma.masked_invalid(ncvar1[...])
            values2 = numpy.ma.masked_invalid(ncvar2[...])
            mismatch = numpy.ma.filled(values1 != values2, False)
            mismatch_count = mismatch.sum(axis=tuple(range(1, mismatch.ndim)))
            if mismatch_count.sum() > 0:
                mismatches.append((name, mismatch_count, mismatch_count.nonzero()[0][0]))
    for name, mismatch_count, ifirst in sorted(mismatches, key=lambda x: x[2]):
        print('%s: mismatches at %i of %i times, first at time=%i with %i mismatches' % (name, (mismatch_count > 0).sum(), mismatch_count.size, ifirst, mismatch_count[ifirst]))
    return len(mismatches) == 0

def compare_command():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file1')
    parser.add_argument('file2')    
    args = parser.parse_args()
    if not compare(args.file1, args.file2):
        sys.exit(1)