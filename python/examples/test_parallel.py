import numpy
from mpi4py import MPI
import pygetm.parallel

size = MPI.COMM_WORLD.Get_size()
assert size == 4, 'This example must be run on 4 nodes (found %i)' % size

halo = 1
plusint = 2
t = pygetm.parallel.Tiling(2, 2)
f = numpy.empty((2, 3 * halo + plusint, 3 * halo + plusint))
f[...] = -1
f[..., halo:-halo, halo:-halo] = t.rank
exf = t.wrap(f, halo)
exf.update_halos()
print(t.rank, f)

print(t.rank, exf.gather())