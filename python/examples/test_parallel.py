import pygetm.parallel

import numpy

halo = 2
plusint = 1
t = pygetm.parallel.Tiling(2, 2)
f = numpy.empty((3 * halo + plusint, 3 * halo + plusint))
f[...] = -1
f[..., halo:-halo, halo:-halo] = t.rank
exf = t.setup(f, halo)
exf.exchange()
print(t.rank, f)