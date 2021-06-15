from typing import Optional
import functools

from mpi4py import MPI
import numpy

Waitall = MPI.Request.Waitall

class Tiling:
    def __init__(self, nrow: int, ncol: int, comm=MPI.COMM_WORLD, periodic_x: bool=False, periodic_y: bool=False):
        def find_neighbor(i, j):
            if periodic_x:
                j = j % ncol
            if periodic_y:
                i = i % nrow
            if i >= 0 and i < nrow and j >= 0 and j < ncol:
                return self.map[i, j]

        self.comm = comm
        assert self.comm.Get_size() == nrow * ncol, 'number of subdomains (%i rows * %i columns = %i) does not match group size of MPI communicator (%i)' % (nrow, ncol, nrow * ncol, self.comm.Get_size())
        self.rank = self.comm.Get_rank()
        self.map = numpy.arange(nrow * ncol).reshape(nrow, ncol)

        self.irow, self.icol = divmod(self.rank, ncol)
        self.nrow, self.ncol = nrow, ncol

        self.top = find_neighbor(self.irow + 1, self.icol)
        self.bottom = find_neighbor(self.irow - 1, self.icol)
        self.left = find_neighbor(self.irow, self.icol - 1)
        self.right = find_neighbor(self.irow, self.icol + 1)
        self.topleft = find_neighbor(self.irow + 1, self.icol - 1)
        self.topright = find_neighbor(self.irow + 1, self.icol + 1)
        self.bottomleft = find_neighbor(self.irow - 1, self.icol - 1)
        self.bottomright = find_neighbor(self.irow - 1, self.icol + 1)

    def wrap(self, field, halo: int):
        return DistributedArray(self, field, halo)

    def describe(self):
        p = lambda x: '-' if x is None else x
        print('{:^3}  {:^3}  {:^3}'.format(p(self.topleft), p(self.top), p(self.topright)))
        print('{:^3} [{:^3}] {:^3}'.format(p(self.left), self.rank, p(self.right)))
        print('{:^3}  {:^3}  {:^3}'.format(p(self.bottomleft), p(self.bottom), p(self.bottomright)))

    def plot(self, ax=None):
        x = numpy.linspace(0., self.ncol, self.ncol + 1)
        y = numpy.linspace(0., self.nrow, self.nrow + 1)
        if ax is None:
            import matplotlib.pyplot
            fig, ax = matplotlib.pyplot.subplots()
        shape = self.map.shape
        edge_shape = (shape[0] + 1, shape[1] + 1)
        ax.pcolormesh(numpy.broadcast_to(x[numpy.newaxis, :], edge_shape), numpy.broadcast_to(y[:, numpy.newaxis], edge_shape), numpy.ones_like(self.map), edgecolors='k')
        for i in range(self.nrow):
            for j in range(self.ncol):
                ax.text(0.5 * (x[j] + x[j + 1]), 0.5 * (y[i] + y[i + 1]), '%i' % self.map[i, j], horizontalalignment='center', verticalalignment='center')

class DistributedArray:
    __slots__ = ['tiling', 'comm', 'sendtasks', 'recvtasks', 'f_', 'f', 'halo']
    def __init__(self, tiling: Tiling, field: numpy.ndarray, halo: int):
        self.tiling = tiling
        self.comm = tiling.comm
        self.f_ = field
        self.f = field[..., halo:-halo, halo:-halo]
        self.halo = halo

        self.sendtasks = []
        self.recvtasks = []
        def add_task(neighbor, sendtag, recvtag, inner, outer):
            assert inner.shape == outer.shape
            if neighbor is not None:
                inner_cache = numpy.empty_like(inner)
                outer_cache = numpy.empty_like(outer)
                self.sendtasks.append((functools.partial(self.comm.Isend, inner_cache, neighbor, sendtag), inner, inner_cache))
                self.recvtasks.append((functools.partial(self.comm.Irecv, outer_cache, neighbor, recvtag), outer, outer_cache))

        add_task(tiling.left, 0, 1, field[..., halo:-halo, halo:halo*2], field[..., halo:-halo, :halo])
        add_task(tiling.right, 1, 0, field[..., halo:-halo, -halo*2:-halo], field[..., halo:-halo, -halo:])
        add_task(tiling.top, 2, 3, field[..., -halo*2:-halo, halo:-halo], field[..., -halo:, halo:-halo])
        add_task(tiling.bottom, 3, 2, field[..., halo:halo*2, halo:-halo], field[..., :halo, halo:-halo])
        add_task(tiling.topleft, 4, 7, field[..., -halo*2:-halo, halo:halo*2], field[..., -halo:, :halo])
        add_task(tiling.topright, 5, 6, field[..., -halo*2:-halo, -halo*2:-halo], field[..., -halo:, -halo:])
        add_task(tiling.bottomleft, 6, 5, field[..., halo:halo*2, halo:halo*2], field[..., :halo, :halo])
        add_task(tiling.bottomright, 7, 4, field[..., halo:halo*2, -halo*2:-halo], field[..., :halo, -halo:])

    def update_halos(self):
        recreqs = [fn() for fn, _, _ in self.recvtasks]
        sendreqs = []
        for fn, inner, cache in self.sendtasks:
            cache[...] = inner
            sendreqs.append(fn())
        Waitall(recreqs)
        for _, outer, cache in self.recvtasks:
            outer[...] = cache
        Waitall(sendreqs)

class Gather:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, root: int=0):
        self.rankmap = tiling.map
        self.comm = tiling.comm
        self.root = root
        self.field = field
        self.recvbuf = None
        if self.comm.Get_rank() == self.root:
            self.recvbuf = numpy.empty((self.rankmap.size,) + self.field.shape, dtype=self.field.dtype)

    def __call__(self, out: Optional[numpy.ndarray]=None) -> numpy.ndarray:
        sendbuf = numpy.ascontiguousarray(self.field)
        self.comm.Gather(sendbuf, self.recvbuf, root=self.root)
        if self.recvbuf is not None:
            nrow, ncol = self.rankmap.shape
            ny, nx = self.recvbuf.shape[-2:]
            if out is None:
                out = numpy.empty(self.recvbuf.shape[1:-2] + (nrow * ny, ncol * nx), dtype=self.recvbuf.dtype)
            for row in range(nrow):
                for col in range(ncol):
                    out[..., row * ny:(row + 1) * ny, col * nx:(col + 1) * nx] = self.recvbuf[self.rankmap[row, col], ...]
            return out

class Scatter:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, halo: int, share: int=0, root: int=0):
        self.field = field
        self.recvbuf = numpy.ascontiguousarray(field)
        self.rankmap = tiling.map
        self.halo = halo
        self.share = share
        self.sendbuf = None
        if tiling.comm.Get_rank() == root:
            self.sendbuf = numpy.zeros((self.rankmap.size,) + self.recvbuf.shape, dtype=self.recvbuf.dtype)
        self.mpi_scatter = functools.partial(tiling.comm.Scatter, self.sendbuf, self.recvbuf, root=root)

    def __call__(self, global_data: Optional[numpy.ndarray]):
        if self.sendbuf is not None:
            ny, nx = self.field.shape[-2:]
            xspacing, yspacing = nx - 2 * self.halo - self.share, ny - 2 * self.halo - self.share
            nrow, ncol = self.rankmap.shape
            assert nrow * yspacing + 2 * self.halo + self.share == global_data.shape[-2] and ncol * xspacing + 2 * self.halo + self.share == global_data.shape[-1], '%s, %i, %i' % (global_data.shape, nrow * yspacing + 2 * self.halo + self.share, ncol * xspacing + 2 * self.halo + self.share)
            for row in range(nrow):
                for col in range(ncol):
                    self.sendbuf[self.rankmap[row, col], ...] = global_data[..., row * yspacing:row * yspacing + ny, col * xspacing:col * xspacing + nx]
        self.mpi_scatter()
        if self.recvbuf is not self.field:
            self.field[...] = self.recvbuf

