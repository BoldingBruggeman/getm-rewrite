from typing import Optional
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
                self.sendtasks.append((neighbor, sendtag, inner, numpy.empty_like(inner)))
                self.recvtasks.append((neighbor, recvtag, outer, numpy.empty_like(outer)))

        add_task(tiling.left, 0, 1, field[..., halo:-halo, halo:halo*2], field[..., halo:-halo, :halo])
        add_task(tiling.right, 1, 0, field[..., halo:-halo, -halo*2:-halo], field[..., halo:-halo, -halo:])
        add_task(tiling.top, 2, 3, field[..., -halo*2:-halo, halo:-halo], field[..., -halo:, halo:-halo])
        add_task(tiling.bottom, 3, 2, field[..., halo:halo*2, halo:-halo], field[..., :halo, halo:-halo])
        add_task(tiling.topleft, 4, 7, field[..., -halo*2:-halo, halo:halo*2], field[..., -halo:, :halo])
        add_task(tiling.topright, 5, 6, field[..., -halo*2:-halo, -halo*2:-halo], field[..., -halo:, -halo:])
        add_task(tiling.bottomleft, 6, 5, field[..., halo:halo*2, halo:halo*2], field[..., :halo, :halo])
        add_task(tiling.bottomright, 7, 4, field[..., halo:halo*2, -halo*2:-halo], field[..., :halo, -halo:])

    def update_halos(self):
        recreqs = [self.comm.Irecv(cache, neigbor, tag) for neigbor, tag, _, cache in self.recvtasks]
        sendreqs = []
        for neigbor, tag, inner, cache in self.sendtasks:
            cache[...] = inner
            sendreqs.append(self.comm.Isend(cache, neigbor, tag))
        Waitall(recreqs)
        for _, _, outer, cache in self.recvtasks:
            outer[...] = cache
        Waitall(sendreqs)

    def gather(self, root: int=0, out: Optional[numpy.ndarray]=None):
        rankmap = self.tiling.map
        rank = self.comm.Get_rank()
        sendbuf = numpy.ascontiguousarray(self.f)
        recvbuf = None
        if rank == root:
            recvbuf = numpy.empty((rankmap.size,) + self.f.shape, dtype=self.f.dtype)
        self.comm.Gather(sendbuf, recvbuf, root=root)
        if rank == root:
            nrow, ncol = rankmap.shape
            ny, nx = recvbuf.shape[-2:]
            if out is None:
                out = numpy.empty(recvbuf.shape[1:-2] + (nrow * ny, ncol * nx), dtype=recvbuf.dtype)
            for i in range(nrow):
                for j in range(ncol):
                    out[..., i * ny:(i + 1) * ny, j * nx:(j + 1) * nx] = recvbuf[rankmap[i, j], ...]
            return out

    def scatter(self, data: numpy.ndarray, root: int=0):
        rankmap = self.tiling.map
        rank = self.comm.Get_rank()
        recvbuf = numpy.empty_like(self.f_)
        sendbuf = None
        if rank == root:
            sendbuf = numpy.zeros((rankmap.size,) + recvbuf.shape, dtype=recvbuf.dtype)
            ny, nx = self.f.shape[-2:]
            nrow, ncol = rankmap.shape
            halo = self.halo
            assert nrow * ny == data.shape[-2] and ncol * nx == data.shape[-1], '%s, %i, %i' % (data.shape, nrow * ny, ncol * nx)
            for i in range(nrow):
                for j in range(ncol):
                    imin_off = 0 if i == 0 else halo
                    imax_off = 0 if i == nrow - 1 else halo
                    jmin_off = 0 if j == 0 else halo
                    jmax_off = 0 if j == ncol - 1 else halo
                    sendbuf[rankmap[i, j], ..., halo - imin_off:halo + ny + imax_off, halo - jmin_off:halo + nx + jmax_off] = data[..., i * ny - imin_off:(i + 1) * ny + imax_off, j * nx - jmin_off:(j + 1) * nx + jmax_off]
        self.comm.Scatter(sendbuf, recvbuf, root=root)
        self.f_[...] = recvbuf

