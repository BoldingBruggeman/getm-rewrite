from typing import Optional, Tuple
import functools

from mpi4py import MPI
import numpy

Waitall = MPI.Request.Waitall

class Tiling:
    def __init__(self, nrow: Optional[int]=None, ncol: Optional[int]=None, comm=MPI.COMM_WORLD, periodic_x: bool=False, periodic_y: bool=False):
        self.n_neigbors = 0
        def find_neighbor(i, j):
            if periodic_x:
                j = j % ncol
            if periodic_y:
                i = i % nrow
            if i >= 0 and i < nrow and j >= 0 and j < ncol:
                self.n_neigbors += 1
                return self.map[i, j]

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.n: int = self.comm.Get_size()
        if nrow is None and ncol is None:
            # Temporary algorithm for determining subdomain decomposition
            best = None
            for nrow in range(1, self.n + 1):
                if self.n % nrow == 0:
                    ncol = self.n // nrow
                    n_exchange = nrow * ncol * 8 - nrow * 2 - ncol * 2 - 4
                    if best is None or n_exchange < best[0]:
                        best = (n_exchange, nrow, ncol)
            nrow, ncol = best[1:]
            if self.n > 1 and self.rank == 0:
                print('Using subdomain decomposition %i x %i' % (nrow, ncol))
        assert nrow * ncol in (1, self.n), 'number of subdomains (%i rows * %i columns = %i) does not match group size of MPI communicator (%i)' % (nrow, ncol, nrow * ncol, self.n)
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

        self.caches = {}

    def __bool__(self) -> bool:
        return self.n_neigbors > 0

    def wrap(self, field: numpy.ndarray, halo: int):
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

ALL = 0
TOP_BOTTOM = 1
LEFT_RIGHT = 2

class DistributedArray:
    __slots__ = ['rank', 'group2task', 'neighbor2name']
    def __init__(self, tiling: Tiling, field: numpy.ndarray, halo: int):
        self.rank = tiling.rank
        self.group2task = {ALL: ([], []), TOP_BOTTOM: ([], []), LEFT_RIGHT: ([], [])}
        self.neighbor2name = {}

        key = (field.shape, halo, field.dtype)
        caches = tiling.caches.get(key)
        owncaches = []

        def add_task(name: str, sendtag: int, recvtag: int, inner: numpy.ndarray, outer: numpy.ndarray, groups: Tuple[int, ...]):
            neighbor = getattr(tiling, name)
            assert inner.shape == outer.shape
            if neighbor is not None:
                if caches:
                    inner_cache, outer_cache = caches[len(owncaches)]
                else:
                    inner_cache, outer_cache = numpy.empty_like(inner), numpy.empty_like(outer)
                owncaches.append((inner_cache, outer_cache))
                sendtask = (functools.partial(tiling.comm.Isend, inner_cache, neighbor, sendtag), inner, inner_cache)
                recvtask = (functools.partial(tiling.comm.Irecv, outer_cache, neighbor, recvtag), outer, outer_cache)
                self.neighbor2name[neighbor] = name
                for group in groups:
                    sendtasks, recvtasks = self.group2task[group]
                    sendtasks.append(sendtask)
                    recvtasks.append(recvtask)

        add_task('bottomleft',  6, 5, field[...,  halo  :halo*2,  halo  :halo*2], field[...,      :halo,       :halo ], groups=(ALL,))
        add_task('bottom',      3, 2, field[...,  halo  :halo*2,  halo  :-halo],  field[...,      :halo,   halo:-halo], groups=(ALL, TOP_BOTTOM))
        add_task('bottomright', 7, 4, field[...,  halo  :halo*2, -halo*2:-halo],  field[...,      :halo,  -halo:     ], groups=(ALL,))
        add_task('left',        0, 1, field[...,  halo  :-halo,   halo  :halo*2], field[...,  halo:-halo,      :halo ], groups=(ALL, LEFT_RIGHT))
        add_task('right',       1, 0, field[...,  halo  :-halo,  -halo*2:-halo],  field[...,  halo:-halo, -halo:     ], groups=(ALL, LEFT_RIGHT))
        add_task('topleft',     4, 7, field[..., -halo*2:-halo,   halo  :halo*2], field[..., -halo:,           :halo ], groups=(ALL,))
        add_task('top',         2, 3, field[..., -halo*2:-halo,   halo  :-halo],  field[..., -halo:,       halo:-halo], groups=(ALL, TOP_BOTTOM))
        add_task('topright',    5, 6, field[..., -halo*2:-halo,  -halo*2:-halo],  field[..., -halo:,      -halo:     ], groups=(ALL,))
        if caches is None:
            tiling.caches[key] = owncaches

    def update_halos(self, group: int=ALL):
        sendtasks, recvtasks = self.group2task[group]
        recreqs = [fn() for fn, _, _ in recvtasks]
        sendreqs = []
        for fn, inner, cache in sendtasks:
            cache[...] = inner
            sendreqs.append(fn())
        Waitall(recreqs)
        for _, outer, cache in recvtasks:
            outer[...] = cache
        Waitall(sendreqs)

    def compare_halos(self, group: int=ALL) -> bool:
        sendtasks, recvtasks = self.group2task[group]
        recreqs = [fn() for fn, _, _ in recvtasks]
        sendreqs = []
        for fn, inner, cache in sendtasks:
            cache[...] = inner
            sendreqs.append(fn())
        Waitall(recreqs)
        match = True
        for fn, outer, cache in recvtasks:
            if not (outer == cache).all():
                delta = outer - cache
                print('Rank %i: mismatch in %s halo! Maximum absolute difference: %s. Values: %s' % (self.rank, self.neighbor2name[fn.args[1]], numpy.abs(delta).max(), delta))
                match = False
        Waitall(sendreqs)
        return match

class Sum:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, root: int=0):
        self.comm = tiling.comm
        self.root = root
        self.field = field
        self.result = None if tiling.rank != self.root else numpy.empty_like(field)

    def __call__(self):
        self.comm.Reduce(self.field, self.result, op=MPI.SUM, root=self.root)
        return self.result

class Gather:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, root: int=0):
        self.rankmap = tiling.map
        self.comm = tiling.comm
        self.root = root
        self.field = field
        self.recvbuf = None
        if tiling.rank == self.root:
            self.recvbuf = numpy.empty((self.rankmap.size,) + self.field.shape, dtype=self.field.dtype)

    def __call__(self, out: Optional[numpy.ndarray]=None, slice_spec=()) -> numpy.ndarray:
        sendbuf = numpy.ascontiguousarray(self.field)
        self.comm.Gather(sendbuf, self.recvbuf, root=self.root)
        if self.recvbuf is not None:
            nrow, ncol = self.rankmap.shape
            ny, nx = self.recvbuf.shape[-2:]
            if out is None:
                out = numpy.empty(self.recvbuf.shape[1:-2] + (nrow * ny, ncol * nx), dtype=self.recvbuf.dtype)
            for row in range(nrow):
                for col in range(ncol):
                    s = slice_spec + (Ellipsis, slice(row * ny, (row + 1) * ny), slice(col * nx, (col + 1) * nx))
                    out[s] = self.recvbuf[self.rankmap[row, col], ...]
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

