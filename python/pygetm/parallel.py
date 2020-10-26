from mpi4py import MPI
import numpy

class Tiling:
    def __init__(self, nrow, ncol, comm=MPI.COMM_WORLD):
        def get_rank(i, j):
            if i >= 0 and i < nrow and j >= 0 and j < ncol:
                return self.map[i, j]

        self.comm = comm
        assert self.comm.Get_size() == nrow * ncol
        self.rank = self.comm.Get_rank()
        self.map = numpy.arange(nrow * ncol).reshape(nrow, ncol)

        self.irow, self.icol = divmod(self.rank, ncol)
        self.nrow, self.ncol = nrow, ncol

        self.top = get_rank(self.irow + 1, self.icol)
        self.bottom = get_rank(self.irow - 1, self.icol)
        self.left = get_rank(self.irow, self.icol - 1)
        self.right = get_rank(self.irow, self.icol + 1)
        self.topleft = get_rank(self.irow + 1, self.icol - 1)
        self.topright = get_rank(self.irow + 1, self.icol + 1)
        self.bottomleft = get_rank(self.irow - 1, self.icol - 1)
        self.bottomright = get_rank(self.irow - 1, self.icol + 1)

    def wrap(self, field, halo):
        return DistributedArray(self, field, halo)

    def describe(self):
        p = lambda x: '-' if x is None else x
        print('{:^3}  {:^3}  {:^3}'.format(p(self.topleft), p(self.top), p(self.topright)))
        print('{:^3} [{:^3}] {:^3}'.format(p(self.left), self.rank, p(self.right)))
        print('{:^3}  {:^3}  {:^3}'.format(p(self.bottomleft), p(self.bottom), p(self.bottomright)))

class DistributedArray:
    def __init__(self, tiling, field, halo):
        self.tiling = tiling
        self.comm = tiling.comm

        self.tasks = []
        def add_task(neighbor, inner, outer):
            assert inner.shape == outer.shape
            if neighbor is not None:
                self.tasks.append((neighbor, inner, outer, numpy.empty_like(outer)))

        add_task(tiling.left, field[..., halo:-halo, halo:halo*2], field[..., halo:-halo, :halo])
        add_task(tiling.right, field[..., halo:-halo, -halo*2:-halo], field[..., halo:-halo, -halo:])
        add_task(tiling.top, field[..., -halo*2:-halo, halo:-halo], field[..., -halo:, halo:-halo])
        add_task(tiling.bottom, field[..., halo:halo*2, halo:-halo], field[..., :halo, halo:-halo])
        add_task(tiling.topleft, field[..., -halo*2:-halo, halo:halo*2], field[..., -halo:, :halo])
        add_task(tiling.topright, field[..., -halo*2:-halo, -halo*2:-halo], field[..., -halo:, -halo:])
        add_task(tiling.bottomleft, field[..., halo:halo*2, halo:halo*2], field[..., :halo, :halo])
        add_task(tiling.bottomright, field[..., halo:halo*2, -halo*2:-halo], field[..., :halo, -halo:])
        self.field = field
        self.interior = field[..., halo:-halo, halo:-halo]
        self.halo = halo

    def update_halos(self):
        Irecv, Isend = self.comm.Irecv, self.comm.Isend
        reqs = []
        for neigbor, inner, outer, cache in self.tasks:
            reqs.append(Irecv(cache, neigbor))
        for neigbor, inner, outer, cache in self.tasks:
            reqs.append(Isend(numpy.ascontiguousarray(inner), neigbor))
        MPI.Request.Waitall(reqs)
        for neigbor, inner, outer, cache in self.tasks:
            outer[...] = cache

    def gather(self, root=0):
        rankmap = self.tiling.map
        rank = self.comm.Get_rank()
        sendbuf = numpy.ascontiguousarray(self.interior)
        recvbuf = None
        if rank == root:
            recvbuf = numpy.empty((rankmap.size,) + self.interior.shape, dtype=self.interior.dtype)
        self.comm.Gather(sendbuf, recvbuf, root=root)
        if rank == root:
            nrow, ncol = rankmap.shape
            ny, nx = recvbuf.shape[-2:]
            out = numpy.empty(recvbuf.shape[1:-2] + (nrow * ny, ncol * nx), dtype=recvbuf.dtype)
            for i in range(nrow):
                for j in range(ncol):
                    out[..., i * ny:(i + 1) * ny, j * nx:(j + 1) * nx] = recvbuf[rankmap[i, j], ...]
            return out

    def scatter(self, data, root=0):
        rankmap = self.tiling.map
        rank = self.comm.Get_rank()
        recvbuf = numpy.empty_like(self.field)
        sendbuf = None
        if rank == root:
            sendbuf = numpy.zeros((rankmap.size,) + recvbuf.shape, dtype=recvbuf.dtype)
            ny, nx = self.interior.shape[-2:]
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
        self.field[...] = recvbuf

