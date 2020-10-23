from mpi4py import MPI
import numpy

class Tiling:
    def __init__(self, nx, ny, comm=MPI.COMM_WORLD):
        def get_rank(i, j):
            if i >= 0 and i < nx and j >= 0 and j < ny:
                return self.map[i, j]

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.map = numpy.arange(nx * ny).reshape(nx, ny)

        irow, icol = divmod(self.rank, nx)

        self.top = get_rank(irow - 1, icol)
        self.bottom = get_rank(irow + 1, icol)
        self.left = get_rank(irow, icol - 1)
        self.right = get_rank(irow, icol + 1)
        self.topleft = get_rank(irow - 1, icol - 1)
        self.topright = get_rank(irow - 1, icol + 1)
        self.bottomleft = get_rank(irow + 1, icol - 1)
        self.bottomright = get_rank(irow + 1, icol + 1)

    def setup(self, field, halo):
        return DistributedArray(self, field, halo)

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
        add_task(tiling.top, field[..., halo:halo*2, halo:-halo], field[..., :halo, halo:-halo])
        add_task(tiling.bottom, field[..., -halo*2:-halo, halo:-halo], field[..., -halo:, halo:-halo])
        add_task(tiling.topleft, field[..., halo:halo*2, halo:halo*2], field[..., :halo, :halo])
        add_task(tiling.topright, field[..., halo:halo*2, -halo*2:-halo], field[..., :halo, -halo:])
        add_task(tiling.bottomleft, field[..., -halo*2:-halo, halo:halo*2], field[..., -halo:, :halo])
        add_task(tiling.bottomright, field[..., -halo*2:-halo, -halo*2:-halo], field[..., -halo:, -halo:])
        self.interior = field[..., halo:-halo, halo:-halo]

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
            nx, ny = rankmap.shape
            ni, nj = recvbuf.shape[-2:]
            out = numpy.empty(recvbuf.shape[1:-2] + (nx * ni, ny * nj), dtype=recvbuf.dtype)
            for i in range(nx):
                for j in range(ny):
                    out[..., i * ni:(i + 1) * ni, j * nj:(j + 1) * nj] = recvbuf[rankmap[i, j], ...]
            return out

