from mpi4py import MPI
import numpy

class ExchangedArray:
    def __init__(self, tiling, field, halo):
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

    def exchange(self):
        Irecv, Isend = self.comm.Irecv, self.comm.Isend
        reqs = []
        for neigbor, inner, outer, cache in self.tasks:
            reqs.append(Irecv(cache, neigbor))
        for neigbor, inner, outer, cache in self.tasks:
            reqs.append(Isend(inner.reshape(-1, order='A'), neigbor))
        MPI.Request.Waitall(reqs)
        for neigbor, inner, outer, cache in self.tasks:
            outer[...] = cache

class Tiling:
    def __init__(self, nx, ny):
        def get_rank(i, j):
            if i >= 0 and i < nx and j >= 0 and j < ny:
                return i * ny + j

        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
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
        return ExchangedArray(self, field, halo)