from typing import Iterable, Mapping, MutableMapping, Optional, Tuple, Union, Any, List
import logging
import functools
import pickle
import enum

from mpi4py import MPI
import numpy
import numpy.typing

from . import _pygetm

Waitall = MPI.Request.Waitall
Startall = MPI.Prequest.Startall

def iterate_rankmap(rankmap):
    for irow in range(rankmap.shape[0]):
        for icol in range(rankmap.shape[1]):
            yield irow, icol, rankmap[irow, icol]

def getLogger(log_level=logging.INFO, comm=MPI.COMM_WORLD):
    rank = comm.Get_rank()
    ncpus = comm.Get_size()
    handlers = []
    if rank == 0:
        handlers.append(logging.StreamHandler())
    if ncpus > 1:
        file_handler = logging.FileHandler('getm-%04i.log' % rank, mode='w')
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(message)s'))
        handlers.append(file_handler)
    logging.basicConfig(level=log_level, handlers=handlers)

    return logging.getLogger()

class Tiling:
    @staticmethod
    def autodetect(mask, ncpus: Optional[int]=None, logger: Optional[logging.Logger]=None, max_protrude: float=0.5, **kwargs) -> 'Tiling':
        comm = kwargs.get('comm', MPI.COMM_WORLD)
        if comm.Get_rank() == 0:
            solution = find_optimal_divison(mask, ncpus, max_protrude=max_protrude, logger=logger)
        else:
            solution = None
        solution = comm.bcast(solution, root=0)
        counts = solution['map']
        rank_map = numpy.full(counts.shape, -1, dtype=int)
        rank = 0
        for irow, icol, count in iterate_rankmap(counts):
            if count > 0:
                rank_map[irow, icol] = rank
                rank += 1
        tiling = Tiling(map=rank_map, ncpus=solution['ncpus'], **kwargs)
        tiling.set_extent(mask.shape[1], mask.shape[0], solution['nx'], solution['ny'], yoffset_global=solution['yoffset'], xoffset_global=solution['xoffset'])
        return tiling

    @staticmethod
    def load(path: str) -> 'Tiling':
        with open(path, 'rb') as f:
            map, periodic_x, periodic_y, nx_glob, ny_glob, nx_sub, ny_sub, xoffset_global, yoffset_global = pickle.load(f)
        tiling = Tiling(map=map, periodic_x=periodic_x, periodic_y=periodic_y, ncpus=(map != -1).sum())
        tiling.set_extent(nx_glob, ny_glob, nx_sub, ny_sub, xoffset_global, yoffset_global)
        return tiling

    def __init__(self, nrow: Optional[int]=None, ncol: Optional[int]=None, map: Optional[numpy.ndarray]=None, comm=MPI.COMM_WORLD, periodic_x: bool=False, periodic_y: bool=False, ncpus: Optional[int]=None):
        if nrow is None or ncol is None:
            assert map is not None, 'If the number of rows and column in the subdomain decomposition is not provided, the rank map must be provided instead.'
            self.map = map
        else:
            self.map = numpy.arange(nrow * ncol).reshape(nrow, ncol)
        self.nrow, self.ncol = self.map.shape

        self.n_neigbors = 0
        def find_neighbor(i: int, j: int) -> int:
            if periodic_x:
                j = j % self.ncol
            if periodic_y:
                i = i % self.nrow
            if i >= 0 and i < self.nrow and j >= 0 and j < self.ncol:
                self.n_neigbors += 1
                return int(self.map[i, j])
            return -1

        self.comm = comm
        self.rank: int = self.comm.Get_rank()
        self.n = ncpus if ncpus is not None else self.comm.Get_size()
        nactive = (self.map[:, :] != -1).sum()
        assert nactive == self.n, 'number of active subdomains (%i) does not match group size of MPI communicator (%i). Map: %s' % (nactive, self.n, self.map)

        self.periodic_x = periodic_x
        self.periodic_y = periodic_y

        # Determine own row and column in subdomain decomposition
        for self.irow, self.icol, r in iterate_rankmap(self.map):
            if r == self.rank:
                break

        self.top = find_neighbor(self.irow + 1, self.icol)
        self.bottom = find_neighbor(self.irow - 1, self.icol)
        self.left = find_neighbor(self.irow, self.icol - 1)
        self.right = find_neighbor(self.irow, self.icol + 1)
        self.topleft = find_neighbor(self.irow + 1, self.icol - 1)
        self.topright = find_neighbor(self.irow + 1, self.icol + 1)
        self.bottomleft = find_neighbor(self.irow - 1, self.icol - 1)
        self.bottomright = find_neighbor(self.irow - 1, self.icol + 1)
        #print('Rank %i: top=%i, bottom=%i, left=%i, right=%i, topleft=%i, topright=%i, bottomleft=%i, bottomright=%i' % (self.rank, self.top, self.bottom, self.left, self.right, self.topleft, self.topright, self.bottomleft, self.bottomright))

        self.caches = {}
        self.nx_glob = None

    def dump(self, path: str):
        with open(path, 'wb') as f:
            pickle.dump((self.map, self.periodic_x, self.periodic_y, self.nx_glob, self.ny_glob, self.nx_sub, self.ny_sub, self.xoffset_global, self.yoffset_global), f)

    def set_extent(self, nx_glob: int, ny_glob: int, nx_sub: Optional[int]=None, ny_sub: Optional[int]=None, xoffset_global: int=0, yoffset_global: int=0):
        if nx_sub is None:
            nx_sub = int(numpy.ceil(nx_glob / self.ncol))
        if ny_sub is None:
            ny_sub = int(numpy.ceil(ny_glob / self.nrow))

        assert self.nx_glob is None, 'Domain extent has already been set.'
        assert isinstance(nx_glob, (int, numpy.integer))
        assert isinstance(ny_glob, (int, numpy.integer))
        assert isinstance(nx_sub, (int, numpy.integer))
        assert isinstance(ny_sub, (int, numpy.integer))
        assert isinstance(xoffset_global, (int, numpy.integer))
        assert isinstance(yoffset_global, (int, numpy.integer))

        self.nx_glob, self.ny_glob = nx_glob, ny_glob
        self.nx_sub, self.ny_sub = nx_sub, ny_sub
        self.xoffset_global = xoffset_global
        self.yoffset_global = yoffset_global
        self.xoffset = self.icol * self.nx_sub + self.xoffset_global
        self.yoffset = self.irow * self.ny_sub + self.yoffset_global

    def report(self, logger: Optional[logging.Logger]=None):
        if logger and (self.nrow > 1 or self.ncol > 1):
            logger.info('Using subdomain decomposition %i x %i (%i active nodes)' % (self.nrow, self.ncol, (self.map[:, :] != -1).sum()))
            logger.info('Global domain shape %i x %i, subdomain shape %i x %i, global offsets x=%i, y=%i' % (self.nx_glob, self.ny_glob, self.nx_sub, self.ny_sub, self.xoffset_global, self.yoffset_global))
            logger.info('I am rank %i at subdomain row %i, column %i, with offset x=%i, y=%i' % (self.rank, self.irow, self.icol, self.xoffset, self.yoffset))

    def __bool__(self) -> bool:
        return self.n_neigbors > 0

    def wrap(self, *args, **kwargs) -> 'DistributedArray':
        return DistributedArray(self, *args, **kwargs)

    def subdomain2slices(self, irow: Optional[int]=None, icol: Optional[int]=None, halo_sub: int=0, halo_glob: int=0, scale: int=1, share: int=0, exclude_halos: bool=True, exclude_global_halos: bool=False) -> Tuple[Tuple[Any, slice, slice], Tuple[Any, slice, slice], Tuple[int, int], Tuple[int, int]]:
        if irow is None:
            irow = self.irow
        if icol is None:
            icol = self.icol

        assert isinstance(share, int)
        assert isinstance(scale, int)
        assert isinstance(halo_sub, int)
        assert isinstance(halo_glob, int)
        assert self.map[irow, icol] >= 0

        # Global start and stop of local subdomain (no halos yet)
        xstart_glob = scale * ( icol      * self.nx_sub + self.xoffset_global)
        xstop_glob =  scale * ((icol + 1) * self.nx_sub + self.xoffset_global)
        ystart_glob = scale * ( irow      * self.ny_sub + self.yoffset_global)
        ystop_glob =  scale * ((irow + 1) * self.ny_sub + self.yoffset_global)

        # Adjust for halos and any overlap ("share")
        xstart_glob += - halo_sub + halo_glob
        xstop_glob  +=   halo_sub + halo_glob + share
        ystart_glob += - halo_sub + halo_glob
        ystop_glob  +=   halo_sub + halo_glob + share

        # Local start and stop
        xstart_loc = 0
        xstop_loc = xstop_glob - xstart_glob
        ystart_loc = 0
        ystop_loc = ystop_glob - ystart_glob

        # Shapes of local and global arrays (including halos and inactive strips)
        local_shape  = (scale * self.ny_sub  + 2 * halo_sub  + share, scale * self.nx_sub  + 2 * halo_sub  + share)
        global_shape = (scale * self.ny_glob + 2 * halo_glob + share, scale * self.nx_glob + 2 * halo_glob + share)

        # Calculate offsets based on limits of the global domain
        extra_offset = halo_sub if exclude_halos else 0
        margin_glob = halo_glob if exclude_global_halos else 0
        xstart_offset = max(xstart_glob + extra_offset, margin_glob                  ) - xstart_glob
        xstop_offset  = min(xstop_glob  - extra_offset, global_shape[1] - margin_glob) - xstop_glob
        ystart_offset = max(ystart_glob + extra_offset, margin_glob                  ) - ystart_glob
        ystop_offset  = min(ystop_glob  - extra_offset, global_shape[0] - margin_glob) - ystop_glob
        assert xstart_offset >= 0 and xstop_offset <= 0
        assert ystart_offset >= 0 and ystop_offset <= 0

        global_slice = (Ellipsis, slice(ystart_glob + ystart_offset, ystop_glob + ystop_offset), slice(xstart_glob + xstart_offset, xstop_glob + xstop_offset))
        local_slice  = (Ellipsis, slice(ystart_loc  + ystart_offset, ystop_loc  + ystop_offset), slice(xstart_loc  + xstart_offset, xstop_loc  + xstop_offset))
        return local_slice, global_slice, local_shape, global_shape

    def describe(self):
        p = lambda x: '-' if x is None else x
        print('{:^3}  {:^3}  {:^3}'.format(p(self.topleft), p(self.top), p(self.topright)))
        print('{:^3} [{:^3}] {:^3}'.format(p(self.left), self.rank, p(self.right)))
        print('{:^3}  {:^3}  {:^3}'.format(p(self.bottomleft), p(self.bottom), p(self.bottomright)))

    def plot(self, ax=None, background: Optional[numpy.ndarray]=None):
        x = self.xoffset_global + numpy.arange(self.ncol + 1) * self.nx_sub
        y = self.yoffset_global + numpy.arange(self.nrow + 1) * self.ny_sub
        if ax is None:
            import matplotlib.pyplot
            fig, ax = matplotlib.pyplot.subplots()
        import matplotlib.patches
        if background is not None:
            assert background.shape == (self.ny_glob, self.nx_glob), 'Argument background has incorrrect shape %s. Expected %s' % (background.shape, (self.ny_glob, self.nx_glob))
            ax.pcolormesh(numpy.arange(background.shape[1] + 1), numpy.arange(background.shape[0] + 1), background, alpha=0.5)
        else:
            ax.add_patch(matplotlib.patches.Rectangle((0, 0), self.nx_glob, self.ny_glob, edgecolor='None', facecolor='C0', zorder=-1))
        ax.pcolormesh(x, y, numpy.empty((self.nrow, self.ncol)), edgecolors='k', facecolor='none')
        for i in range(self.nrow):
            for j in range(self.ncol):
                if self.map[i, j] >= 0:
                    ax.text(0.5 * (x[j] + x[j + 1]), 0.5 * (y[i] + y[i + 1]), '%i' % self.map[i, j], horizontalalignment='center', verticalalignment='center')
                else:
                    ax.plot([x[j], x[j + 1]], [y[i], y[i + 1]], '-k')
                    ax.plot([x[j], x[j + 1]], [y[i + 1], y[i]], '-k')
        return ax

    def get_work_array(self, shape: Tuple[int, ...], dtype: numpy.typing.DTypeLike, fill_value=None):
        key = (shape, dtype, fill_value)
        if key not in self.caches:
            self.caches[key] = numpy.empty(shape, dtype)
            if fill_value is not None:
                self.caches[key].fill(fill_value)
        return self.caches[key]

@enum.unique
class Neighbor(enum.IntEnum):
    # Specific neighbors
    BOTTOMLEFT = 1
    BOTTOM = 2
    BOTTOMRIGHT = 3
    LEFT = 4
    RIGHT = 5
    TOPLEFT = 6
    TOP = 7
    TOPRIGHT = 8

    # Groups of neighbors (for update_halos command)
    ALL = 0
    TOP_AND_BOTTOM = 9
    LEFT_AND_RIGHT = 10
    TOP_AND_RIGHT = 11     # for T to U/V grid interpolation
    LEFT_AND_RIGHT_AND_TOP_AND_BOTTOM = 12  # for T to U/V grid 2nd order interpolation

GROUP2PARTS = {
    Neighbor.TOP_AND_BOTTOM: (Neighbor.TOP, Neighbor.BOTTOM),
    Neighbor.LEFT_AND_RIGHT: (Neighbor.LEFT, Neighbor.RIGHT),
    Neighbor.TOP_AND_RIGHT: (Neighbor.TOP, Neighbor.RIGHT),
    Neighbor.LEFT_AND_RIGHT_AND_TOP_AND_BOTTOM: (Neighbor.LEFT, Neighbor.RIGHT, Neighbor.TOP, Neighbor.BOTTOM),
}

class DistributedArray:
    __slots__ = ['rank', 'group2task', 'halo2name']
    def __init__(self, tiling: Tiling, field: numpy.ndarray, halo: int, overlap: int=0, share_caches: bool=False):
        self.rank = tiling.rank
        self.group2task: List[Tuple[List[MPI.Prequest], List[MPI.Prequest], List[Tuple[numpy.ndarray, numpy.ndarray]], List[Tuple[numpy.ndarray, numpy.ndarray]]]] = [([], [], [], []) for _ in range(max(Neighbor) + 1)]
        self.halo2name = {}

        key = (field.shape, halo, overlap, field.dtype)
        caches = None if not share_caches else tiling.caches.get(key)
        owncaches = []

        def add_task(name: str, recvtag: Neighbor, sendtag: Neighbor, outer: numpy.ndarray, inner: numpy.ndarray):
            neighbor = getattr(tiling, name)
            assert isinstance(neighbor, int), 'Wrong type for neighbor %s: %s (type %s)' % (name, neighbor, type(neighbor))
            assert inner.shape == outer.shape
            if neighbor != -1:
                if caches:
                    inner_cache, outer_cache = caches[len(owncaches)]
                    assert inner.shape == inner_cache.shape
                    assert outer.shape == outer_cache.shape
                else:
                    inner_cache, outer_cache = numpy.empty_like(inner), numpy.empty_like(outer)
                owncaches.append((inner_cache, outer_cache))
                send_req = tiling.comm.Send_init(inner_cache, neighbor, sendtag)
                recv_req = tiling.comm.Recv_init(outer_cache, neighbor, recvtag)
                self.halo2name[id(outer)] = name
                for group in [Neighbor.ALL, sendtag] + [group for (group, parts) in GROUP2PARTS.items() if sendtag in parts]:
                    send_reqs, recv_reqs, send_data, recv_data = self.group2task[group]
                    send_reqs.append(send_req)
                    send_data.append((inner, inner_cache))
                for group in [Neighbor.ALL, recvtag] + [group for (group, parts) in GROUP2PARTS.items() if recvtag in parts]:
                    send_reqs, recv_reqs, send_data, recv_data = self.group2task[group]
                    recv_reqs.append(recv_req)
                    recv_data.append((outer, outer_cache))

        in_start = halo + overlap
        in_stop = in_start + halo
        add_task('bottomleft',  Neighbor.BOTTOMLEFT,  Neighbor.TOPRIGHT,    field[...,      :halo,       :halo ], field[...,  in_start: in_stop,  in_start: in_stop ])
        add_task('bottom',      Neighbor.BOTTOM,      Neighbor.TOP,         field[...,      :halo,   halo:-halo], field[...,  in_start: in_stop,  halo    :-halo    ])
        add_task('bottomright', Neighbor.BOTTOMRIGHT, Neighbor.TOPLEFT,     field[...,      :halo,  -halo:     ], field[...,  in_start: in_stop, -in_stop :-in_start])
        add_task('left',        Neighbor.LEFT,        Neighbor.RIGHT,       field[...,  halo:-halo,      :halo ], field[...,  halo    :-halo,     in_start: in_stop ])
        add_task('right',       Neighbor.RIGHT,       Neighbor.LEFT,        field[...,  halo:-halo, -halo:     ], field[...,  halo    :-halo,    -in_stop :-in_start])
        add_task('topleft',     Neighbor.TOPLEFT,     Neighbor.BOTTOMRIGHT, field[..., -halo:,           :halo ], field[..., -in_stop :-in_start, in_start: in_stop ])
        add_task('top',         Neighbor.TOP,         Neighbor.BOTTOM,      field[..., -halo:,       halo:-halo], field[..., -in_stop :-in_start, halo    :-halo    ])
        add_task('topright',    Neighbor.TOPRIGHT,    Neighbor.BOTTOMLEFT,  field[..., -halo:,      -halo:     ], field[..., -in_stop :-in_start,-in_stop :-in_start])
        if caches is None and share_caches:
            tiling.caches[key] = owncaches

    def update_halos(self, group: Neighbor=Neighbor.ALL):
        send_reqs, recv_reqs, send_data, recv_data = self.group2task[group]
        Startall(recv_reqs)
        for inner, cache in send_data:
            cache[...] = inner
        Startall(send_reqs)
        Waitall(recv_reqs)
        for outer, cache in recv_data:
            outer[...] = cache
        Waitall(send_reqs)

    def update_halos_start(self, group: Neighbor=Neighbor.ALL):
        send_reqs, recv_reqs, send_data, _ = self.group2task[group]
        for inner, cache in send_data:
            cache[...] = inner
        Startall(send_reqs + recv_reqs)

    def update_halos_finish(self, group: Neighbor=Neighbor.ALL):
        send_reqs, recv_reqs, _, recv_data = self.group2task[group]
        Waitall(send_reqs + recv_reqs)
        for outer, cache in recv_data:
            outer[...] = cache

    def compare_halos(self, group: Neighbor=Neighbor.ALL) -> bool:
        send_reqs, recv_reqs, send_data, recv_data = self.group2task[group]
        Startall(recv_reqs)
        for inner, cache in send_data:
            cache[...] = inner
        Startall(send_reqs)
        Waitall(recv_reqs)
        match = True
        for outer, cache in recv_data:
            if not numpy.array_equal(outer, cache, equal_nan=True):
                delta = outer - cache
                print('Rank %i: mismatch in %s halo! Maximum absolute difference: %s. Current %s vs. received %s' % (self.rank, self.halo2name[id(outer)], numpy.abs(delta).max(), outer, cache))
                match = False
        Waitall(send_reqs)
        return match

class Sum:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, root: int=0):
        self.comm = tiling.comm
        self.root = root
        self.field = field
        self.result = None if tiling.rank != self.root else numpy.empty_like(field)

    def __call__(self) -> Optional[numpy.ndarray]:
        self.comm.Reduce(self.field, self.result, op=MPI.SUM, root=self.root)
        return self.result

class Gather:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, fill_value, root: int=0):
        self.comm = tiling.comm
        self.root = root
        self.field = field
        self.recvbuf = None
        self.fill_value = numpy.nan if fill_value is None else fill_value
        if tiling.rank == self.root:
            self.buffers = []
            self.recvbuf = tiling.get_work_array((tiling.n,) + self.field.shape, self.field.dtype)
            for irow, icol, rank in iterate_rankmap(tiling.map):
                if rank >= 0:
                    local_slice, global_slice, local_shape, global_shape = tiling.subdomain2slices(irow, icol, exclude_halos=True)
                    assert field.shape[-2:] == local_shape, 'Field shape %s differs from expected %s' % (field.shape[-2:], local_shape)
                    self.buffers.append((self.recvbuf[(rank,) + local_slice], global_slice))
            self.global_shape = global_shape

    def __call__(self, out: Optional[numpy.ndarray]=None, slice_spec=()) -> Optional[numpy.ndarray]:
        sendbuf = numpy.ascontiguousarray(self.field)
        self.comm.Gather(sendbuf, self.recvbuf, root=self.root)
        if self.recvbuf is not None:
            if out is None:
                out = numpy.empty(self.recvbuf.shape[1:-2] + self.global_shape, dtype=self.recvbuf.dtype)
                if self.fill_value is not None:
                    out.fill(self.fill_value)
            assert out.shape[-2:] == self.global_shape, 'Global shape %s differs from expected %s' % (out.shape[-2:], self.global_shape)
            for source, global_slice in self.buffers:
                out[slice_spec + global_slice] = source
            return out

class Scatter:
    def __init__(self, tiling: Tiling, field: numpy.ndarray, halo: int, share: int=0, fill_value=0., scale: int=1, exclude_halos: bool=False, root: int=0):
        self.field = field
        self.recvbuf = numpy.ascontiguousarray(field)
        self.rankmap = tiling.map
        self.halo = halo
        self.share = share
        self.sendbuf = None
        if tiling.comm.Get_rank() == root:
            self.buffers = []
            self.sendbuf = tiling.get_work_array((tiling.n,) + self.recvbuf.shape, self.recvbuf.dtype, fill_value=fill_value)
            for irow, icol, rank in iterate_rankmap(tiling.map):
                if rank >= 0:
                    local_slice, global_slice, local_shape, global_shape = tiling.subdomain2slices(irow, icol, halo_sub=halo, halo_glob=halo, share=share, scale=scale, exclude_halos=exclude_halos)
                    assert field.shape[-2:] == local_shape, 'Field shape %s differs from expected %s' % (field.shape[-2:], local_shape)
                    self.buffers.append((self.sendbuf[(rank,) + local_slice], global_slice))
            self.global_shape = global_shape
        self.mpi_scatter = functools.partial(tiling.comm.Scatter, self.sendbuf, self.recvbuf, root=root)

    def __call__(self, global_data: Optional[numpy.ndarray]):
        if self.sendbuf is not None:
            # we are root and have to send the global field
            assert global_data.shape[-2:] == self.global_shape, 'Global shape %s differs from expected %s' % (global_data.shape[-2:], self.global_shape)
            for sendbuf, global_slice in self.buffers:
                sendbuf[...] = global_data[global_slice]
        self.mpi_scatter()
        if self.recvbuf is not self.field:
            self.field[...] = self.recvbuf

def find_optimal_divison(mask: numpy.typing.ArrayLike, ncpus: Optional[int]=None, max_aspect_ratio: int=2, weight_unmasked: int=2, weight_any: int=1, weight_halo: int=10, max_protrude: float=0.5, logger: Optional[logging.Logger]=None) -> Optional[Mapping[str, Any]]:
    if ncpus is None:
        ncpus = MPI.COMM_WORLD.Get_size()
    if ncpus == 1:
        return {'ncpus': 1, 'nx': mask.shape[1], 'ny': mask.shape[0], 'xoffset': 0, 'yoffset': 0, 'cost': 0, 'map': numpy.ones((1,1), dtype=numpy.intc)}
    cost, solution = None, None
    mask = numpy.asarray(mask)
    mask_x, = mask.any(axis=0).nonzero()
    mask_y, = mask.any(axis=1).nonzero()
    imin, imax = mask_x[0], mask_x[-1]
    jmin, jmax = mask_y[0], mask_y[-1]
    mask = numpy.ascontiguousarray(mask[jmin:jmax + 1, imin:imax + 1], dtype=numpy.intc)
    if logger:
        logger.info('Determining optimal subdomain decomposition of global domain of %i x %i (%i active cells) for %i cores' % (mask.shape[1], mask.shape[0], (mask != 0).sum(), ncpus))
    for ny_sub in range(4, mask.shape[0] + 1):
        if logger and (ny_sub - 3) % 10 == 0:
            logger.info('%.1f %% complete' % (100 * (ny_sub - 4) / (mask.shape[0] + 1 - 4),))
        for nx_sub in range(max(4, ny_sub // max_aspect_ratio), min(max_aspect_ratio * ny_sub, mask.shape[1] + 1)):
            current_solution = _pygetm.find_subdiv_solutions(mask, nx_sub, ny_sub, ncpus, weight_unmasked, weight_any, weight_halo, max_protrude)
            if current_solution:
                xoffset, yoffset, current_cost, submap = current_solution
                if cost is None or current_cost < cost:
                    solution = {'ncpus': ncpus, 'nx': nx_sub, 'ny': ny_sub, 'xoffset': imin + xoffset, 'yoffset': jmin + yoffset, 'cost': current_cost, 'map': submap}
                    cost = current_cost
    if logger:
        logger.info('Optimal subdomain decomposition: %s' % solution)
    return solution

def find_all_optimal_divisons(mask: numpy.typing.ArrayLike, max_ncpus: Union[int, Iterable[int]], **kwargs) -> Optional[Mapping[str, Any]]:
    if isinstance(max_ncpus, int):
        max_ncpus = numpy.arange(1, max_ncpus + 1)
    for ncpus in max_ncpus:
        solution = find_optimal_divison(mask, ncpus, **kwargs)
        if solution:
            print('{ncpus} cores, optimal subdomain size: {nx} x {ny}, {0} land-only subdomains were dropped, max cost per core: {cost}'.format((solution['map'] == 0).sum(), **solution))

def test_scaling_command():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('setup_script', help='path to Python script that starts the simulation')
    parser.add_argument('--nmin', type=int, default=1, help='minimum number of CPUs to test with')
    parser.add_argument('--nmax', type=int, help='maximum number of CPUs to test with')
    parser.add_argument('--plot', action='store_true', help='whether to create scaling plot')
    parser.add_argument('--out', help='path to write scaling figure to')
    parser.add_argument('--compare', action='append', default=[], help='Path of NetCDF file to compare across simulations')
    args, leftover_args = parser.parse_known_args()
    if leftover_args:
        print('Arguments %s will be passed to %s' % (leftover_args, args.setup_script))
    test_scaling(args.setup_script, args.nmax, args.nmin, extra_args=leftover_args, plot=args.plot, out=args.out, compare=args.compare)

def test_scaling(setup_script: str, nmax: Optional[int]=None, nmin: int=1, extra_args: Iterable[str]=[], plot: bool=True, out: Optional[str]=None, compare: Optional[Iterable[str]]=None):
    import subprocess
    import sys
    import os
    import re
    import tempfile
    import shutil
    from .util import compare_nc
    if nmax is None:
        nmax = os.cpu_count()
    ncpus = []
    durations = []
    success = True
    refnames = []
    ntest = list(range(nmin, nmax + 1))
    if nmin > 1 and compare:
        ntest.insert(0, 1)
    for n in ntest:
        print('Running %s with %i CPUs... ' % (os.path.basename(setup_script), n), end='', flush=True)
        p = subprocess.run(['mpiexec', '-n', str(n), sys.executable, setup_script] + extra_args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.returncode != 0:
            success = False
            if n == ntest[0]:
                print('First simulation failed with return code %i - quitting. Last result:' % p.returncode)
                print(p.stdout)
                print()
                break
            log_path = 'scaling-%03i.log' % n
            with open(log_path, 'w') as f:
                f.write(p.stdout)
            print('FAILED with return code %i, log written to %s' % (p.returncode, log_path))
            continue
        m = re.search('Time spent in main loop: ([\d\.]+) s', p.stdout)
        duration = float(m.group(1))
        print('%.3f s in main loop' % (duration,))
        ncpus.append(n)
        durations.append(duration)
        for i, path in enumerate(compare):
            if n == 1:
                print('  copying reference %s to temporary file... ' % path, end='', flush=True)
                with open(path, 'rb') as fin, tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as fout:
                    shutil.copyfileobj(fin, fout)
                    refnames.append(fout.name)
                print('done.')
            else:
                print('  comparing %s wih reference produced with 1 CPUs... ' % (path,), flush=True)
                if not compare_nc.compare(path, refnames[i]):
                    success = False

    for refname in refnames:
        print('Deleting reference result %s... ' % refname, end='')
        os.remove(refname)
        print('done.')

    if plot and success:
        from matplotlib import pyplot
        fig, ax = pyplot.subplots()
        ax.plot(ncpus, durations, '.')
        ax.grid()
        ax.set_xlabel('number of CPUs')
        ax.set_ylabel('duration of %s' % os.path.basename(setup_script))
        ax.set_ylim(0, None)
        if out:
            fig.savefig(out, dpi=300)
        else:
            pyplot.show()

    sys.exit(0 if success else 1)