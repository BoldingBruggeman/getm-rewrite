# cython: language_level=3

cimport cython

from libc.math cimport ceil

cimport numpy
import numpy

cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz) nogil
cdef extern void domain_initialize_open_boundaries(void* domain, int nbdyp, int nwb, int nnb, int neb, int nsb, int* bdy_info, int* bdy_i, int* bdy_j) nogil
cdef extern void* domain_get_grid(void* domain, int grid_type) nogil
cdef extern void domain_initialize(void* grid, int runtype, double* maxdt) nogil
cdef extern void domain_finalize(void* domain) nogil
cdef extern void domain_update_depths(void* domain) nogil
cdef extern void domain_do_vertical(void* domain) nogil
cdef extern void grid_interp_x(void* grid, double* source, double* target, int ioffset) nogil
cdef extern void grid_interp_y(void* grid, double* source, double* target, int joffset) nogil
cdef extern void grid_interp_xy(void* source_grid, double* source, void* target_grid, double* target, int ioffset, int joffset) nogil
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* sub_type, int* data_type, void** p) nogil
cdef extern void* advection_create(int scheme, void* tgrid, void** p) nogil
cdef extern void advection_2d_calculate(int direction, void* advection, void* tgrid, void* ugrid, double* pu, double timestep, double* pvar) nogil
cdef extern void* vertical_diffusion_create(void* tgrid) nogil
cdef extern void vertical_diffusion_calculate(void* diffusion, void* tgrid, double molecular, double* pnuh, double timestep, double cnpar, double* pvar) nogil
cdef extern void* momentum_create(int runtype, void* pdomain, int apply_bottom_friction) nogil
cdef extern void momentum_uv_momentum_2d(void* momentum, int runtype, double timestep, double* ptausx, double* ptausy, double* pdpdx, double* pdpdy) nogil
cdef extern void* pressure_create(int runtype, void* pdomain) nogil
cdef extern void pressure_surface(void* pressure, double* pz, double* psp) nogil
cdef extern void* sealevel_create(void* pdomain) nogil
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV) nogil
#cdef extern void sealevel_update_uvx(void* sealevel) nogil
cdef extern void sealevel_boundaries(void* sealevel, void* momentum, double timestep) nogil

cpdef enum:
    TGRID = 1
    UGRID = 2
    VGRID = 3
    XGRID = 4
    WGRID = 5
    UUGRID = -1
    VVGRID = -2
    UVGRID = -3
    VUGRID = -4

cdef class Array:
    cdef void* p
    cdef readonly numpy.ndarray all_values
    cdef readonly Grid grid

    def __init__(self, Grid grid=None):
        if grid is not None:
            self.grid = grid

    cdef wrap_c_array(self, Domain domain, int source, void* obj, bytes name):
        cdef int grid_type
        cdef int sub_type
        cdef int data_type
        get_array(source, obj, name, &grid_type, &sub_type, &data_type, &self.p)
        if self.p == NULL:
            return
        self.grid = domain.grids[grid_type]
        if sub_type == 1:
            # Horizontal-only array on open boundaries
            if data_type == 0:
                self.all_values = numpy.asarray(<double[:self.grid.nbdyp:1]> self.p)
            else:
                self.all_values = numpy.asarray(<int[:self.grid.nbdyp:1]> self.p)
            return self.all_values
        elif sub_type == 2:
            # Depth-explicit array on normal grid
            if data_type == 0:
                self.all_values = numpy.asarray(<double[:self.grid.nz_, :self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self.all_values = numpy.asarray(<int[:self.grid.nz_, :self.grid.ny_, :self.grid.nx_:1]> self.p)
        else:
            # Horizontal-only array on normal grid
            assert sub_type == 0, 'Subtypes other than 0,1,2 not yet implemented'
            if data_type == 0:
                self.all_values = numpy.asarray(<double[:self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self.all_values = numpy.asarray(<int[:self.grid.ny_, :self.grid.nx_:1]> self.p)
        self.finish_initialization()
        self.register()
        return self

    def wrap_ndarray(self, numpy.ndarray data):
        assert data.ndim in (2, 3) and data.flags['C_CONTIGUOUS']
        self.all_values = data
        self.p = self.all_values.data
        self.finish_initialization()

cdef class Grid:
    cdef void* p
    cdef readonly int nx, ny, nz
    cdef readonly int nx_, ny_, nz_
    cdef readonly Domain domain

    def __init__(self, Domain domain, int grid_type):
        self.domain = domain
        self.p = domain_get_grid(domain.p, grid_type)
        self.nx, self.ny, self.nz = domain.nx, domain.ny, domain.nz
        if grid_type == XGRID:
            self.nx += 1
            self.ny += 1
        if grid_type == WGRID:
            self.nz += 1
        self.nx_, self.ny_, self.nz_ = self.nx + 2 * domain.halox, self.ny + 2 * domain.haloy, self.nz + 2 * domain.haloz
        domain.grids[grid_type] = self

    def wrap(self, Array ar, bytes name):
        return ar.wrap_c_array(self.domain, 0, self.p, name)

    def interp_x(self, Array source, Array target, int offset):
        grid_interp_x(self.p, <double *>source.p, <double *>target.p, offset)

    def interp_y(self, Array source, Array target, int offset):
        grid_interp_y(self.p, <double *>source.p, <double *>target.p, offset)

    def interp_xy(self, Array source, Array target, int ioffset, int joffset):
        grid_interp_xy(self.p, <double *>source.p, <double *>target.grid.p, <double *>target.p, ioffset, joffset)

cdef class Domain:
    cdef void* p
    cdef readonly int nx, ny, nz
    cdef readonly int halox, haloy, haloz
    cdef readonly double maxdt

    def __init__(self, int imin, int imax, int jmin, int jmax, int kmin, int kmax):
        self.p = domain_create(imin, imax, jmin, jmax, kmin, kmax, &self.halox, &self.haloy, &self.haloz)
        self.nx = imax - imin + 1
        self.ny = jmax - jmin + 1
        self.nz = kmax - kmin + 1
        self.grids = {}

    def __dealloc__(self):
        if self.p != NULL:
            domain_finalize(self.p)

    def update_depths(self):
        domain_update_depths(self.p)

    def do_vertical(self):
        domain_do_vertical(self.p)

    def initialize(self, int runtype):
        domain_initialize(self.p, runtype, &self.maxdt)

    def initialize_open_boundaries(self, int nwb, int nnb, int neb, int nsb, int nbdyp, int[:,::1] bdy_info, int[::1] bdy_i, int[::1] bdy_j):
        assert bdy_info.shape[0] == 6, 'bdy_info should have 6 rows'
        assert bdy_info.shape[1] == nwb + nnb + neb + nsb, 'bdy_info should have as many columns as the number of boundaries'
        assert bdy_i.shape[0] == nbdyp, 'bdy_i should have a length equal to the number of boundary points'
        assert bdy_j.shape[0] == nbdyp, 'bdy_j should have a length equal to the number of boundary points'
        domain_initialize_open_boundaries(self.p, nbdyp, nwb, nnb, neb, nsb, &bdy_info[0,0], &bdy_i[0], &bdy_j[0])

cdef class Advection:
    cdef void* p
    cdef Grid tgrid
    cdef Grid ugrid
    cdef Grid vgrid
    cdef readonly numpy.ndarray D

    def __init__(self, Grid grid, int scheme):
        cdef void* pD
        self.tgrid = grid
        self.ugrid = grid.ugrid
        self.vgrid = grid.vgrid
        self.p = advection_create(scheme, self.tgrid.p, &pD)
        self.D = numpy.asarray(<double[:self.tgrid.ny_, :self.tgrid.nx_:1]> pD)

    def calculate(self, Array u not None, Array v not None, double timestep, Array var not None):
        self.D[...] = self.tgrid.D.all_values
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, 0.5 * timestep, <double *>var.p)
        var.update_halos(2)
        advection_2d_calculate(1, self.p, self.tgrid.p, self.ugrid.p, <double *>u.p, timestep, <double *>var.p)
        var.update_halos(1)
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, 0.5 * timestep, <double *>var.p)

cdef class VerticalDiffusion:
    cdef void* p
    cdef Grid tgrid
    cdef double cnpar

    def __init__(self, Grid grid, double cnpar=1.):
        self.tgrid = grid
        self.p = vertical_diffusion_create(self.tgrid.p)
        self.cnpar = cnpar

    def calculate(self, Array nuh not None, double timestep, Array var not None, double molecular=0.):
        vertical_diffusion_calculate(self.p, self.tgrid.p, molecular, <double *>nuh.p, timestep, self.cnpar, <double *>var.p)

cdef class Simulation:
    cdef readonly Domain domain
    cdef readonly int runtype
    cdef void* pmomentum
    cdef void* ppressure
    cdef void* psealevel
    cdef readonly int nx, ny

    def __init__(self, Domain domain, int runtype, int apply_bottom_friction):
        self.domain = domain
        domain.initialize(runtype)
        self.runtype = runtype
        self.pmomentum = momentum_create(runtype, domain.p, apply_bottom_friction)
        self.ppressure = pressure_create(runtype, domain.p)
        self.psealevel = sealevel_create(domain.p)
        self.nx, self.ny = domain.nx, domain.ny

    def uv_momentum_2d(self, double timestep, Array tausx not None, Array tausy not None, Array dpdx not None, Array dpdy not None):
        momentum_uv_momentum_2d(self.pmomentum, self.runtype, timestep, <double *>tausx.p, <double *>tausy.p, <double *>dpdx.p, <double *>dpdy.p)

    def update_surface_pressure_gradient(self, Array z not None, Array sp not None):
        pressure_surface(self.ppressure, <double *>z.p, <double *>sp.p)

    def update_sealevel(self, double timestep, Array U not None, Array V not None):
        sealevel_update(self.psealevel, timestep, <double *>U.p, <double *>V.p)

    #def update_sealevel_uvx(self):
    #    sealevel_update_uvx(self.psealevel)

    def update_sealevel_boundaries(self, double timestep):
        sealevel_boundaries(self.psealevel, self.pmomentum, timestep)

    def wrap(self, Array ar not None, bytes name, int source):
        cdef void* obj = self.pmomentum
        if (source == 2):
            obj = self.ppressure
        elif (source == 3):
            obj = self.psealevel
        return ar.wrap_c_array(self.domain, source, obj, name)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int[:, ::1] get_map(const int[:, ::1] mask, int nx, int ny, int ioffset, int joffset, int ncpus):
    cdef int nrow = <int>ceil((mask.shape[0] - joffset) / float(ny))
    cdef int ncol = <int>ceil((mask.shape[1] - ioffset) / float(nx))
    cdef int[:, ::1] map
    cdef int current_ncpus = 0
    cdef int row, col, i, j, n
    if ncpus != -1 and nrow * ncol < ncpus:
        return None
    map = numpy.empty((nrow, ncol), dtype=int)
    for row in range(nrow):
        for col in range(ncol):
            n = 0
            for j in range(max(0, joffset + row * ny), min(mask.shape[0], joffset + (row + 1) * ny)):
                for i in range(max(0, ioffset + col * nx), min(mask.shape[1], ioffset + (col + 1) * nx)):
                    if mask[j, i] != 0: n += 1
            if n > 0:
                current_ncpus += 1
                if ncpus != -1 and current_ncpus > ncpus:
                    return None
            map[row, col] = n
    if ncpus == -1 or current_ncpus == ncpus:
        return map
    return None

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int get_from_map(const int[:, ::1] map, int row, int col):
    if row < 0 or col < 0 or row >= map.shape[0] or col >= map.shape[1]:
        return 0
    return map[row, col]

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int get_cost(const int[:, ::1] map, int nx, int ny):
    cdef int row, col, nint, nout, max_cost
    cdef int halo = 2
    max_cost = 0
    for row in range(map.shape[0]):
        for col in range(map.shape[1]):
            nint = map[row, col]
            nout = 0
            if get_from_map(map, row - 1, col - 1) > 0: nout += halo * halo
            if get_from_map(map, row + 1, col - 1) > 0: nout += halo * halo
            if get_from_map(map, row - 1, col + 1) > 0: nout += halo * halo
            if get_from_map(map, row + 1, col + 1) > 0: nout += halo * halo
            if get_from_map(map, row - 1, col) > 0: nout += halo * ny
            if get_from_map(map, row + 1, col) > 0: nout += halo * ny
            if get_from_map(map, row, col - 1) > 0: nout += halo * nx
            if get_from_map(map, row, col + 1) > 0: nout += halo * nx
            max_cost = max(max_cost, nint + 10 * nout)  # for now assume halo cells are 10x as expensive as interior cells
    return max_cost

def find_subdiv_solutions(const int[:, ::1] mask not None, int nx, int ny, int ncpus):
    cdef int[:, ::1] map
    cost = -1
    solution = None
    for ioffset in range(1 - nx, 1):
        for joffset in range(1 - ny, 1):
            map = get_map(mask, nx, ny, ioffset, joffset, ncpus)
            if map is not None:
                current_cost = get_cost(map, nx, ny)
                if cost == -1 or current_cost < cost:
                    cost = current_cost
                    solution = (ioffset, joffset, cost, numpy.asarray(map))
    return solution