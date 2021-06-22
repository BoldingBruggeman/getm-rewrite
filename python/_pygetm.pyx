# cython: language_level=3

cimport numpy
import numpy

cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz) nogil
cdef extern domain_initialize_open_boundaries(void* domain, int nwb, int nnb, int neb, int nsb, int nbdyp) nogil
cdef extern void* domain_get_grid(void* domain, int grid_type) nogil
cdef extern void domain_initialize(void* grid, int runtype, double* maxdt) nogil
cdef extern void domain_finalize(void* domain) nogil
cdef extern void domain_update_depths(void* domain) nogil
cdef extern void grid_interp_x(void* grid, double* source, double* target, int ioffset) nogil
cdef extern void grid_interp_y(void* grid, double* source, double* target, int joffset) nogil
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* data_type, void** p) nogil
cdef extern void* advection_create(int scheme, void* tgrid, void** p) nogil
cdef extern void advection_2d_calculate(int direction, void* advection, void* tgrid, void* ugrid, double* pu, double timestep, double* pvar) nogil
cdef extern void* momentum_create(int runtype, void* pdomain, int apply_bottom_friction) nogil
cdef extern void momentum_uv_momentum_2d(void* momentum, int runtype, double timestep, double* ptausx, double* ptausy, double* pdpdx, double* pdpdy) nogil
cdef extern void* pressure_create(int runtype, void* pdomain) nogil
cdef extern void pressure_surface(void* pressure, double* pz, double* psp) nogil
cdef extern void* sealevel_create(void* pdomain) nogil
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV) nogil
cdef extern void sealevel_update_uvx(void* sealevel) nogil

cpdef enum:
    TGRID = 1
    UGRID = 2
    VGRID = 3
    XGRID = 4
    UUGRID = -1
    VVGRID = -2
    UVGRID = -3
    VUGRID = -4

cdef class Array:
    cdef void* p
    cdef readonly numpy.ndarray all_values
    cdef readonly Grid grid

    cdef wrap_c_array(self, Domain domain, int source, void* obj, bytes name):
        cdef int data_type
        cdef int grid_type
        get_array(source, obj, name, &grid_type, &data_type, &self.p)
        self.grid = domain.grids[grid_type]
        if data_type == 0:
            self.all_values = numpy.asarray(<double[:self.grid.ny, :self.grid.nx:1]> self.p)
        else:
            self.all_values = numpy.asarray(<int[:self.grid.ny, :self.grid.nx:1]> self.p)
        self.finish_initialization()

    def wrap_ndarray(self, Grid grid, numpy.ndarray data):
        self.grid = grid
        self.all_values = numpy.ascontiguousarray(data)
        self.p = self.all_values.data
        self.finish_initialization()

cdef class Grid:
    cdef void* p
    cdef readonly int nx, ny
    cdef readonly Domain domain

    def __init__(self, Domain domain, int grid_type):
        self.domain = domain
        self.p = domain_get_grid(domain.p, grid_type)
        self.nx, self.ny = domain.nx, domain.ny
        if grid_type == XGRID:
            self.nx += 1
            self.ny += 1
        domain.grids[grid_type] = self

    def wrap(self, Array ar, bytes name):
        ar.wrap_c_array(self.domain, 0, self.p, name)
        return ar

    def interp_x(self, Array source, Array target, int offset):
        grid_interp_x(self.p, <double *>source.p, <double *>target.p, offset)

    def interp_y(self, Array source, Array target, int offset):
        grid_interp_y(self.p, <double *>source.p, <double *>target.p, offset)

cdef class Domain:
    cdef void* p
    cdef readonly int halox, haloy, haloz
    cdef readonly double maxdt
    cdef int nx, ny

    def __init__(self, int imin, int imax, int jmin, int jmax, int kmin, int kmax):
        self.p = domain_create(imin, imax, jmin, jmax, kmin, kmax, &self.halox, &self.haloy, &self.haloz)
        self.nx = imax - imin + 1 + 2 * self.halox
        self.ny = jmax - jmin + 1 + 2 * self.haloy
        self.grids = {}

    def __dealloc__(self):
        if self.p != NULL:
            domain_finalize(self.p)

    def update_depths(self):
        domain_update_depths(self.p)

    def initialize(self, int runtype):
        domain_initialize(self.p, runtype, &self.maxdt)

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
        self.D = numpy.asarray(<double[:self.tgrid.ny, :self.tgrid.nx:1]> pD)

    def calculate(self, Array u not None, Array v not None, double timestep, Array var not None):
        self.D[...] = self.tgrid.D.all_values
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, 0.5 * timestep, <double *>var.p)
        var.update_halos(2)
        advection_2d_calculate(1, self.p, self.tgrid.p, self.ugrid.p, <double *>u.p, timestep, <double *>var.p)
        var.update_halos(1)
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, 0.5 * timestep, <double *>var.p)

cdef class Simulation:
    cdef readonly Domain domain
    cdef readonly int runtype
    cdef void* pmomentum
    cdef void* ppressure
    cdef void* psealevel
    cdef int nx, ny

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

    def update_sealevel_uvx(self):
        sealevel_update_uvx(self.psealevel)

    def wrap(self, Array ar not None, bytes name, int source):
        cdef void* obj = self.pmomentum
        if (source == 2): obj = self.ppressure
        ar.wrap_c_array(self.domain, source, obj, name)
        return ar
