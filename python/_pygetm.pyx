# cython: language_level=3

cimport numpy
import numpy

cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz)
cdef extern domain_initialize_open_boundaries(void* domain, int nwb, int nnb, int neb, int nsb, int nbdyp)
cdef extern void* domain_get_grid(void* domain, char grid_type)
cdef extern void domain_initialize(void* grid, int runtype, double* maxdt)
cdef extern void domain_finalize(void* domain)
cdef extern void domain_update_depths(void* domain)
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* data_type, void** p)
cdef extern void* advection_create()
cdef extern void advection_calculate(void* advection, int scheme, void* domain, double* pu, double* pv, double timestep, double* pvar)
cdef extern void* momentum_create(int runtype, void* pdomain, void* padvection, int advection_scheme, int apply_bottom_friction)
cdef extern void momentum_uv_momentum_2d(void* momentum, int runtype, double timestep, double* ptausx, double* ptausy, double* pdpdx, double* pdpdy)
cdef extern void* pressure_create(int runtype, void* pdomain)
cdef extern void pressure_surface(void* pressure, double* pz, double* psp)
cdef extern void* sealevel_create(void* pdomain)
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV)
cdef extern void sealevel_update_uvx(void* sealevel)

cpdef enum:
    TGRID = 1
    UGRID = 2
    VGRID = 3
    XGRID = 4

cdef class Array:
    cdef void* p
    cdef readonly numpy.ndarray all_values
    cdef readonly Grid grid

    cdef wrap(self, Domain domain, int source, void* obj, bytes name):
        cdef int data_type
        cdef int grid_type
        get_array(source, obj, name, &grid_type, &data_type, &self.p)
        self.grid = domain.grids[grid_type]
        if data_type == 0:
            self.all_values = numpy.asarray(<double[:self.grid.ny, :self.grid.nx:1]> self.p)
        else:
            self.all_values = numpy.asarray(<int[:self.grid.ny, :self.grid.nx:1]> self.p)
        self.finish_initialization()

    def empty(self, Grid grid, dtype):
        self.grid = grid
        self.all_values = numpy.empty((self.grid.ny, self.grid.nx), dtype=dtype)
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
        ar.wrap(self.domain, 0, self.p, name)
        return ar

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
        domain_finalize(self.p)

    def get_grid(self, int grid_type):
        return Grid(self, grid_type)

    def update_depths(self):
        domain_update_depths(self.p)

    def initialize(self, int runtype):
        domain_initialize(self.p, runtype, &self.maxdt)

cdef class Advection:
    cdef void* p
    cdef void* pdomain
    cdef int scheme

    def __init__(self, Domain domain, int scheme):
        self.p = advection_create()
        self.scheme = scheme
        self.pdomain = domain.p

    def calculate(self, Array u not None, Array v not None, double timestep, Array var not None):
        advection_calculate(self.p, self.scheme, self.pdomain, <double *>u.p, <double *>v.p, timestep, <double *>var.p)

cdef class Simulation:
    cdef readonly Domain domain
    cdef readonly int runtype
    cdef void* pmomentum
    cdef void* ppressure
    cdef void* psealevel
    cdef int nx, ny

    def __init__(self, Domain domain, int runtype, int advection_scheme, int apply_bottom_friction):
        self.domain = domain
        domain.initialize(runtype)
        padvection = advection_create()
        self.runtype = runtype
        self.pmomentum = momentum_create(runtype, domain.p, padvection, advection_scheme, apply_bottom_friction)
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

    def wrap(self, Array ar, bytes name, int source):
        cdef void* obj = self.pmomentum
        if (source == 2): obj = self.ppressure
        ar.wrap(self.domain, source, obj, name)
        return ar
