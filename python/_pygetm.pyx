# cython: language_level=3

cimport numpy
import numpy

cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz)
cdef extern domain_initialize_open_boundaries(void* domain, int nwb, int nnb, int neb, int nsb, int nbdyp)
cdef extern void* domain_get_grid(void* domain, char grid_type)
cdef extern void domain_initialize(void* grid, int runtype)
cdef extern void domain_finalize(void* domain)
cdef extern void domain_update_depths(void* domain)
cdef extern double* grid_get_array(void* grid, char* name)
cdef extern void* advection_create()
cdef extern void advection_calculate(void* advection, int scheme, void* domain, double* pu, double* pv, double timestep, double* pvar)
cdef extern void* momentum_create(int runtype, void* pdomain, void* padvection, int advection_scheme, int apply_bottom_friction)
cdef extern double* momentum_get_array(void* momentum, char* name)
cdef extern void momentum_uv_momentum_2d(void* momentum, int runtype, double timestep, double* ptausx, double* ptausy, double* pdpdx, double* pdpdy)
cdef extern void* pressure_create(int runtype, void* pdomain)
cdef extern double* pressure_get_array(void* pressure, char* name)
cdef extern void pressure_surface(void* pressure, double* pz, double* psp)
cdef extern void* sealevel_create(void* pdomain)
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV)
cdef extern void sealevel_update_uvx(void* sealevel)

cdef class Grid:
    cdef void* p
    cdef int nx, ny

    def __init__(self, Domain domain, int grid_type):
        self.p = domain_get_grid(domain.p, grid_type)
        self.nx, self.ny = domain.nx, domain.ny
        if grid_type == ord('X'):
            self.nx += 1
            self.ny += 1

    def get_array(self, bytes name, dtype):
        cdef void* p = grid_get_array(self.p, name)
        if dtype == 0:
            return numpy.asarray(<double[:self.ny, :self.nx:1]> p)
        else:
            return numpy.asarray(<int[:self.ny, :self.nx:1]> p)

cdef class Domain:
    cdef void* p
    cdef public int halox, haloy, haloz
    cdef int nx, ny

    def __init__(self, int imin, int imax, int jmin, int jmax, int kmin, int kmax):
        self.p = domain_create(imin, imax, jmin, jmax, kmin, kmax, &self.halox, &self.haloy, &self.haloz)
        self.nx = imax - imin + 1 + 2 * self.halox
        self.ny = jmax - jmin + 1 + 2 * self.haloy

    def __dealloc__(self):
        domain_finalize(self.p)

    def get_grid(self, int grid_type):
        return Grid(self, grid_type)

    def update_depths(self):
        domain_update_depths(self.p)

    def initialize(self, int runtype):
        domain_initialize(self.p, runtype)

cdef class Advection:
    cdef void* p
    cdef void* pdomain
    cdef int scheme

    def __init__(self, Domain domain, int scheme):
        self.p = advection_create()
        self.scheme = scheme
        self.pdomain = domain.p

    def calculate(self, double [:, ::1] pu not None, double [:, ::1] pv not None, double timestep, double [:, ::1] pvar):
        advection_calculate(self.p, self.scheme, self.pdomain, &pu[0,0], &pv[0,0], timestep, &pvar[0,0])

cdef class Simulation:
    cdef public int runtype
    cdef void* pmomentum
    cdef void* ppressure
    cdef void* psealevel
    cdef int nx, ny

    def __init__(self, Domain domain, int runtype, int advection_scheme, int apply_bottom_friction):
        domain.initialize(runtype)
        padvection = advection_create()
        self.runtype = runtype
        self.pmomentum = momentum_create(runtype, domain.p, padvection, advection_scheme, apply_bottom_friction)
        self.ppressure = pressure_create(runtype, domain.p)
        self.psealevel = sealevel_create(domain.p)
        self.nx, self.ny = domain.nx, domain.ny

    def uv_momentum_2d(self, double timestep, double [:, ::1] tausx not None, double [:, ::1] tausy not None, double [:, ::1] dpdx not None, double [:, ::1] dpdy not None):
        momentum_uv_momentum_2d(self.pmomentum, self.runtype, timestep, &tausx[0,0], &tausy[0,0], &dpdx[0,0], &dpdy[0,0])

    def update_surface_pressure_gradient(self, double [:, ::1] pz not None, double [:, ::1] psp not None):
        pressure_surface(self.ppressure, &pz[0,0], &psp[0,0])

    def update_sealevel(self, double timestep, double [:, ::1] pU not None, double [:, ::1] pV not None):
        sealevel_update(self.psealevel, timestep, &pU[0,0], &pV[0,0])

    def update_sealevel_uvx(self):
        sealevel_update_uvx(self.psealevel)

    def get_array(self, bytes name, int source, int dtype=0):
        cdef double* p
        if source == 0:
            p = momentum_get_array(self.pmomentum, name)
        else:
            p = pressure_get_array(self.ppressure, name)
        return numpy.asarray(<double[:self.ny, :self.nx:1]> p)