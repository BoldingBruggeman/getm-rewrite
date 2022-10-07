# cython: language_level=3
# cython: profile=True

cimport cython

from libc.math cimport ceil

cimport numpy
import numpy

from pygetm.constants import *

cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz) nogil
cdef extern void domain_initialize_open_boundaries(void* domain, int nbdyp, int nwb, int nnb, int neb, int nsb, int* bdy_info, int* bdy_i, int* bdy_j) nogil
cdef extern void* domain_get_grid(void* domain, int grid_type, int imin, int imax, int jmin, int jmax, int kmin, int kmax, int halox, int haloy, int haloz) nogil
cdef extern void domain_initialize(void* domain, int runtype, double Dmin, int method_vertical_coordinates, double ddl, double ddu, double Dgamma, int gamma_surf, double* maxdt) nogil
cdef extern void domain_finalize(void* domain) nogil
cdef extern void domain_do_vertical(void* domain, double timestep) nogil
cdef extern void grid_interp_x(int nx, int ny, int nz, double* source, double* target, int ioffset) nogil
cdef extern void grid_interp_y(int nx, int ny, int nz, double* source, double* target, int joffset) nogil
cdef extern void grid_interp_z(int nx, int ny, int nz1, int nz2, double* source, double* target, int koffset) nogil
cdef extern void grid_interp_xy(int nx1, int ny1, int nx2, int ny2, int nz, double* source, double* target, int ioffset, int joffset) nogil
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* sub_type, int* data_type, void** p) nogil
cdef extern void* advection_create(int scheme, void* tgrid) nogil
cdef extern void advection_uv_calculate(int direction, int nk, void* advection, void* tgrid, void* ugrid, double* pu, double* Ah, double timestep, double* pD, double* pDU, double* pvar) nogil
cdef extern void advection_w_calculate(void* padvection, void* tgrid, double* pw, double* pw_var, double timestep, double* ph, double* pvar)
cdef extern void* vertical_diffusion_create(void* tgrid) nogil
cdef extern void vertical_diffusion_calculate(void* diffusion, void* tgrid, double molecular, double* pnuh, double timestep, double cnpar, double* pho, double* phn, double* pvar, double* pea2, double* pea4) nogil
cdef extern void* momentum_create(int runtype, void* pdomain, double Am0, double cnpar, int coriolis_scheme) nogil
cdef extern void momentum_u_2d(int direction, void* momentum, double timestep, double* ptausx, double* pdpdx) nogil
cdef extern void momentum_u_3d(int direction, void* momentum, double timestep, double* ptausx, double* pdpdx, double* pidpdx, double* pviscosity) nogil
cdef extern void momentum_w_3d(void* momentum, double timestep) nogil
cdef extern void momentum_uv_coriolis(int direction, void* momentum) nogil
cdef extern void momentum_uv_coriolis_3d(int direction, void* momentum) nogil
cdef extern void momentum_bottom_friction_2d(void* momentum, int runtype) nogil
cdef extern void momentum_bottom_friction_3d(void* momentum) nogil
cdef extern void momentum_shear_frequency(void* momentum, double* pviscosity) nogil
cdef extern void momentum_stresses(void* momentum, double* tausx, double* tausy) nogil
cdef extern void momentum_diffusion_driver(void* momentum, int nk, double* h, double* hx, double* u, double* v, double* diffu, double* )
cdef extern void* pressure_create(int runtype, void* pdomain, int method_internal_pressure) nogil
cdef extern void pressure_surface(void* pressure, double* pz, double* psp) nogil
cdef extern void pressure_internal(void* pressure, double* buoy) nogil
cdef extern void* sealevel_create(void* pdomain) nogil
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV, double* pfwf) nogil
cdef extern void sealevel_boundaries(void* sealevel, void* momentum, double timestep) nogil
cdef extern void c_exponential_profile_1band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* initial, int up, double* out) nogil
cdef extern void c_exponential_profile_1band_centers(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* top, double* out) nogil
cdef extern void c_exponential_profile_2band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* f, double* k1, double* k2, double* initial, int up, double* out) nogil
cdef extern void c_thickness2center_depth(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* out) nogil
cdef extern void c_thickness2vertical_coordinates(int nx, int ny, int nz, int* mask, double* bottom_depth, double* h, double* zc, double* zf) nogil
cdef extern void c_alpha(int n, double* D, double Dmin, double Dcrit, int* mask, double* alpha)
cdef extern void c_clip_z(int n, double* z, double* H, double Dmin, int* mask)
cdef extern void c_horizontal_diffusion(int imin,int imax,int jmin,int jmax,int halox,int haloy,int* umask,int* vmask,double* idxu,double* dyu,double* idyv,double* dxv,double* Ah_u,double* Ah_v,int* tmask,double* iA,double dt,double* f,double* out)

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
    cdef numpy.ndarray _array
    cdef readonly Grid grid
    cdef readonly int on_boundary

    def __init__(self, Grid grid=None):
        if grid is not None:
            self.grid = grid

    @property
    def all_values(self):
        return self._array

    @all_values.setter
    def all_values(self, value):
        assert value is self._array

    cdef wrap_c_array(self, Domain domain, int source, void* obj, bytes name, register=True):
        cdef int grid_type
        cdef int sub_type
        cdef int data_type
        get_array(source, obj, name, &grid_type, &sub_type, &data_type, &self.p)
        if self.p == NULL:
            return
        self.on_boundary = False
        self.grid = domain.grids[grid_type]
        if sub_type == 0:
            # Horizontal-only array on normal grid
            if data_type == 0:
                self._array = numpy.asarray(<double[:self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:self.grid.ny_, :self.grid.nx_:1]> self.p)
        elif sub_type == 1:
            # Horizontal-only array on open boundaries
            if data_type == 0:
                self._array = numpy.asarray(<double[:self.grid.nbdyp:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:self.grid.nbdyp:1]> self.p)
            self.on_boundary = True
        elif sub_type == 2:
            # Depth-explicit array on normal grid, layer centers
            if data_type == 0:
                self._array = numpy.asarray(<double[:self.grid.nz_, :self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:self.grid.nz_, :self.grid.ny_, :self.grid.nx_:1]> self.p)
        elif sub_type == 3:
            # Depth-explicit array on normal grid, layer interfaces
            if data_type == 0:
                self._array = numpy.asarray(<double[:self.grid.nz_ + 1, :self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:self.grid.nz_ + 1, :self.grid.ny_, :self.grid.nx_:1]> self.p)
        else:
            # Depth-explicit array on open boundaries, layer centers
            assert sub_type == 4, 'Subtypes other than 0,1,2,3,4 not yet implemented'
            if data_type == 0:
                self._array = numpy.asarray(<double[:self.grid.nbdyp, :self.grid.nz_:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:self.grid.nbdyp, :self.grid.nz_:1]> self.p)
            self.on_boundary = True
        if self._fill_value is None:
            self._fill_value = self._array.flat[0]
        else:
            self._array[...] = self._fill_value
        self.finish_initialization()
        if register:
            self.register()
        return self

    def wrap_ndarray(self, numpy.ndarray data not None, on_boundary=False, register=True):
        self.on_boundary = on_boundary
        if self.on_boundary:
            assert data.ndim in (1, 2) and data.flags['C_CONTIGUOUS'], 'Invalid array properties for wrapping: %i dimensions, flags %s' % (data.ndim, data.flags)
            assert data.shape[0] == self.grid.domain.open_boundaries.np, 'Incorrect shape of first dimension (number of boundary points): expected %i, got %i' % (self.grid.domain.open_boundaries.np, data.shape[0])
            assert data.ndim == 1 or data.shape[1] == self.grid.nz_ or data.shape[1] == self.grid.nz_ + 1, 'Incorrect shape of second dimension (number of layers): expected %i or %i, got %i' % (self.grid.nz_, self.grid.nz_ + 1, data.shape[1])
        elif data.ndim != 0:
            assert data.ndim in (2, 3) and data.flags['C_CONTIGUOUS'], 'Invalid array properties for wrapping: %i dimensions, flags %s' % (data.ndim, data.flags)
            assert data.shape[data.ndim - 1] == self.grid.nx_ and data.shape[data.ndim - 2] == self.grid.ny_, 'Incorrect horizontal extent: expected (ny=%i,nx=%i), got (ny=%i,nx=%i)' % (self.grid.ny_, self.grid.nx_, data.shape[data.ndim - 2], data.shape[data.ndim - 1])
            assert data.ndim == 2 or data.shape[0] == self.grid.nz_ or data.shape[0] == self.grid.nz_ + 1, 'Incorrect vertical extent: expected %i or %i, got %i' % (self.grid.nz_, self.grid.nz_ + 1, data.shape[0])
        self._array = data
        self.p = self._array.data
        self.finish_initialization()
        if register:
            self.register()

cdef class Grid:
    cdef void* p
    cdef readonly int nx, ny, nz
    cdef readonly int nx_, ny_, nz_
    cdef readonly Domain domain

    def __init__(self, Domain domain, int grid_type):
        self.domain = domain
        self.nx, self.ny, self.nz = domain.nx, domain.ny, domain.nz
        if grid_type == XGRID:
            self.nx += 1
            self.ny += 1
        self.p = domain_get_grid(domain.p, grid_type, 1, self.nx, 1, self.ny, 1, self.nz, domain.halox, domain.haloy, domain.haloz)
        self.nx_, self.ny_, self.nz_ = self.nx + 2 * domain.halox, self.ny + 2 * domain.haloy, self.nz + 2 * domain.haloz
        domain.grids[grid_type] = self

    def wrap(self, Array ar, bytes name, register=True):
        return ar.wrap_c_array(self.domain, 0, self.p, name, register)

def interp_x(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_x(source.shape[2], source.shape[1], source.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_y(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_y(source.shape[2], source.shape[1], source.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_z(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_z(source.shape[2], source.shape[1], source.shape[0], target.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_xy(double[:,:,::1] source not None, double[:,:,::1] target not None, int ioffset, int joffset):
    grid_interp_xy(source.shape[2], source.shape[1], target.shape[2], target.shape[1], source.shape[0], &source[0,0,0], &target[0,0,0], ioffset, joffset)

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

    def do_vertical(self, double timestep):
        domain_do_vertical(self.p, timestep)

    def initialize(self, int runtype, double Dmin, int method_vertical_coordinates, double ddl=0., double ddu=0., double Dgamma=0., gamma_surf=True):
        domain_initialize(self.p, runtype, Dmin, method_vertical_coordinates, ddl, ddu, Dgamma, 1 if gamma_surf else 0, &self.maxdt)

    def initialize_open_boundaries(self, int nwb, int nnb, int neb, int nsb, int nbdyp, int[:,::1] bdy_info, int[::1] bdy_i, int[::1] bdy_j):
        assert bdy_info.shape[0] == 6, 'bdy_info should have 6 rows'
        assert bdy_info.shape[1] == nwb + nnb + neb + nsb, 'bdy_info should have as many columns as the number of boundaries'
        assert bdy_i.shape[0] == nbdyp, 'bdy_i should have a length equal to the number of boundary points'
        assert bdy_j.shape[0] == nbdyp, 'bdy_j should have a length equal to the number of boundary points'
        domain_initialize_open_boundaries(self.p, nbdyp, nwb, nnb, neb, nsb, &bdy_info[0,0], &bdy_i[0], &bdy_j[0])

cdef class Advection:
    cdef void* p
    cdef readonly Grid grid
    cdef readonly Grid ugrid
    cdef readonly Grid vgrid
    cdef readonly numpy.ndarray D, h
    cdef double* ph
    cdef double* phu
    cdef double* phv
    cdef double* pDU
    cdef double* pDV

    def __init__(self, Grid grid, int scheme):
        self.grid = grid
        self.ugrid = grid.ugrid
        self.vgrid = grid.vgrid
        self.p = advection_create(scheme, self.grid.p)

        # Work array for layer heights that change during 3d advection
        # This is shared by the temporary column height that evolves similarly during 2d advection
        self.h = numpy.empty((self.grid.nz_, self.grid.ny_, self.grid.nx_))
        self.D = self.h[0,...]

        # Store references to column height and layer thicknesses
        self.pDU = <double*>(<Array>self.ugrid.D).p
        self.pDV = <double*>(<Array>self.vgrid.D).p
        self.phu = <double*>(<Array>self.ugrid.hn).p
        self.phv = <double*>(<Array>self.vgrid.hn).p
        self.ph = <double*>self.h.data

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def u_2d(self, Array u not None, double timestep, Array var not None, Array Ah=None):
        cdef double * pAh = <double *>Ah.p if Ah is not None else NULL
        advection_uv_calculate(1, 1, self.p, self.grid.p, self.ugrid.p, <double *>u.p, pAh, timestep, self.ph, self.pDU, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def v_2d(self, Array v not None, double timestep, Array var not None, Array Ah=None):
        cdef double * pAh = <double *>Ah.p if Ah is not None else NULL
        advection_uv_calculate(2, 1, self.p, self.grid.p, self.vgrid.p, <double *>v.p, pAh, timestep, self.ph, self.pDV, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def u_3d(self, Array u not None, double timestep, Array var not None, Array Ah=None):
        cdef double * pAh = <double *>Ah.p if Ah is not None else NULL
        advection_uv_calculate(1, var._array.shape[0], self.p, self.grid.p, self.ugrid.p, <double *>u.p, pAh, timestep, self.ph, self.phu, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def v_3d(self, Array v not None, double timestep, Array var not None, Array Ah=None):
        cdef double * pAh = <double *>Ah.p if Ah is not None else NULL
        advection_uv_calculate(2, var._array.shape[0], self.p, self.grid.p, self.vgrid.p, <double *>v.p, pAh, timestep, self.ph, self.phv, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def w_3d(self, Array w not None, Array w_var not None, double timestep, Array var not None):
        advection_w_calculate(self.p, self.grid.p, <double *>w.p, <double *>w_var.p, timestep, self.ph, <double *>var.p)

cdef class VerticalDiffusion:
    cdef void* p
    cdef Grid tgrid
    cdef double cnpar
    cdef Array ho, hn

    def __init__(self, Grid grid, double cnpar=1.):
        self.tgrid = grid
        self.p = vertical_diffusion_create(self.tgrid.p)
        self.cnpar = cnpar
        self.ho = self.tgrid.ho
        self.hn = self.tgrid.hn

    def __call__(self, Array nuh not None, double timestep, Array var not None, double molecular=0., Array ea2=None, Array ea4=None, bint use_ho=False):
        cdef double* pea2 = NULL
        cdef double* pea4 = NULL
        cdef double* pho = <double *>(self.ho if use_ho else self.hn).p
        if ea2 is not None:
            pea2 = <double *>ea2.p
        if ea4 is not None:
            pea4 = <double *>ea4.p
        vertical_diffusion_calculate(self.p, self.tgrid.p, molecular, <double *>nuh.p, timestep, self.cnpar, pho, <double *>self.hn.p, <double *>var.p, pea2, pea4)

cdef class Simulation:
    cdef readonly Domain domain
    cdef readonly int runtype
    cdef void* ppressure
    cdef void* psealevel
    cdef readonly int nx, ny

    def __init__(self, Domain domain, int runtype, int internal_pressure_method=1):
        self.domain = domain
        domain.initialize(runtype)
        self.runtype = runtype
        self.ppressure = pressure_create(runtype, domain.p, internal_pressure_method)
        self.psealevel = sealevel_create(domain.p)
        self.nx, self.ny = domain.nx, domain.ny

    def update_surface_pressure_gradient(self, Array z not None, Array sp not None):
        assert z.grid is self.domain.T
        assert sp.grid is self.domain.T
        pressure_surface(self.ppressure, <double *>z.p, <double *>sp.p)

    def update_internal_pressure_gradient(self, Array buoy not None):
        assert buoy.grid is self.domain.T and buoy.z == CENTERS
        pressure_internal(self.ppressure, <double *>buoy.p)

    def advance_surface_elevation(self, double timestep, Array U not None, Array V not None, Array fwf not None):
        assert U.grid is self.domain.U
        assert V.grid is self.domain.V
        sealevel_update(self.psealevel, timestep, <double *>U.p, <double *>V.p, <double *>fwf.p)

    def update_surface_elevation_boundaries(self, double timestep, Momentum momentum):
        sealevel_boundaries(self.psealevel, momentum.p, timestep)

    def wrap(self, Array ar not None, bytes name, int source):
        assert source in (2, 3)
        cdef void* obj = NULL
        if (source == 2):
            obj = self.ppressure
        elif (source == 3):
            obj = self.psealevel
        return ar.wrap_c_array(self.domain, source, obj, name)

cdef class Momentum:
    cdef void* p
    cdef readonly Domain domain

    def __init__(self, Domain domain, int runtype, double Am0, double cnpar, int coriolis_scheme):
        self.domain = domain
        self.p = momentum_create(runtype, domain.p, Am0, cnpar, coriolis_scheme)

    def u_2d(self, double timestep, Array tausx not None, Array dpdx not None):
        assert tausx.grid is self.domain.U, 'grid mismatch for tausx: expected %s, got %s' % (self.domain.U.postfix, tausx.grid.postfix)
        assert dpdx.grid is self.domain.U, 'grid mismatch for dpdx: expected %s, got %s' % (self.domain.U.postfix, dpdx.grid.postfix)
        momentum_u_2d(1, self.p, timestep, <double *>tausx.p, <double *>dpdx.p)

    def v_2d(self, double timestep, Array tausy not None, Array dpdy not None):
        assert tausy.grid is self.domain.V, 'grid mismatch for tausy: expected %s, got %s' % (self.domain.V.postfix, tausy.grid.postfix)
        assert dpdy.grid is self.domain.V, 'grid mismatch for dpdy: expected %s, got %s' % (self.domain.V.postfix, dpdy.grid.postfix)
        momentum_u_2d(2, self.p, timestep, <double *>tausy.p, <double *>dpdy.p)

    def coriolis_fu(self):
        momentum_uv_coriolis(1, self.p)

    def coriolis_fv(self):
        momentum_uv_coriolis(2, self.p)

    def pk_3d(self, double timestep, Array tausx not None, Array dpdx not None, Array idpdx not None, Array viscosity_u not None):
        assert tausx.grid is self.domain.U, 'grid mismatch for tausx: expected %s, got %s' % (self.domain.U.postfix, tausx.grid.postfix)
        assert dpdx.grid is self.domain.U, 'grid mismatch for dpdx: expected %s, got %s' % (self.domain.U.postfix, dpdx.grid.postfix)
        assert idpdx.grid is self.domain.U, 'grid mismatch for idpdx: expected %s, got %s' % (self.domain.U.postfix, idpdx.grid.postfix)
        assert viscosity_u.grid is self.domain.U and viscosity_u.z == INTERFACES, 'grid mismatch for viscosity_u: expected %s, got %s' % (self.domain.U.postfix, viscosity_u.grid.postfix)
        momentum_u_3d(1, self.p, timestep, <double *>tausx.p, <double *>dpdx.p, <double *>idpdx.p, <double *>viscosity_u.p)

    def qk_3d(self, double timestep, Array tausy not None, Array dpdy not None, Array idpdy not None, Array viscosity_v not None):
        assert tausy.grid is self.domain.V, 'grid mismatch for tausy: expected %s, got %s' % (self.domain.V.postfix, tausy.grid.postfix)
        assert dpdy.grid is self.domain.V, 'grid mismatch for dpdy: expected %s, got %s' % (self.domain.V.postfix, dpdy.grid.postfix)
        assert idpdy.grid is self.domain.V, 'grid mismatch for idpdy: expected %s, got %s' % (self.domain.V.postfix, idpdy.grid.postfix)
        assert viscosity_v.grid is self.domain.V and viscosity_v.z == INTERFACES, 'grid mismatch for viscosity_v: expected %s, got %s' % (self.domain.V.postfix, viscosity_v.grid.postfix)
        momentum_u_3d(2, self.p, timestep, <double *>tausy.p, <double *>dpdy.p, <double *>idpdy.p, <double *>viscosity_v.p)

    def coriolis_fpk(self):
        momentum_uv_coriolis_3d(1, self.p)

    def coriolis_fqk(self):
        momentum_uv_coriolis_3d(2, self.p)

    def bottom_friction_2d(self, int update_z0b):
        momentum_bottom_friction_2d(self.p, 1 if update_z0b else 4)

    def bottom_friction_3d(self):
        momentum_bottom_friction_3d(self.p)

    def w_3d(self, double timestep):
        momentum_w_3d(self.p, timestep)

    def momentum_diffusion_driver(self, Array h not None, Array hx not None, Array u not None, Array v not None, Array diffu not None, Array diffv not None):
        assert h.grid is self.domain.T
        assert hx.grid is self.domain.X
        assert u.grid is self.domain.U
        assert v.grid is self.domain.V
        assert diffu.grid is self.domain.U
        assert diffv.grid is self.domain.V
        cdef int nk = h._array.shape[0] if h.z else 1
        momentum_diffusion_driver(self.p, nk, <double*> h.p, <double*> hx.p, <double*> u.p, <double*> v.p, <double*> diffu.p, <double*> diffv.p)

    def update_stresses(self, Array tausx not None, Array tausy not None):
        assert tausx.grid is self.domain.T, 'grid mismatch for tausx: expected %s, got %s' % (self.domain.T.postfix, tausx.grid.postfix)
        assert tausy.grid is self.domain.T, 'grid mismatch for tausy: expected %s, got %s' % (self.domain.T.postfix, tausy.grid.postfix)
        momentum_stresses(self.p, <double *>tausx.p, <double *>tausy.p)

    def update_shear_frequency(self, Array viscosity not None):
        assert viscosity.grid is self.domain.T and viscosity.z == INTERFACES, 'grid mismatch for viscosity: expected %s, got %s' % (self.domain.T.postfix, viscosity.grid.postfix)
        momentum_shear_frequency(self.p, <double *>viscosity.p)

    def wrap(self, Array ar not None, bytes name):
        return ar.wrap_c_array(self.domain, 1, self.p, name)

def exponential_profile_1band_interfaces(Array mask not None, Array h not None, Array k not None, Array initial not None, bint up=False, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is k.grid and k.z == CENTERS
    assert mask.grid is initial.grid and not initial.z
    assert mask.grid is out.grid and out.z == INTERFACES
    c_exponential_profile_1band_interfaces(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>k.p, <double *>initial.p, up, <double *>out.p)

def exponential_profile_1band_centers(Array mask not None, Array h not None, Array k not None, Array top not None, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is k.grid and k.z == CENTERS
    assert mask.grid is top.grid and not top.z
    assert mask.grid is out.grid and out.z == CENTERS
    c_exponential_profile_1band_centers(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>k.p, <double *>top.p, <double *>out.p)

def exponential_profile_2band_interfaces(Array mask not None, Array h not None, Array f not None, Array k1 not None, Array k2 not None, Array initial not None, bint up=False, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is f.grid and not f.z
    assert mask.grid is k1.grid and not k1.z
    assert mask.grid is k2.grid and not k2.z
    assert mask.grid is initial.grid and not initial.z
    assert mask.grid is out.grid and out.z == INTERFACES
    c_exponential_profile_2band_interfaces(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>f.p, <double *>k1.p, <double *>k2.p, <double *>initial.p, up, <double *>out.p)

def thickness2center_depth(Array mask not None, Array h not None, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    if out is None:
        out = h.grid.array(z=CENTERS)
    assert mask.grid is out.grid and out.z == CENTERS
    c_thickness2center_depth(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>out.p)
    return out

def thickness2vertical_coordinates(Array mask not None, Array H not None, Array h not None, Array zc not None, Array zf not None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is H.grid and not H.z
    assert mask.grid is zc.grid and zc.z == CENTERS
    assert mask.grid is zf.grid and zf.z == INTERFACES
    c_thickness2vertical_coordinates(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, <int *>mask.p, <double *>H.p, <double *>h.p, <double *>zc.p, <double *>zf.p)

def alpha(Array D not None, double Dmin, double Dcrit, Array out not None):
    cdef Array mask
    mask = D.grid.mask
    c_alpha(D._array.size, <double*>D.p, Dmin, Dcrit, <int*>mask.p, <double*>out.p)

def clip_z(Array z not None, double Dmin):
    cdef Array mask, H
    mask = z.grid.mask
    H = z.grid.H
    c_clip_z(z._array.size, <double*>z.p, <double*>H.p, Dmin, <int*>mask.p)

def horizontal_diffusion(Array f not None, Array Ah_u not None, Array Ah_v not None, double dt, Array out not None):
    cdef int halox = f.grid.domain.halox
    cdef int haloy = f.grid.domain.haloy
    cdef int nx = f.grid.nx
    cdef int ny = f.grid.ny
    cdef Grid ugrid = f.grid.ugrid
    cdef Grid vgrid = f.grid.vgrid
    cdef Array umask = ugrid.mask
    cdef Array vmask = vgrid.mask
    cdef Array mask = f.grid.mask
    cdef Array iarea = f.grid.iarea
    cdef Array idx = ugrid.idx
    cdef Array dy = ugrid.dy
    cdef Array idy = vgrid.idy
    cdef Array dx = vgrid.dx
    assert not f.z
    assert Ah_u.grid is ugrid and not Ah_u.z
    assert Ah_v.grid is vgrid and not Ah_v.z
    assert out.grid is f.grid and not out.z
    c_horizontal_diffusion(1, nx, 1, ny, halox, haloy,
        <int*>umask.p, <int*>vmask.p, <double*>idx.p, <double*>dy.p, <double*>idy.p, <double*>dx.p,
        <double*>Ah_u.p, <double*>Ah_v.p, <int*>mask.p, <double*>iarea.p, dt, <double*>f.p, <double*>out.p)

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.profile(False)
cdef int subdomain_used(const int[:, ::1] mask, int istart, int istop, int jstart, int jstop) nogil:
    for j in range(jstart, jstop):
        for i in range(istart, istop):
            if mask[j, i] != 0:
                return True
    return False

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.profile(False)
cdef int subdomain_count_cells(const int[:, ::1] mask, int istart, int istop, int jstart, int jstop) nogil:
    cdef int count = 0
    for j in range(jstart, jstop):
        for i in range(istart, istop):
            if mask[j, i] != 0:
                count += 1
    return count

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True)
cdef int decomposition_is_valid(const int[:, ::1] mask, int nx, int ny, int xoffset, int yoffset, int ncpus, double max_protrude) nogil:
    cdef int nrow = <int>ceil((mask.shape[0] - yoffset) / float(ny))
    cdef int ncol = <int>ceil((mask.shape[1] - xoffset) / float(nx))
    cdef int free_cpus = ncpus
    cdef int row, col, i, j, n
    if xoffset <= -nx * max_protrude or xoffset + ncol * nx - mask.shape[1] >= nx * max_protrude: return False
    if yoffset <= -ny * max_protrude or yoffset + nrow * ny - mask.shape[0] >= ny * max_protrude: return False
    if (nrow == 1 and yoffset != 0) or (ncol == 1 and xoffset != 0): return False
    if (nrow == 2 and yoffset != 0 and yoffset + nrow * ny > mask.shape[0]): return False
    if (ncol == 2 and xoffset != 0 and xoffset + ncol * nx > mask.shape[1]): return False
    if nrow * ncol < ncpus:
        # impossible to use that many cores with the current number of rows and columns
        return False
    for row in range(nrow):
        for col in range(ncol):
            if subdomain_used(mask, max(0, xoffset + col * nx), min(mask.shape[1], xoffset + (col + 1) * nx), max(0, yoffset + row * ny), min(mask.shape[0], yoffset + (row + 1) * ny)):
                if free_cpus == 0: return False
                free_cpus -= 1
    return free_cpus == 0

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.cdivision(True)
@cython.initializedcheck(False)
cdef void get_map(const int[:, ::1] mask, int nx, int ny, int xoffset, int yoffset, int[::1] map, int* nrow, int* ncol) nogil:
    cdef int row, col, i
    nrow[0] = <int>ceil((mask.shape[0] - yoffset) / float(ny))
    ncol[0] = <int>ceil((mask.shape[1] - xoffset) / float(nx))
    i = 0
    for row in range(nrow[0]):
        for col in range(ncol[0]):
            map[i] = subdomain_count_cells(mask, max(0, xoffset + col * nx), min(mask.shape[1], xoffset + (col + 1) * nx), max(0, yoffset + row * ny), min(mask.shape[0], yoffset + (row + 1) * ny))
            i += 1

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.profile(False)
@cython.initializedcheck(False)
cdef int get_from_map(const int[::1] map, int nrow, int ncol, int row, int col) nogil:
    if row < 0 or col < 0 or row >= nrow or col >= ncol:
        return 0
    return map[row * ncol + col]

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.initializedcheck(False)
cdef int get_cost(const int[::1] map, int nrow, int ncol, int nx, int ny, int weight_unmasked, int weight_any, int weight_halo) nogil:
    cdef int row, col, nint, nout, max_cost
    cdef int halo = 2
    max_cost = 0
    for row in range(nrow):
        for col in range(ncol):
            nint = map[row * ncol + col]  # unmasked cells only
            nout = 0
            if get_from_map(map, nrow, ncol, row - 1, col - 1) > 0: nout += halo * halo
            if get_from_map(map, nrow, ncol, row + 1, col - 1) > 0: nout += halo * halo
            if get_from_map(map, nrow, ncol, row - 1, col + 1) > 0: nout += halo * halo
            if get_from_map(map, nrow, ncol, row + 1, col + 1) > 0: nout += halo * halo
            if get_from_map(map, nrow, ncol, row - 1, col) > 0: nout += halo * nx
            if get_from_map(map, nrow, ncol, row + 1, col) > 0: nout += halo * nx
            if get_from_map(map, nrow, ncol, row, col - 1) > 0: nout += halo * ny
            if get_from_map(map, nrow, ncol, row, col + 1) > 0: nout += halo * ny
            max_cost = max(max_cost, weight_unmasked * nint + weight_halo * nout)
    if nx % 4 != 0: max_cost *= 2  # penalize non-alignment
    return max_cost + weight_any * nx * ny     # add an overhead for all cells - unmasked or not

def find_subdiv_solutions(const int[:, ::1] mask not None, int nx, int ny, int ncpus, int weight_unmasked, int weight_any, int weight_halo, double max_protrude):
    cdef int[::1] map
    cdef int nrow, ncol
    cost = -1
    map = numpy.empty(((mask.shape[1] // nx + 2) * (mask.shape[0] // ny + 2),), dtype=numpy.intc)
    for yoffset in range(1 - ny, 1):
        for xoffset in range(1 - nx, 1):
            if decomposition_is_valid(mask, nx, ny, xoffset, yoffset, ncpus, max_protrude):
                get_map(mask, nx, ny, xoffset, yoffset, map, &nrow, &ncol)
                current_cost = get_cost(map, nrow, ncol, nx, ny, weight_unmasked, weight_any, weight_halo)
                if cost == -1 or current_cost < cost:
                    cost = current_cost
                    best_xoffset = xoffset
                    best_yoffset = yoffset
    if cost != -1:
        get_map(mask, nx, ny, best_xoffset, best_yoffset, map, &nrow, &ncol)
        return (best_xoffset, best_yoffset, cost, numpy.asarray(map[:nrow * ncol]).reshape(nrow, ncol))
