# cython: language_level=3
# cython: profile=True

cimport cython
include "version.pxi"

from libc.math cimport ceil

cimport numpy
import numpy

from pygetm.constants import *

cdef extern void c_allocate_array(int n, int dtype, void** ptype, void** pdata) nogil
cdef extern void c_deallocate_array(void* ptype) nogil
cdef extern void* domain_create(int imin, int imax, int jmin, int jmax, int kmin, int kmax, int* halox, int* haloy, int* haloz) nogil
cdef extern void* domain_get_grid(void* domain, int grid_type, int imin, int imax, int jmin, int jmax, int kmin, int kmax, int halox, int haloy, int haloz) nogil
cdef extern void domain_initialize(void* domain, int runtype, double Dmin, int method_vertical_coordinates, double ddl, double ddu, double Dgamma, int gamma_surf, double* maxdt) nogil
cdef extern void grid_finalize(void* grid) nogil
cdef extern void domain_finalize(void* domain) nogil
cdef extern void domain_do_vertical(void* domain, double timestep) nogil
cdef extern void grid_interp_x(int nx, int ny, int nz, double* source, double* target, int ioffset) nogil
cdef extern void grid_interp_y(int nx, int ny, int nz, double* source, double* target, int joffset) nogil
cdef extern void grid_interp_z(int nx, int ny, int nz1, int nz2, double* source, double* target, int koffset) nogil
cdef extern void grid_interp_xy(int nx1, int ny1, int nx2, int ny2, int nz, double* source, double* target, int ioffset, int joffset) nogil
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* sub_type, int* data_type, void** p) nogil
cdef extern void* advection_create(int scheme, void* tgrid) nogil
cdef extern void advection_finalize(void* advection) nogil
cdef extern void advection_uv_calculate(int direction, int nk, void* advection, void* tgrid, void* ugrid, double* pu, double* Ah, double timestep, double* pD, double* pDU, double* pvar) nogil
cdef extern void advection_w_calculate(void* padvection, void* tgrid, double* pw, double* pw_var, double timestep, double* ph, double* pvar)
cdef extern void* vertical_diffusion_create(void* tgrid) nogil
cdef extern void vertical_diffusion_finalize(void* diffusion) nogil
cdef extern void c_vertical_diffusion_prepare(void* diffusion, int nx, int ny, int nz, double molecular, double* pnuh, double timestep, double cnpar, int* pmask, double* pho, double* phn) nogil
cdef extern void c_vertical_diffusion_apply(void* diffusion, int nx, int ny, int nz, int* pmask, double* pho, double* phn, double* pvar, double* pea2, double* pea4) nogil
cdef extern void* momentum_create(int runtype, void* pdomain, double Am0, double cnpar, int coriolis_scheme) nogil
cdef extern void momentum_finalize(void* momentum) nogil
cdef extern void momentum_diffusion_driver(void* momentum, int nk, double* h, double* hx, double* u, double* v, double* diffu, double* diffv) nogil
cdef extern void c_exponential_profile_1band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* initial, int up, double* out) nogil
cdef extern void c_exponential_profile_1band_centers(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* top, double* out) nogil
cdef extern void c_exponential_profile_2band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* f, double* k1, double* k2, double* initial, int up, double* out) nogil
cdef extern void c_thickness2center_depth(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* out) nogil
cdef extern void c_thickness2vertical_coordinates(int nx, int ny, int nz, int* mask, double* bottom_depth, double* h, double* zc, double* zf) nogil
cdef extern void c_alpha(int n, double* D, double Dmin, double Dcrit, int* mask, double* alpha) nogil
cdef extern void c_clip_z(int n, double* z, double* H, double Dmin, int* mask) nogil
cdef extern void c_horizontal_diffusion(int imin,int imax,int jmin,int jmax,int halox,int haloy,int* umask,int* vmask,double* idxu,double* dyu,double* idyv,double* dxv,double* Ah_u,double* Ah_v,int* tmask,double* iA,double dt,double* f,double* out) nogil
cdef extern void c_bottom_friction(int nx, int ny, int* mask, double* u, double* v, double* D, double* z0b, double* z0b_in, double avmmol, double* ru, int iterate) nogil
cdef extern void c_collect_3d_momentum_sources(int nx, int ny, int nz, int halox, int haloy, int* mask, double* alpha, double* ho, double* hn, double* dp, double* cor, double* adv, double* diff, double* idp, double* taus, double* rr, double dt, double* ea2, double* ea4) nogil
cdef extern void c_advance_2d_transport(int ny, int ny, int halox, int haloy, int* mask, double* alpha, double* D, double* dp, double* taus, double* cor, double* adv, double* diff, double* damp, double* SA, double* SB, double* SD, double* SF, double* r, double dt, double* U) nogil
cdef extern void c_w_momentum_3d(int nx, int ny, int nz, int imin, int imax, int jmin, int jmax, const int* mask, const double* dyu, const double* dxv, const double* iarea, const double* ho, const double* hn, const double* pk, const double* qk, double dt, double* w) nogil
cdef extern void c_shear_frequency(int nx, int ny, int nz, int imin, int imax, int jmin, int jmax, const int* mask, const double* hu, const double* hv, const double* uk, const double* vk, double* SS) nogil
cdef extern void c_shear_frequency2(int nx, int ny, int nz, int imin, int imax, int jmin, int jmax, const int* mask, const double* h, const double* hu, const double* hv, const double* uk, const double* vk, const double* num, double* SS) nogil
cdef extern void c_surface_shear_velocity(int nx, int ny, int imin, int imax, int jmin, int jmax, const int* mask, const double* tausx, const double* tausy, double* ustar2_s) nogil
cdef extern void c_bottom_shear_velocity(int nx, int ny, int imin, int imax, int jmin, int jmax, const int* mask, const int* umask, const int* vmask, const double* uk_bot, const double* vk_bot, const double* rru, const double* rrv, double* taubx, double* tauby, double* ustar2_b) nogil
cdef extern void c_multiply_add(int n, double* tgt, double* add, double scale_factor) nogil
cdef extern void c_advance_surface_elevation(int nx, int ny, int halox, int haloy, int* mask, double* dyu, double* dxv, double* iarea, double* z, double* U, double* V, double* fwf, double dt) nogil
cdef extern void c_surface_pressure_gradient(int nx, int ny, int imin, int imax, int jmin, int jmax, int* umask, int* vmask, double* idxu, double* idyv, double* z, double* sp, double* H, double* D, double Dmin, double* dpdx, double* dpdy) nogil
cdef extern void c_blumberg_mellor(int nx, int ny, int nz, int imin, int imax, int jmin, int jmax, const int* umask, const int* vmask, const double* idxu, const double* idyv, const double* hu, const double* hv, const double* zf, const double* buoy, double* idpdx, double* idpdy) nogil
cdef extern void c_shchepetkin_mcwilliams(int nx, int ny, int nz, int imin, int imax, int jmin, int jmax, const int* mask, const int* umask, const int* vmask, const double* idxu, const double* idyv, const double* h, const double* z, const double* zc, const double* buoy, double* idpdx, double* idpdy) nogil


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
    cdef void* ptype
    cdef numpy.ndarray _array
    cdef readonly Grid grid
    cdef readonly bint on_boundary

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
        else:
            # Depth-explicit array on normal grid, layer centers or interface
            assert sub_type == 1 or sub_type == 2, 'Subtypes other than 0,1,2 not yet implemented'
            nz = self.grid.nz_ if sub_type == 1 else self.grid.nz_ + 1
            if data_type == 0:
                self._array = numpy.asarray(<double[:nz, :self.grid.ny_, :self.grid.nx_:1]> self.p)
            else:
                self._array = numpy.asarray(<int[:nz, :self.grid.ny_, :self.grid.nx_:1]> self.p)
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

    def allocate(self, shape, dtype):
        cdef int n = 0
        cdef bint is_float = dtype == float
        if len(shape) > 0:
            n = 1
            for l in shape:
                n *= l
        if (is_float or dtype == int) and n > 0:
            c_allocate_array(n, 0 if is_float else 1, &self.ptype, &self.p)
            if is_float:
                self._array = numpy.asarray(<double[:n:1]> self.p).reshape(shape)
            else:
                self._array = numpy.asarray(<int[:n:1]> self.p).reshape(shape)
        else:
            self._array = numpy.empty(shape, dtype=dtype)
        return self._array

    def __dealloc__(self):
        if self.ptype != NULL:
            c_deallocate_array(self.ptype)
        self.ptype = NULL
        self.p = NULL

cdef class Grid:
    cdef void* p
    cdef readonly int nx, ny, nz
    cdef readonly int nx_, ny_, nz_
    cdef readonly Domain domain
    cdef bint owned

    def __init__(self, Domain domain, int grid_type):
        self.domain = domain
        self.nx, self.ny, self.nz = domain.nx, domain.ny, domain.nz
        if grid_type == XGRID:
            self.nx += 1
            self.ny += 1
        self.p = domain_get_grid(domain.p, grid_type, 1, self.nx, 1, self.ny, 1, self.nz, domain.halox, domain.haloy, domain.haloz)
        self.owned = grid_type < 1 or grid_type > 4
        self.nx_, self.ny_, self.nz_ = self.nx + 2 * domain.halox, self.ny + 2 * domain.haloy, self.nz + 2 * domain.haloz
        domain.grids[grid_type] = self

    def __dealloc__(self):
        if self.owned and self.p != NULL:
            grid_finalize(self.p)

    def wrap(self, Array ar, bytes name, register=True):
        return ar.wrap_c_array(self.domain, 0, self.p, name, register)

def interp_x(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_x(<int>source.shape[2], <int>source.shape[1], <int>source.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_y(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_y(<int>source.shape[2], <int>source.shape[1], <int>source.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_z(double[:,:,::1] source not None, double[:,:,::1] target not None, int offset):
    grid_interp_z(<int>source.shape[2], <int>source.shape[1], <int>source.shape[0], <int>target.shape[0], &source[0,0,0], &target[0,0,0], offset)

def interp_xy(double[:,:,::1] source not None, double[:,:,::1] target not None, int ioffset, int joffset):
    grid_interp_xy(<int>source.shape[2], <int>source.shape[1], <int>target.shape[2], <int>target.shape[1], <int>source.shape[0], &source[0,0,0], &target[0,0,0], ioffset, joffset)

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

    def __dealloc__(self):
        if self.p != NULL:
            advection_finalize(self.p)

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
        advection_uv_calculate(1, <int>var._array.shape[0], self.p, self.grid.p, self.ugrid.p, <double *>u.p, pAh, timestep, self.ph, self.phu, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def v_3d(self, Array v not None, double timestep, Array var not None, Array Ah=None):
        cdef double * pAh = <double *>Ah.p if Ah is not None else NULL
        advection_uv_calculate(2, <int>var._array.shape[0], self.p, self.grid.p, self.vgrid.p, <double *>v.p, pAh, timestep, self.ph, self.phv, <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def w_3d(self, Array w not None, Array w_var not None, double timestep, Array var not None):
        advection_w_calculate(self.p, self.grid.p, <double *>w.p, <double *>w_var.p, timestep, self.ph, <double *>var.p)

cdef class VerticalDiffusion:
    cdef void* p
    cdef double cnpar

    def __init__(self, Grid grid, double cnpar=1.):
        self.p = vertical_diffusion_create(grid.p)
        self.cnpar = cnpar

    def __dealloc__(self):
        if self.p != NULL:
            vertical_diffusion_finalize(self.p)

    def prepare(self, Array nuh not None, double timestep, double molecular=0., bint use_ho=False):
        cdef Array mask = nuh.grid.mask
        cdef Array hn = nuh.grid.hn
        cdef Array ho = nuh.grid.ho if use_ho else hn
        c_vertical_diffusion_prepare(self.p, nuh.grid.nx_, nuh.grid.ny_, nuh.grid.nz_, molecular, <double *>nuh.p, timestep, self.cnpar, <int *>mask.p, <double *>ho.p, <double *>hn.p)

    def apply(self, Array var not None, Array ea2=None, Array ea4=None, bint use_ho=False):
        cdef double* pea2 = NULL
        cdef double* pea4 = NULL
        if ea2 is not None:
            pea2 = <double *>ea2.p
        if ea4 is not None:
            pea4 = <double *>ea4.p
        cdef Array mask = var.grid.mask
        cdef Array hn = var.grid.hn
        cdef Array ho = var.grid.ho if use_ho else hn
        c_vertical_diffusion_apply(self.p, var.grid.nx_, var.grid.ny_, var.grid.nz_, <int *>mask.p, <double *>ho.p, <double *>hn.p, <double *>var.p, pea2, pea4)

cdef class Momentum:
    cdef void* p
    cdef readonly Domain domain

    def __init__(self, Domain domain, int runtype, double Am0, double cnpar, int coriolis_scheme):
        self.domain = domain
        self.p = momentum_create(runtype, domain.p, Am0, cnpar, coriolis_scheme)

    def __dealloc__(self):
        if self.p != NULL:
            momentum_finalize(self.p)

    def momentum_diffusion_driver(self, Array h not None, Array hx not None, Array u not None, Array v not None, Array diffu not None, Array diffv not None):
        assert h.grid is self.domain.T
        assert hx.grid is self.domain.X
        assert u.grid is self.domain.U
        assert v.grid is self.domain.V
        assert diffu.grid is self.domain.U
        assert diffv.grid is self.domain.V
        cdef int nk = <int>h._array.shape[0] if h.z else 1
        momentum_diffusion_driver(self.p, nk, <double*> h.p, <double*> hx.p, <double*> u.p, <double*> v.p, <double*> diffu.p, <double*> diffv.p)

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

def bottom_friction(Array u not None, Array v not None, Array D not None, double avmmol, Array out not None, bint update_z0b):
    cdef int nx = u.grid.nx_
    cdef int ny = u.grid.ny_
    cdef Array mask = u.grid.mask
    cdef Array z0b = u.grid.z0b
    cdef Array z0b_min = u.grid.z0b_min
    assert not u.z
    assert v.grid is u.grid and not v.z
    assert D.grid is u.grid and not D.z
    assert out.grid is u.grid and not out.z
    c_bottom_friction(nx, ny, <int*>mask.p, <double*>u.p, <double*>v.p, <double*>D.p, <double*>z0b.p, <double*>z0b_min.p, avmmol, <double*>out.p, update_z0b)

def collect_3d_momentum_sources(Array dp, Array cor, Array adv, Array diff, Array idp, Array taus, Array rr, double dt, Array ea2, Array ea4):
    cdef Grid grid = dp.grid
    cdef Array mask = grid.mask
    cdef Array alpha = grid.alpha
    cdef Array ho = grid.ho
    cdef Array hn = grid.hn
    assert dp.grid is grid and not dp.z, 'dp'
    assert cor.grid is grid and cor.z == CENTERS, 'cor'
    assert adv.grid is grid and adv.z == CENTERS, 'adv'
    assert diff.grid is grid and diff.z == CENTERS, 'diff'
    assert idp.grid is grid and idp.z == CENTERS, 'idp'
    assert taus.grid is grid and not taus.z, 'taus'
    assert rr.grid is grid and not rr.z, 'rr'
    assert ea2.z == CENTERS
    assert ea4.z == CENTERS
    c_collect_3d_momentum_sources(grid.nx_, grid.ny_, grid.nz, grid.domain.halox, grid.domain.haloy,
        <int*>mask.p, <double*>alpha.p, <double*>ho.p, <double*>hn.p,
        <double*>dp.p, <double*>cor.p, <double*>adv.p, <double*>diff.p, <double*>idp.p, <double*>taus.p, <double*>rr.p, dt, <double*>ea2.p, <double*>ea4.p)

def advance_2d_transport(Array U, Array dp, Array taus, Array cor, Array adv, Array diff, Array damp, Array SA, Array SB, Array SD, Array SF, Array r, double dt):
    cdef Grid grid = U.grid
    cdef Array mask = grid.mask
    cdef Array alpha = grid.alpha
    cdef Array D = grid.D
    assert dp.grid is grid and not dp.z, 'dp'
    assert taus.grid is grid and not taus.z, 'taus'
    assert cor.grid is grid and not cor.z, 'cor'
    assert adv.grid is grid and not adv.z, 'adv'
    assert diff.grid is grid and not diff.z, 'diff'
    assert damp.grid is grid and not damp.z, 'damp'
    assert SA.grid is grid and not SA.z, 'SA'
    assert SB.grid is grid and not SB.z, 'SB'
    assert SD.grid is grid and not SD.z, 'SD'
    assert SF.grid is grid and not SF.z, 'SF'
    assert r.grid is grid and not r.z, 'rr'
    c_advance_2d_transport(grid.nx_, grid.ny_, grid.domain.halox, grid.domain.haloy, <int*>mask.p, <double*>alpha.p, <double*>D.p,
       <double*>dp.p, <double*>taus.p, <double*>cor.p, <double*>adv.p, <double*>diff.p, <double*>damp.p,
       <double*>SA.p, <double*>SB.p, <double*>SD.p, <double*>SF.p, <double*>r.p, dt, <double*>U.p)

def w_momentum_3d(Array pk, Array qk, double dt, Array w):
    cdef Grid grid = w.grid
    cdef Array mask = grid.mask
    assert w.z == INTERFACES, 'w'
    assert pk.grid is grid.ugrid and pk.z == CENTERS, 'pk'
    assert qk.grid is grid.vgrid and qk.z == CENTERS, 'qk'
    cdef Array dyu = pk.grid.dy
    cdef Array dxv = qk.grid.dx
    cdef Array iarea = grid.iarea
    cdef Array ho = grid.ho
    cdef Array hn = grid.hn
    c_w_momentum_3d(grid.nx_, grid.ny_, grid.nz, grid.domain.halox + 1, grid.domain.halox + grid.nx + 1, grid.domain.haloy + 1, grid.domain.haloy + grid.ny + 1,
        <int*>mask.p, <double*>dyu.p, <double*>dxv.p, <double*>iarea.p, <double*>ho.p, <double*>hn.p, <double*>pk.p, <double*>qk.p, dt, <double*>w.p)

def shear_frequency(Array uk, Array vk, Array SS):
    cdef Grid grid = SS.grid
    cdef Array mask = grid.mask
    assert SS.z == INTERFACES, 'w'
    assert uk.grid is grid.ugrid and uk.z == CENTERS, 'uk'
    assert vk.grid is grid.vgrid and vk.z == CENTERS, 'vk'
    cdef Array hu = uk.grid.hn
    cdef Array hv = vk.grid.hn
    c_shear_frequency(grid.nx_, grid.ny_, grid.nz, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny,
        <int*>mask.p, <double*>hu.p, <double*>hv.p, <double*>uk.p, <double*>vk.p, <double*>SS.p)

def shear_frequency2(Array uk, Array vk, Array num, Array SS):
    cdef Grid grid = SS.grid
    cdef Array mask = grid.mask
    assert SS.z == INTERFACES, 'w'
    assert uk.grid is grid.ugrid and uk.z == CENTERS, 'uk'
    assert vk.grid is grid.vgrid and vk.z == CENTERS, 'vk'
    assert num.grid is grid and num.z == INTERFACES, 'num'
    cdef Array h = grid.hn
    cdef Array hu = uk.grid.hn
    cdef Array hv = vk.grid.hn
    c_shear_frequency2(grid.nx_, grid.ny_, grid.nz, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny,
        <int*>mask.p, <double*>h.p, <double*>hu.p, <double*>hv.p, <double*>uk.p, <double*>vk.p, <double*>num.p, <double*>SS.p)

def surface_shear_velocity(Array taux, Array tauy, Array ustar):
    cdef Grid grid = ustar.grid
    assert not ustar.z, 'ustar'
    assert taux.grid is grid and not taux.z, 'taux'
    assert tauy.grid is grid and not tauy.z, 'tauy'
    cdef Array mask = grid.mask
    c_surface_shear_velocity(grid.nx_, grid.ny_, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny,
        <int*>mask.p, <double*>taux.p, <double*>tauy.p, <double*>ustar.p)

def bottom_shear_velocity(Array uk, Array vk, Array rru, Array rrv, Array ustar2_x, Array ustar2_y, Array ustar):
    cdef Grid grid = ustar.grid
    assert not ustar.z, 'ustar'
    assert uk.grid is grid.ugrid, 'uk'
    assert vk.grid is grid.vgrid, 'vk'
    assert rru.grid is grid.ugrid and not rru.z, 'rru'
    assert rrv.grid is grid.vgrid and not rrv.z, 'rrv'
    assert ustar2_x.grid is grid.ugrid and not ustar2_x.z, 'ustar2_x'
    assert ustar2_y.grid is grid.vgrid and not ustar2_y.z, 'ustar2_y'
    cdef Array mask = grid.mask
    cdef Array umask = grid.ugrid.mask
    cdef Array vmask = grid.vgrid.mask
    c_bottom_shear_velocity(grid.nx_, grid.ny_, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny,
        <int*>mask.p, <int*>umask.p, <int*>vmask.p, <double*>uk.p, <double*>vk.p, <double*>rru.p, <double*>rrv.p,
        <double*>ustar2_x.p, <double*>ustar2_y.p, <double*>ustar.p)

def advance_surface_elevation(double timestep, Array z, Array U, Array V, Array fwf):
    cdef Grid grid = z.grid
    cdef Array mask = grid.mask
    cdef Array dyu = U.grid.dy
    cdef Array dxv = V.grid.dx
    cdef Array iarea = grid.iarea
    assert not z.z, 'z'
    assert U.grid is grid.ugrid and not U.z, 'U'
    assert V.grid is grid.vgrid and not V.z, 'V'
    assert fwf.grid is grid and not fwf.z, 'fwf'
    c_advance_surface_elevation(grid.nx_, grid.ny_, grid.domain.halox, grid.domain.haloy, <int*>mask.p, <double*>dyu.p, <double*>dxv.p, <double*>iarea.p,
        <double*>z.p, <double*>U.p, <double*>V.p, <double*>fwf.p, timestep)

def surface_pressure_gradient(Array z, Array sp, Array dpdx, Array dpdy):
    cdef Grid grid = z.grid
    assert not z.z, 'z'
    assert sp.grid is grid and not sp.z, 'sp'
    assert dpdx.grid is grid.ugrid and not dpdx.z, 'dpdx'
    assert dpdy.grid is grid.vgrid and not dpdy.z, 'dpdy'
    cdef Array umask = dpdx.grid.mask
    cdef Array vmask = dpdy.grid.mask
    cdef Array idxu = dpdx.grid.idx
    cdef Array idyv = dpdy.grid.idy
    cdef Array H = grid.H
    cdef Array D = grid.D
    cdef double Dmin = grid.domain.Dmin
    c_surface_pressure_gradient(grid.nx_, grid.ny_, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny,
        <int*>umask.p, <int*>vmask.p, <double*>idxu.p, <double*>idyv.p,
        <double*>z.p, <double*>sp.p, <double*>H.p, <double*>D.p, Dmin, <double*>dpdx.p, <double*>dpdy.p)

def blumberg_mellor(Array buoy, Array idpdx, Array idpdy):
    cdef Grid grid = buoy.grid
    assert buoy.z == CENTERS, 'buoy'
    assert idpdx.grid is grid.ugrid and idpdx.z == CENTERS, 'idpdx'
    assert idpdy.grid is grid.vgrid and idpdy.z == CENTERS, 'idpdy'
    cdef Array umask = idpdx.grid.mask
    cdef Array vmask = idpdy.grid.mask
    cdef Array idxu = idpdx.grid.idx
    cdef Array idyv = idpdy.grid.idy
    cdef Array hu = idpdx.grid.hn
    cdef Array hv = idpdy.grid.hn
    cdef Array zf = grid.zf
    c_blumberg_mellor(grid.nx_, grid.ny_, grid.nz_, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny, <int*>umask.p, <int*>vmask.p, <double*>idxu.p, <double*>idyv.p, <double*>hu.p, <double*>hv.p, <double*>zf.p, <double*>buoy.p, <double*>idpdx.p, <double*>idpdy.p)

def shchepetkin_mcwilliams(Array buoy, Array idpdx, Array idpdy):
    cdef Grid grid = buoy.grid
    assert buoy.z == CENTERS, 'buoy'
    assert idpdx.grid is grid.ugrid and idpdx.z == CENTERS, 'idpdx'
    assert idpdy.grid is grid.vgrid and idpdy.z == CENTERS, 'idpdy'
    cdef Array mask = grid.mask
    cdef Array umask = idpdx.grid.mask
    cdef Array vmask = idpdy.grid.mask
    cdef Array idxu = idpdx.grid.idx
    cdef Array idyv = idpdy.grid.idy
    cdef Array h = grid.hn
    cdef Array z = grid.zin
    cdef Array zc = grid.zc
    c_shchepetkin_mcwilliams(grid.nx_, grid.ny_, grid.nz_, grid.domain.halox + 1, grid.domain.halox + grid.nx, grid.domain.haloy + 1, grid.domain.haloy + grid.ny, <int*>mask.p, <int*>umask.p, <int*>vmask.p, <double*>idxu.p, <double*>idyv.p, <double*>h.p, <double*>z.p, <double*>zc.p, <double*>buoy.p, <double*>idpdx.p, <double*>idpdy.p)

def multiply_add(double[::1] tgt, double[::1] add, double scale_factor):
    if tgt.shape[0] != 0:
        c_multiply_add(<int>tgt.shape[0], &tgt[0], &add[0], scale_factor)

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
            if subdomain_used(mask, max(0, xoffset + col * nx), min(<int>mask.shape[1], xoffset + (col + 1) * nx), max(0, yoffset + row * ny), min(<int>mask.shape[0], yoffset + (row + 1) * ny)):
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
            map[i] = subdomain_count_cells(mask, max(0, xoffset + col * nx), min(<int>mask.shape[1], xoffset + (col + 1) * nx), max(0, yoffset + row * ny), min(<int>mask.shape[0], yoffset + (row + 1) * ny))
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
