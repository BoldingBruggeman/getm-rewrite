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
cdef extern void domain_initialize(void* grid, int runtype, double Dmin, double* maxdt) nogil
cdef extern void domain_finalize(void* domain) nogil
cdef extern void domain_update_depths(void* domain) nogil
cdef extern void domain_do_vertical(void* domain) nogil
cdef extern void domain_tracer_bdy(void* domain, void* grid, int nz, double* field, int bdytype, double* bdy)
cdef extern void grid_interp_x(int nx, int ny, int nz, double* source, double* target, int ioffset) nogil
cdef extern void grid_interp_y(int nx, int ny, int nz, double* source, double* target, int joffset) nogil
cdef extern void grid_interp_z(int nx, int ny, int nz1, int nz2, double* source, double* target, int koffset) nogil
cdef extern void grid_interp_xy(int nx1, int ny1, int nx2, int ny2, int nz, double* source, double* target, int ioffset, int joffset) nogil
cdef extern void get_array(int source_type, void* grid, const char* name, int* grid_type, int* sub_type, int* data_type, void** p) nogil
cdef extern void* advection_create(int scheme, void* tgrid) nogil
cdef extern void advection_2d_calculate(int direction, void* advection, void* tgrid, void* ugrid, double* pu, double Ah, double timestep, double* pD, double* pDU, double* pvar) nogil
cdef extern void advection_w_calculate(void* padvection, void* tgrid, double* pw, double* pw_var, double timestep, double* ph, double* pvar)
cdef extern void* vertical_diffusion_create(void* tgrid) nogil
cdef extern void vertical_diffusion_calculate(void* diffusion, void* tgrid, double molecular, double* pnuh, double timestep, double cnpar, double* pho, double* phn, double* pvar, double* pea2, double* pea4) nogil
cdef extern void* momentum_create(int runtype, void* pdomain, int apply_bottom_friction) nogil
cdef extern void momentum_u_2d(int direction, void* momentum, double timestep, double* ptausx, double* pdpdx) nogil
cdef extern void momentum_u_3d(int direction, void* momentum, double timestep, double* ptausx, double* pdpdx, double* pidpdx, double* pviscosity) nogil
cdef extern void momentum_w_3d(void* momentum, double timestep) nogil
cdef extern void momentum_uv_coriolis(int direction, void* momentum) nogil
cdef extern void momentum_uv_coriolis_3d(int direction, void* momentum) nogil
cdef extern void momentum_bottom_friction_2d(void* momentum, int runtype) nogil
cdef extern void momentum_bottom_friction_3d(void* momentum) nogil
cdef extern void momentum_shear_frequency(void* momentum, double* pviscosity) nogil
cdef extern void* pressure_create(int runtype, void* pdomain) nogil
cdef extern void pressure_surface(void* pressure, double* pz, double* psp) nogil
cdef extern void* sealevel_create(void* pdomain) nogil
cdef extern void sealevel_update(void* sealevel, double timestep, double* pU, double* pV, double* pfwf) nogil
#cdef extern void sealevel_update_uvx(void* sealevel) nogil
cdef extern void sealevel_boundaries(void* sealevel, void* momentum, double timestep) nogil
cdef extern void c_exponential_profile_1band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* top, double* out) nogil
cdef extern void c_exponential_profile_1band_centers(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* k, double* top, double* out) nogil
cdef extern void c_exponential_profile_2band_interfaces(int nx, int ny, int nz, int istart, int istop, int jstart, int jstop, int* mask, double* h, double* f, double* k1, double* k2, double* top, double* out) nogil

cpdef enum:
    TGRID = 1
    UGRID = 2
    VGRID = 3
    XGRID = 4
    UUGRID = -1
    VVGRID = -2
    UVGRID = -3
    VUGRID = -4

# JB: for now, this is copy of the constants in parallel. Values must saty in sync
cpdef enum:
    TOP_BOTTOM = 1
    LEFT_RIGHT = 2

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

    cdef wrap_c_array(self, Domain domain, int source, void* obj, bytes name):
        cdef int grid_type
        cdef int sub_type
        cdef int data_type
        self.on_boundary = False
        get_array(source, obj, name, &grid_type, &sub_type, &data_type, &self.p)
        if self.p == NULL:
            return
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
        self.register()
        return self

    def wrap_ndarray(self, numpy.ndarray data not None, on_boundary=False):
        self.on_boundary = on_boundary
        if self.on_boundary:
            assert data.ndim in (1, 2) and data.flags['C_CONTIGUOUS'], 'Invalid array properties for wrapping: %i dimensions, flags %s' % (data.ndim, data.flags)
            assert data.shape[0] == self.grid.nbdyp, 'Incorrect shape of first dimension (number of boundary points): expected %i, got %i' % (self.grid.nbdyp, data.shape[0])
            assert data.ndim == 1 or data.shape[1] == self.grid.nz_ or data.shape[1] == self.grid.nz_ + 1, 'Incorrect shape of second dimension (number of layers): expected %i or %i, got %i' % (self.grid.nz_, self.grid.nz_ + 1, data.shape[1])
        elif data.ndim != 0:
            assert data.ndim in (2, 3) and data.flags['C_CONTIGUOUS'], 'Invalid array properties for wrapping: %i dimensions, flags %s' % (data.ndim, data.flags)
            assert data.shape[data.ndim - 1] == self.grid.nx_ and data.shape[data.ndim - 2] == self.grid.ny_, 'Incorrect horizontal extent: expected (ny=%i,nx=%i), got (ny=%i,nx=%i)' % (self.grid.ny_, self.grid.nx_, data.shape[data.ndim - 2], data.shape[data.ndim - 1])
            assert data.ndim == 2 or data.shape[0] == self.grid.nz_ or data.shape[0] == self.grid.nz_ + 1, 'Incorrect vertical extent: expected %i or %i, got %i' % (self.grid.nz_, self.grid.nz_ + 1, data.shape[0])
        self._array = data
        self.p = self._array.data
        self.finish_initialization()

    def update_boundary(self, int bdytype, Array bdy=None):
        assert not self.on_boundary, 'update_boundary cannot be called on boundary arrays.'
        assert self.ndim == 2 or bdy is None or self._array.shape[0] == bdy._array.shape[1]
        domain_tracer_bdy(self.grid.domain.p, self.grid.p, 1 if self.ndim == 2 else self._array.shape[0], <double*>self.p, bdytype, NULL if bdy is None else <double*>bdy.p)

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

    def wrap(self, Array ar, bytes name):
        return ar.wrap_c_array(self.domain, 0, self.p, name)

    def interp_x(self, Array source not None, Array target not None, int offset):
        grid_interp_x(source._array.shape[source._array.ndim - 1], source._array.shape[source._array.ndim - 2], 1 if source._array.ndim == 2 else source._array.shape[0], <double *>source.p, <double *>target.p, offset)

    def interp_y(self, Array source not None, Array target not None, int offset):
        grid_interp_y(source._array.shape[source._array.ndim - 1], source._array.shape[source._array.ndim - 2], 1 if source._array.ndim == 2 else source._array.shape[0], <double *>source.p, <double *>target.p, offset)

    def interp_z(self, Array source not None, Array target not None, int offset):
        grid_interp_z(source._array.shape[2], source._array.shape[1], source._array.shape[0], target._array.shape[0], <double *>source.p, <double *>target.p, offset)

    def interp_xy(self, Array source not None, Array target not None, int ioffset, int joffset):
        grid_interp_xy(source._array.shape[source._array.ndim - 1], source._array.shape[source._array.ndim - 2], target._array.shape[target._array.ndim - 1], target._array.shape[target._array.ndim - 2], 1 if source._array.ndim == 2 else source._array.shape[0], <double *>source.p, <double *>target.p, ioffset, joffset)

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

    def initialize(self, int runtype, double Dmin):
        domain_initialize(self.p, runtype, Dmin, &self.maxdt)

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
    cdef readonly numpy.ndarray D, h
    cdef double[:, :, ::1] h_work, hu, hv
    cdef double[:, ::1] D_ref, DU, DV

    def __init__(self, Grid grid, int scheme):
        self.tgrid = grid
        self.ugrid = grid.ugrid
        self.vgrid = grid.vgrid
        self.p = advection_create(scheme, self.tgrid.p)

        # Work array for layer heights that change during 3d advection
        # This is shared by the temeprary column height that evolves similarly during 2d advection
        self.h = numpy.empty((self.tgrid.nz_, self.tgrid.ny_, self.tgrid.nx_))
        self.D = self.h[0,...]

        # Store references to column height and layer thicknesses
        self.D_ref = self.tgrid.D.all_values
        self.DU = self.ugrid.D.all_values
        self.DV = self.vgrid.D.all_values
        self.hu = self.ugrid.hn.all_values
        self.hv = self.vgrid.hn.all_values

        # Efficient reference to work array
        self.h_work = self.h

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def __call__(self, Array u not None, Array v not None, double timestep, Array var not None, double Ah = 0, skip_initial_halo_exchange=True):
        assert u.grid is self.ugrid
        assert v.grid is self.vgrid
        assert var.grid is self.tgrid
        self.h_work[0,:,:] = self.D_ref
        if not skip_initial_halo_exchange:
            var.update_halos(TOP_BOTTOM)
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, Ah, 0.5 * timestep, &self.h_work[0,0,0], &self.DV[0,0], <double *>var.p)
        var.update_halos(LEFT_RIGHT)
        advection_2d_calculate(1, self.p, self.tgrid.p, self.ugrid.p, <double *>u.p, Ah, timestep, &self.h_work[0,0,0], &self.DU[0,0], <double *>var.p)
        var.update_halos(TOP_BOTTOM)
        advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, <double *>v.p, Ah, 0.5 * timestep, &self.h_work[0,0,0], &self.DV[0,0], <double *>var.p)

    @cython.initializedcheck(False)
    @cython.boundscheck(False) # turn off bounds-checking for entire function
    def apply_3d(self, Array u not None, Array v not None, Array w not None, double timestep, Array var not None, double Ah = 0, new_h=False, skip_initial_halo_exchange=False, Array w_var=None):
        if w_var is None:
            w_var = w
        assert u.grid is self.ugrid, 'grid mismatch for u: expected %s, got %s' % (self.ugrid.postfix, u.grid.postfix)
        assert v.grid is self.vgrid, 'grid mismatch for v: expected %s, got %s' % (self.vgrid.postfix, v.grid.postfix)
        assert w.grid is self.tgrid, 'grid mismatch for w: expected %s, got %s' % (self.tgrid.postfix, w.grid.postfix)
        assert w.z == INTERFACES, 'grid mismatch for w: expected values at layer interfaces'
        assert w_var.grid is self.tgrid, 'grid mismatch for w_var: expected %s, got %s' % (self.tgrid.postfix, w_var.grid.postfix)
        assert w_var.z == INTERFACES, 'grid mismatch for w_var: expected values at layer interfaces'
        assert var.grid is self.tgrid, 'grid mismatch for advected quantity: expected %s, got %s' % (self.tgrid.postfix, var.grid.postfix)
        assert var.ndim == 3 and var.z != INTERFACES, 'grid mismatch for advected quantity: expected 3D variable defined at layer centers'
        cdef double[:, :, ::1] avar, au, av, h
        cdef int k
        avar = <double[:var.grid.nz_, :var.grid.ny_, :var.grid.nx_:1]> var.p
        au = <double[:u.grid.nz_, :u.grid.ny_, :u.grid.nx_:1]> u.p
        av = <double[:v.grid.nz_, :v.grid.ny_, :v.grid.nx_:1]> v.p
        h = (self.tgrid.hn if new_h else self.tgrid.ho).all_values
        self.h_work[:,:,:] = h
        if not skip_initial_halo_exchange:
            var.update_halos(LEFT_RIGHT)
        for k in range(var._array.shape[0]):
            advection_2d_calculate(1, self.p, self.tgrid.p, self.ugrid.p, &au[k,0,0], Ah, 0.5 * timestep, &self.h_work[k,0,0], &self.hu[k,0,0], &avar[k,0,0])
        var.update_halos(TOP_BOTTOM)
        for k in range(var._array.shape[0]):
            advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, &av[k,0,0], Ah, 0.5 * timestep, &self.h_work[k,0,0], &self.hv[k,0,0], &avar[k,0,0])
        advection_w_calculate(self.p, self.tgrid.p, <double *>w.p, <double *>w_var.p, timestep, &self.h_work[0,0,0], <double *>var.p)
        var.update_halos(TOP_BOTTOM)
        for k in range(var._array.shape[0]):
            advection_2d_calculate(2, self.p, self.tgrid.p, self.vgrid.p, &av[k,0,0], Ah, 0.5 * timestep, &self.h_work[k,0,0], &self.hv[k,0,0], &avar[k,0,0])
        var.update_halos(LEFT_RIGHT)
        for k in range(var._array.shape[0]):
            advection_2d_calculate(1, self.p, self.tgrid.p, self.ugrid.p, &au[k,0,0], Ah, 0.5 * timestep, &self.h_work[k,0,0], &self.hu[k,0,0], &avar[k,0,0])

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

    def __call__(self, Array nuh not None, double timestep, Array var not None, double molecular=0., Array ea2=None, Array ea4=None):
        cdef double* pea2 = NULL
        cdef double* pea4 = NULL
        if ea2 is not None:
            pea2 = <double *>ea2.p
        if ea4 is not None:
            pea4 = <double *>ea4.p
        vertical_diffusion_calculate(self.p, self.tgrid.p, molecular, <double *>nuh.p, timestep, self.cnpar, <double *>self.hn.p, <double *>self.hn.p, <double *>var.p, pea2, pea4)

cdef class Simulation:
    cdef readonly Domain domain
    cdef readonly int runtype
    cdef void* pmomentum
    cdef void* ppressure
    cdef void* psealevel
    cdef readonly int nx, ny
    cdef int apply_bottom_friction
    cdef int ufirst, u3dfirst

    def __init__(self, Domain domain, int runtype, int apply_bottom_friction, int ufirst=False):
        self.domain = domain
        domain.initialize(runtype)
        self.runtype = runtype
        self.pmomentum = momentum_create(runtype, domain.p, apply_bottom_friction)
        self.ppressure = pressure_create(runtype, domain.p)
        self.psealevel = sealevel_create(domain.p)
        self.nx, self.ny = domain.nx, domain.ny
        self.apply_bottom_friction = apply_bottom_friction
        self.ufirst = ufirst
        self.u3dfirst = ufirst

    def uv_momentum_2d(self, double timestep, Array tausx not None, Array tausy not None, Array dpdx not None, Array dpdy not None):
        assert tausx.grid is self.domain.U, 'grid mismatch for tausx: expected %s, got %s' % (self.domain.U.postfix, tausx.grid.postfix)
        assert tausy.grid is self.domain.V, 'grid mismatch for tausy: expected %s, got %s' % (self.domain.V.postfix, tausy.grid.postfix)
        assert dpdx.grid is self.domain.U, 'grid mismatch for dpdx: expected %s, got %s' % (self.domain.U.postfix, dpdx.grid.postfix)
        assert dpdy.grid is self.domain.V, 'grid mismatch for dpdy: expected %s, got %s' % (self.domain.V.postfix, dpdy.grid.postfix)
        if self.apply_bottom_friction:
            momentum_bottom_friction_2d(self.pmomentum, self.runtype)
        if self.ufirst:
            momentum_u_2d(1, self.pmomentum, timestep, <double *>tausx.p, <double *>dpdx.p)
            self.U.update_halos()
            momentum_uv_coriolis(1, self.pmomentum)
            momentum_u_2d(2, self.pmomentum, timestep, <double *>tausy.p, <double *>dpdy.p)
            self.V.update_halos()
            momentum_uv_coriolis(2, self.pmomentum)
        else:
            momentum_u_2d(2, self.pmomentum, timestep, <double *>tausy.p, <double *>dpdy.p)
            self.V.update_halos()
            momentum_uv_coriolis(2, self.pmomentum)
            momentum_u_2d(1, self.pmomentum, timestep, <double *>tausx.p, <double *>dpdx.p)
            self.U.update_halos()
            momentum_uv_coriolis(1, self.pmomentum)
        self.ufirst = not self.ufirst

    def uvw_momentum_3d(self, double timestep, Array tausx not None, Array tausy not None, Array dpdx not None, Array dpdy not None, Array idpdx not None, Array idpdy not None, Array viscosity_u not None, Array viscosity_v not None):
        assert tausx.grid is self.domain.U, 'grid mismatch for tausx: expected %s, got %s' % (self.domain.U.postfix, tausx.grid.postfix)
        assert tausy.grid is self.domain.V, 'grid mismatch for tausy: expected %s, got %s' % (self.domain.V.postfix, tausy.grid.postfix)
        assert dpdx.grid is self.domain.U, 'grid mismatch for dpdx: expected %s, got %s' % (self.domain.U.postfix, dpdx.grid.postfix)
        assert dpdy.grid is self.domain.V, 'grid mismatch for dpdy: expected %s, got %s' % (self.domain.V.postfix, dpdy.grid.postfix)
        assert idpdx.grid is self.domain.U, 'grid mismatch for idpdx: expected %s, got %s' % (self.domain.U.postfix, idpdx.grid.postfix)
        assert idpdy.grid is self.domain.V, 'grid mismatch for idpdy: expected %s, got %s' % (self.domain.V.postfix, idpdy.grid.postfix)
        assert viscosity_u.grid is self.domain.U and viscosity_u.z == INTERFACES, 'grid mismatch for viscosity_u: expected %s, got %s' % (self.domain.U.postfix, viscosity_u.grid.postfix)
        assert viscosity_v.grid is self.domain.V and viscosity_v.z == INTERFACES, 'grid mismatch for viscosity_v: expected %s, got %s' % (self.domain.V.postfix, viscosity_v.grid.postfix)

        if self.apply_bottom_friction:
            momentum_bottom_friction_3d(self.pmomentum)

        if self.u3dfirst:
            momentum_u_3d(1, self.pmomentum, timestep, <double *>tausx.p, <double *>dpdx.p, <double *>idpdx.p, <double *>viscosity_u.p)
            self.pk.update_halos()
            momentum_uv_coriolis_3d(1, self.pmomentum)
            momentum_u_3d(2, self.pmomentum, timestep, <double *>tausy.p, <double *>dpdy.p, <double *>idpdy.p, <double *>viscosity_v.p)
            self.qk.update_halos()
            momentum_uv_coriolis_3d(2, self.pmomentum)
        else:
            momentum_u_3d(2, self.pmomentum, timestep, <double *>tausy.p, <double *>dpdy.p, <double *>idpdy.p, <double *>viscosity_v.p)
            self.qk.update_halos()
            momentum_uv_coriolis_3d(2, self.pmomentum)
            momentum_u_3d(1, self.pmomentum, timestep, <double *>tausx.p, <double *>dpdx.p, <double *>idpdx.p, <double *>viscosity_u.p)
            self.pk.update_halos()
            momentum_uv_coriolis_3d(1, self.pmomentum)
        self.ufirst = not self.u3dfirst

        momentum_w_3d(self.pmomentum, timestep)

#   call self%velocities_3d()
#   call self%uv_advection_3d(dt)
##if 0
#   call self%uv_diffusion_3d(dt) !KB - makes model go wrong
##else
#if (associated(self%logs)) call self%logs%info('*** missing uv_diffusion_3d() ***',level=0)
##endif
#   call self%shear_frequency(viscosity)
#   call self%stresses(tausx,tausy)

    def update_surface_pressure_gradient(self, Array z not None, Array sp not None):
        assert z.grid is self.domain.T
        assert sp.grid is self.domain.T
        pressure_surface(self.ppressure, <double *>z.p, <double *>sp.p)

    def update_sealevel(self, double timestep, Array U not None, Array V not None, Array fwf not None):
        assert U.grid is self.domain.U
        assert V.grid is self.domain.V
        sealevel_update(self.psealevel, timestep, <double *>U.p, <double *>V.p, <double *>fwf.p)

    #def update_sealevel_uvx(self):
    #    sealevel_update_uvx(self.psealevel)

    def update_sealevel_boundaries(self, double timestep):
        sealevel_boundaries(self.psealevel, self.pmomentum, timestep)

    def update_shear_frequency(self, Array viscosity not None):
        assert viscosity.grid is self.domain.T and viscosity.z == INTERFACES, 'grid mismatch for viscosity: expected %s, got %s' % (self.domain.T.postfix, viscosity.grid.postfix)
        momentum_shear_frequency(self.pmomentum, <double *>viscosity.p)

    def wrap(self, Array ar not None, bytes name, int source):
        cdef void* obj = self.pmomentum
        if (source == 2):
            obj = self.ppressure
        elif (source == 3):
            obj = self.psealevel
        return ar.wrap_c_array(self.domain, source, obj, name)

def exponential_profile_1band_interfaces(Array mask not None, Array h not None, Array k not None, Array top not None, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is k.grid and not k.z
    assert mask.grid is top.grid and not top.z
    assert mask.grid is out.grid and out.z == INTERFACES
    c_exponential_profile_1band_interfaces(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>k.p, <double *>top.p, <double *>out.p)

def exponential_profile_1band_centers(Array mask not None, Array h not None, Array k not None, Array top not None, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is k.grid and not k.z
    assert mask.grid is top.grid and not top.z
    assert mask.grid is out.grid and out.z == CENTERS
    c_exponential_profile_1band_centers(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>k.p, <double *>top.p, <double *>out.p)

def exponential_profile_2band_interfaces(Array mask not None, Array h not None, Array f not None, Array k1 not None, Array k2 not None, Array top not None, Array out=None):
    assert mask.grid is h.grid and h.z == CENTERS
    assert mask.grid is f.grid and not f.z
    assert mask.grid is k1.grid and not k1.z
    assert mask.grid is k2.grid and not k2.z
    assert mask.grid is top.grid and not top.z
    assert mask.grid is out.grid and out.z == INTERFACES
    c_exponential_profile_2band_interfaces(mask.grid.nx_, mask.grid.ny_, mask.grid.nz_, mask.grid.domain.halox, mask.grid.nx_ - mask.grid.domain.halox, mask.grid.domain.haloy, mask.grid.ny_ - mask.grid.domain.haloy, <int *>mask.p, <double *>h.p, <double *>f.p, <double *>k1.p, <double *>k2.p, <double *>top.p, <double *>out.p)

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
cdef int decomposition_is_valid(const int[:, ::1] mask, int nx, int ny, int xoffset, int yoffset, int ncpus) nogil:
    cdef int nrow = <int>ceil((mask.shape[0] - yoffset) / float(ny))
    cdef int ncol = <int>ceil((mask.shape[1] - xoffset) / float(nx))
    cdef int free_cpus = ncpus
    cdef int row, col, i, j, n
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
cdef int[:, ::1] get_map(const int[:, ::1] mask, int nx, int ny, int xoffset, int yoffset):
    cdef int nrow = <int>ceil((mask.shape[0] - yoffset) / float(ny))
    cdef int ncol = <int>ceil((mask.shape[1] - xoffset) / float(nx))
    cdef int[:, ::1] map
    cdef int row, col
    map = numpy.empty((nrow, ncol), dtype=numpy.intc)
    for row in range(nrow):
        for col in range(ncol):
            map[row, col] = subdomain_count_cells(mask, max(0, xoffset + col * nx), min(mask.shape[1], xoffset + (col + 1) * nx), max(0, yoffset + row * ny), min(mask.shape[0], yoffset + (row + 1) * ny))
    return map

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.profile(False)
cdef int get_from_map(const int[:, ::1] map, int row, int col) nogil:
    if row < 0 or col < 0 or row >= map.shape[0] or col >= map.shape[1]:
        return 0
    return map[row, col]

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
cdef int get_cost(const int[:, ::1] map, int nx, int ny) nogil:
    cdef int row, col, nint, nout, max_cost
    cdef int halo = 2
    max_cost = 0
    for row in range(map.shape[0]):
        for col in range(map.shape[1]):
            nint = map[row, col]  # unmasked cells only
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
    return max_cost + nx * ny     # add an overhead for all cells - unmasked or not

def find_subdiv_solutions(const int[:, ::1] mask not None, int nx, int ny, int ncpus):
    cdef int[:, ::1] map
    cost = -1
    solution = None
    for yoffset in range(1 - ny, 1):
        for xoffset in range(1 - nx, 1):
            if decomposition_is_valid(mask, nx, ny, xoffset, yoffset, ncpus):
                map = get_map(mask, nx, ny, xoffset, yoffset)
                current_cost = get_cost(map, nx, ny)
                if cost == -1 or current_cost < cost:
                    cost = current_cost
                    solution = (xoffset, yoffset, cost, numpy.asarray(map))
    return solution
