import os
import sys
import ctypes

import numpy

from . import parallel

# Determine potential names of FABM dynamic library.
if os.name == 'nt':
   dllpaths = ('_pygetm.dll', 'lib_pygetm.dll')
elif os.name == 'posix' and sys.platform == 'darwin':
   dllpaths = ('lib_pygetm.dylib',)
else:
   dllpaths = ('lib_pygetm.so',)

def find_library(basedir):
    for dllpath in dllpaths:
        dllpath = os.path.join(basedir, dllpath)
        if os.path.isfile(dllpath):
            return dllpath

# Find dynamic library.
# Look first in local directory, then in Python path.
if 'PYGETM_DLL' in os.environ:
    dllpath = os.environ['PYGETM_DLL']
    print('Using _pygetm library %s as set by environment variable PYGETM_DLL' % dllpath)
    assert os.path.isfile(dllpath), '%s not found' % dllpath
else:
    dllpath = find_library(os.path.dirname(os.path.abspath(__file__)))
    if not dllpath:
        for basedir in sys.path:
            dllpath = find_library(basedir)
            if dllpath:
                break

if not dllpath:
   print('Unable to locate pygetm dynamic library %s.' % (' or '.join(dllpaths),))
   sys.exit(1)

# Load FABM library.
_pygetm = ctypes.CDLL(str(dllpath))

_pygetm.domain_create.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]
_pygetm.domain_create.restype = ctypes.c_void_p

_pygetm.domain_get_grid.argtypes = [ctypes.c_void_p, ctypes.c_char]
_pygetm.domain_get_grid.restype = ctypes.c_void_p

_pygetm.grid_get_array.argtypes = [ctypes.c_void_p,  ctypes.c_char_p]
_pygetm.grid_get_array.restype = ctypes.c_void_p

_pygetm.domain_initialize.argtypes = [ctypes.c_void_p, ctypes.c_int]
_pygetm.domain_initialize.restype = None

_pygetm.domain_depth_update.argtypes = [ctypes.c_void_p]
_pygetm.domain_depth_update.restype = None

_pygetm.advection_create.argtypes = []
_pygetm.advection_create.restype = ctypes.c_void_p

CONTIGUOUS = str('CONTIGUOUS')
arrtype2D = numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, flags=CONTIGUOUS)

_pygetm.advection_calculate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, arrtype2D, arrtype2D, ctypes.c_double, arrtype2D]
_pygetm.advection_calculate.restype = None

_pygetm.momentum_create.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p]
_pygetm.momentum_create.restype = ctypes.c_void_p

_pygetm.momentum_get_array.argtypes = [ctypes.c_void_p,  ctypes.c_char_p]
_pygetm.momentum_get_array.restype = ctypes.c_void_p

_pygetm.momentum_uv_momentum_2d.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_double, arrtype2D, arrtype2D, arrtype2D, arrtype2D]
_pygetm.momentum_uv_momentum_2d.restype = None

_pygetm.pressure_create.argtypes = [ctypes.c_int, ctypes.c_void_p]
_pygetm.pressure_create.restype = ctypes.c_void_p

_pygetm.pressure_surface.argtypes = [ctypes.c_void_p, arrtype2D, arrtype2D]
_pygetm.pressure_surface.restype = None

_pygetm.pressure_get_array.argtypes = [ctypes.c_void_p,  ctypes.c_char_p]
_pygetm.pressure_get_array.restype = ctypes.c_void_p

_pygetm.sealevel_create.argtypes = [ctypes.c_void_p]
_pygetm.sealevel_create.restype = ctypes.c_void_p

_pygetm.sealevel_update.argtypes = [ctypes.c_void_p, ctypes.c_double, arrtype2D, arrtype2D]
_pygetm.sealevel_update.restype = None

class FortranObject:
    def get_arrays(self, get_function=None, names=(), default_shape=None, default_dtype=ctypes.c_double, shapes={}, dtypes={}, halo=None):
        for name in names:
            p = get_function(self.p, name.encode('ascii'))
            assert p, 'Null pointer returned for %s' % name

            # Cast returned pointer to teh correct shape and dtype.
            shape = shapes.get(name, default_shape)
            dtype = dtypes.get(name, default_dtype)
            for l in shape[::-1]:
                dtype = dtype * l
            p = ctypes.cast(p, ctypes.POINTER(dtype)).contents

            # Cast to numpy array and set byte order to "native" rather than the explicit little/big endian.
            # The latter is needed to allow exchange of this field with mpi4py
            data_ = numpy.ctypeslib.as_array(p).newbyteorder('=')

            # Store field with and without halo as class attribute
            setattr(self, name + '_', data_)
            setattr(self, name, data_[tuple([slice(halo, -halo)] * len(shape))])

class Grid(FortranObject):
    def __init__(self, domain, grid_type='T'):
        self.p = _pygetm.domain_get_grid(domain.p, grid_type.encode('ascii'))
        self.halo = domain.halo
        self.domain = domain
        self.get_arrays(_pygetm.grid_get_array, ('c1', 'c2', 'x', 'y', 'dx', 'dy', 'lon', 'lat', 'dlon', 'dlat', 'H', 'D', 'mask', 'z', 'area', 'inv_area', 'cor'), default_shape=domain.shape[1:], shapes={'c1': (domain.shape[-1],), 'c2': (domain.shape[-2],)}, dtypes={'mask': ctypes.c_int}, halo=domain.halo)

    def array(self, fill=None, dtype=float):
        data = numpy.empty(self.H_.shape, dtype=dtype)
        if fill is not None:
            data[...] = fill
        return data[self.halo:-self.halo, self.halo:-self.halo], data

    def update_halos(self):
        self.domain.distribute(self.dx_).update_halos()
        self.domain.distribute(self.dy_).update_halos()
        self.domain.distribute(self.mask_).update_halos()
        self.domain.distribute(self.H_).update_halos()

def find_interfaces(c):
    c_if = numpy.empty((c.size + 1),)
    d = numpy.diff(c)
    c_if[1:-1] = c[:-1] + 0.5 * d
    c_if[0] = c[0] - 0.5 * d[0]
    c_if[-1] = c[-1] + 0.5 * d[-1]
    return c_if

class Simulation:
    def __init__(self, domain, runtype, advection_scheme=4):
        assert not domain.initialized
        self.runtype = runtype
        self.domain = domain
        self.domain.initialize(self.runtype)
        self.advection = Advection(self.domain, advection_scheme)
        self.pressure = Pressure(self.runtype, self.domain)
        self.momentum = Momentum(self.runtype, self.domain, self.advection)
        self.sealevel = Sealevel(self.domain)

class Domain:
    @staticmethod
    def create_cartesian(x, y, nlev, H=0., **kwargs):
        assert x.ndim == 1, 'x coordinate must be one-dimensional'
        assert y.ndim == 1, 'y coordinate must be one-dimensional'
        assert nlev >= 1, 'number of levels must be >= 1'
        tiling = None if not kwargs else parallel.Tiling(1, 1, **kwargs)
        domain = Domain(1, nlev, 1, y.size, 1, x.size, tiling=tiling)
        domain.T.c1[:] = x
        domain.T.c2[:] = y
        domain.T.x[:, :] = x[numpy.newaxis, :]
        domain.T.y[:, :] = y[:, numpy.newaxis]
        domain.T.H[...] = H
        domain.T.xi = numpy.broadcast_to(find_interfaces(x)[numpy.newaxis, :], (y.size + 1, x.size + 1))
        domain.T.yi = numpy.broadcast_to(find_interfaces(y)[:, numpy.newaxis], (y.size + 1, x.size + 1))
        return domain

    @staticmethod
    def partition(tiling, nx, ny, nlev, global_domain, runtype):
        assert nx % tiling.ncol == 0
        assert ny % tiling.nrow == 0
        assert global_domain is None or global_domain.initialized

        domain = Domain(1, nlev, 1, ny // tiling.nrow, 1, nx // tiling.ncol, tiling=tiling)

        # Scatter coordinates, bathymetry and mask to subdomains
        domain.distribute(domain.T.x_).scatter(None if global_domain is None else global_domain.T.x)
        domain.distribute(domain.T.y_).scatter(None if global_domain is None else global_domain.T.y)
        domain.distribute(domain.T.H_).scatter(None if global_domain is None else global_domain.T.H)
        domain.distribute(domain.T.mask_).scatter(None if global_domain is None else global_domain.T.mask)

        # Extract 1D coordinate vectors per subdomain
        domain.T.c1_[:] = domain.T.x_[domain.halo, :]
        domain.T.c2_[:] = domain.T.y_[:, domain.halo]

        domain.initialize(runtype)

        return domain

    def __init__(self, kmin, kmax, jmin, jmax, imin, imax, tiling=None):
        self.imin, self.imax = imin, imax
        self.jmin, self.jmax = jmin, jmax
        self.kmin, self.kmax = kmin, kmax
        halox, haloy, haloz = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
        self.p = _pygetm.domain_create(imin, imax, jmin, jmax, kmin, kmax, halox, haloy, haloz)
        self.halo = halox.value
        self.shape = (kmax - kmin + 1, jmax - jmin + 1 + 2 * self.halo, imax - imin + 1 + 2 * self.halo)
        self.tiling = tiling
        self.T = Grid(self, 'T')
        self.U = Grid(self, 'U')
        self.V = Grid(self, 'V')
        self.X = Grid(self, 'X')
        if self.tiling:
            self.dist_z = self.distribute(self.T.z_)
        self.initialized = False

    def initialize(self, runtype):
        _pygetm.domain_initialize(self.p, runtype)
        if self.tiling:
            self.T.update_halos()
            self.U.update_halos()
            self.V.update_halos()
            self.X.update_halos()
            self.depth_update()
        self.initialized = True

    def depth_update(self):
        if self.tiling:
            self.dist_z.update_halos()
        _pygetm.domain_depth_update(self.p)
        if self.tiling:
            self.distribute(self.U.D_).update_halos()
            self.distribute(self.V.D_).update_halos()
            self.distribute(self.X.D_).update_halos()

    def distribute(self, field):
        return self.tiling.wrap(field, halo=self.halo)

class Advection:
    def __init__(self, domain, scheme):
        self.p = _pygetm.advection_create()
        self.pdomain = domain.p
        self.scheme = scheme

    def calculate(self, u, v, timestep, var):
        _pygetm.advection_calculate(self.p, self.scheme, self.pdomain, u, v, timestep, var)

class Momentum(FortranObject):
    def __init__(self, runtype, domain, advection):
        self.runtype = runtype
        self.p = _pygetm.momentum_create(self.runtype, domain.p, advection.p)
        self.get_arrays(_pygetm.momentum_get_array, ('U', 'V'), default_shape=domain.shape[1:], halo=domain.halo)

    def uv_momentum_2d(self, timestep, tausx, tausy, dpdx, dpdy):
        _pygetm.momentum_uv_momentum_2d(self.p, self.runtype, timestep, tausx, tausy, dpdx, dpdy)

class Pressure(FortranObject):
    def __init__(self, runtype, domain):
        self.p = _pygetm.pressure_create(runtype, domain.p)
        self.get_arrays(_pygetm.pressure_get_array, ('dpdx', 'dpdy'), default_shape=domain.shape[1:], halo=domain.halo)

    def surface(self, z, sp):
        _pygetm.pressure_surface(self.p, z, sp)

class Sealevel:
    def __init__(self, domain):
        self.p = _pygetm.sealevel_create(domain.p)

    def update(self, timestep, U, V):
        _pygetm.sealevel_update(self.p, timestep, U, V)

