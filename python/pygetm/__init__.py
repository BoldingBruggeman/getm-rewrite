import os
import sys
import ctypes

import numpy

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

_pygetm.domain_create.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int]
_pygetm.domain_create.restype = ctypes.c_void_p

_pygetm.domain_get_grid.argtypes = [ctypes.c_void_p, ctypes.c_char]
_pygetm.domain_get_grid.restype = ctypes.c_void_p

_pygetm.grid_get_arrays.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int))]
_pygetm.grid_get_arrays.restype = None

_pygetm.domain_initialize.argtypes = [ctypes.c_void_p]
_pygetm.domain_initialize.restype = None

_pygetm.advection_create.argtypes = []
_pygetm.advection_create.restype = ctypes.c_void_p

CONTIGUOUS = str('CONTIGUOUS')
arrtype2D = numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, ndim=2, flags=CONTIGUOUS)

_pygetm.advection_calculate.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_void_p, arrtype2D, arrtype2D, ctypes.c_double, arrtype2D]
_pygetm.advection_calculate.restype = None

class Grid:
    def __init__(self, domain, grid_type='T'):
        self.p = _pygetm.domain_get_grid(domain.p, grid_type.encode('ascii'))
        pc1 = ctypes.POINTER(ctypes.c_double)()
        pc2 = ctypes.POINTER(ctypes.c_double)()
        pH = ctypes.POINTER(ctypes.c_double)()
        pmask = ctypes.POINTER(ctypes.c_int)()
        _pygetm.grid_get_arrays(self.p, pc1, pc2, pH, pmask)
        halo = domain.halo
        self.c1 = numpy.ctypeslib.as_array(pc1, shape=(domain.shape[2],))[halo:-halo]
        self.c2 = numpy.ctypeslib.as_array(pc2, shape=(domain.shape[1],))[halo:-halo]
        self.H = numpy.ctypeslib.as_array(pH, shape=domain.shape[1:])[halo:-halo,halo:-halo]
        self.mask = numpy.ctypeslib.as_array(pmask, shape=domain.shape[1:])[halo:-halo,halo:-halo]

class Domain:
    def __init__(self, imin, imax, jmin, jmax, kmin, kmax, halo):
        self.imin, self.imax = imin, imax
        self.jmin, self.jmax = jmin, jmax
        self.kmin, self.kmax = kmin, kmax
        self.halo = halo
        self.p = _pygetm.domain_create(imin, imax, jmin, jmax, kmin, kmax, halo)
        self.shape = (kmax - kmin + 1, jmax - jmin + 1 + 2 * halo, imax - imin + 1 + 2 * halo)
        self.T = Grid(self, 'T')

    def initialize(self):
        _pygetm.domain_initialize(self.p)
        self.U = Grid(self, 'U')
        self.V = Grid(self, 'V')

    def array(self, fill=None):
        data = numpy.empty(self.shape[1:])
        if fill is not None:
            data[...] = fill
        return data[self.halo:-self.halo, self.halo:-self.halo], data

class Advection:
    def __init__(self, domain, scheme):
        self.p = _pygetm.advection_create()
        self.pdomain = domain.p
        self.scheme = scheme

    def calculate(self, u, v, dt, var):
        _pygetm.advection_calculate(self.p, self.scheme, self.pdomain, u, v, dt, var)
