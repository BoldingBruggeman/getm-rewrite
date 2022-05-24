from typing import Optional

import numpy
import netCDF4

import pygetm.domain

def domain_from_topo(path: str, nlev: Optional[int]=None, ioffset: int=0, joffset: int=0, nx: Optional[int]=None, ny: Optional[int]=None, **kwargs) -> pygetm.domain.Domain:
    """Create a domain object from a topo.nc file used by legacy GETM.
    
    Args:
        path: NetCDF file with legacy domain topography
        nlev: number of vertical layers
        ioffset: starting point in x direction (default: 0 = start of the domain in the topo file)
        joffset: starting point in y direction (default: 0 = start of the domain in the topo file)
        nx: number of points to read in x direction (default: read till end of the domain in the topo file)
        ny: number of points to read in y direction (default: read till end of the domain in the topo file)
        **kwargs: keyword arguments that are ultimately passed to :class:`pygetm.domain.Domain`

    Bottom roughness ``z0`` can be prescribed with keyword argument ``z0``, which will be passed to :class:`pygetm.domain.Domain`.
    This argument `must` be provided if the topo file does not contain bottom roughness ``z0``.
    """
    lon, lat, x, y, z0 = None, None, None, None, None
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        grid_type = int(numpy.reshape(nc['grid_type'], ()))
        if grid_type == 1:
            # Cartesian
            raise NotImplementedError('No support yet for Cartesian coordinates')
        elif grid_type == 2:
            # spherical
            assert nlev is not None
            latname, lonname = nc['bathymetry'].dimensions
            nclon = nc[lonname]
            nclat = nc[latname]
            if nx is None:
                nx = nclon.size - ioffset
            if ny is None:
                ny = nclat.size - joffset
            dlon = (nclon[-1] - nclon[0]) / (nclon.size - 1)
            dlat = (nclat[-1] - nclat[0]) / (nclat.size - 1)

            # Define lon, lat on supergrid
            lon = nclon[0] + (ioffset - 0.5) * dlon + numpy.arange(2 * nx + 1) * (0.5 * dlon)
            lat = nclat[0] + (joffset - 0.5) * dlat + numpy.arange(2 * ny + 1) * (0.5 * dlat)
            lat = lat[:, numpy.newaxis]

            H = pygetm.domain.read_centers_to_supergrid(nc['bathymetry'], ioffset, joffset, nx, ny)
            if 'z0' not in kwargs:
                if 'z0' not in nc.variables:
                    raise Exception('Bottom roughness z0 is not present in %s; you need to provide it as keyword argument to domain_from_topo instead.')
                kwargs['z0'] = numpy.ma.filled(pygetm.domain.read_centers_to_supergrid(nc['z0'], ioffset, joffset, nx, ny))
            domain = pygetm.domain.Domain.create(nx, ny, nlev, lon=lon, lat=lat, H=numpy.ma.filled(H), spherical=True, mask=numpy.where(numpy.ma.getmaskarray(H), 0, 1), **kwargs)
        elif grid_type == 3:
            # planar curvilinear
            raise NotImplementedError('No support yet for planar curvilinear coordinates')
        elif grid_type == 4:
            # spherical curvilinear
            raise NotImplementedError('No support yet for spherical curvilinear coordinates')
        else:
            raise NotImplementedError('Unknown grid_type %i found' % grid_type)
    return domain

class DatFile:
    """Support for reading GETM dat files with comments indicated by ! or #.
    Whitespace-only lines are skipped."""
    def __init__(self, path: str):
        self.path = path
        self.f = open(path)

    def get_line(self) -> str:
        """Return next non-empty line"""
        l = None
        while not l:
            l = self.f.readline()
            assert l != '', 'End-of-file reached in %s while trying to read next line.' % self.path
            l = l.split('#', 1)[0].split('!', 1)[0].strip()
        return l

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.f.close()

def load_bdyinfo(dom: pygetm.domain.Domain, path: str, type_2d: Optional[int]=None, type_3d: Optional[int]=None):
    """Add open boundaries from bdyinfo.dat to domain.
    
    Args:
        path: data file with open boundary information
        type_2d: type of 2D open boundary condition to use. If provided, this overrides the type configured in the file.
        type_3d: type of 3D open boundary condition to use. If provided, this overrides the type configured in the file.
    """
    with DatFile(path) as f:
        for side in (pygetm.domain.WEST, pygetm.domain.NORTH, pygetm.domain.EAST, pygetm.domain.SOUTH):
            n = int(f.get_line())
            for _ in range(n):
                # Note: for Western and Eastern boundaries, l and m are indices in x and y dimensions, respectively, 
                # but that is the other way around (y and x, respectively) for Northern and Southern boundaries.
                # Note that indices are 1-based as in Fortran. We convert to the Python convention: 0-based indices,
                # with the upper bound being the first index that is EXcluded.
                l, mstart, mstop, type_2d_, type_3d_ = map(int, f.get_line().split())
                dom.open_boundaries.add_by_index(side, l - 1, mstart - 1, mstop, type_2d_ if type_2d is None else type_2d, type_3d_ if type_3d is None else type_3d)

def load_riverinfo(dom: pygetm.domain.Domain, path: str):
    """Add rivers from riverinfo.dat to domain
    
    Args:
        path: data file with river information
    """
    with DatFile(path) as f:
        n = int(f.get_line())
        for _ in range(n):
            items = f.get_line().split()
            assert len(items) in (3, 5)
            i, j, name = int(items[0]), int(items[1]), items[2]
            zl, zu = None, None
            if len(items) == 5:
                zl, zu = float(items[3]), float(items[4])
            dom.rivers.add_by_index(name, i - 1, j - 1, zl=zl, zu=zu)   # Note: we convert from 1-based indices to 0-based indices!
