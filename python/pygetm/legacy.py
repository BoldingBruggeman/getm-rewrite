from typing import Optional

import numpy
import netCDF4

from . import domain

def domain_from_topo(path: str, nlev: Optional[int]=None, ioffset: int=0, joffset: int=0, nx: Optional[int]=None, ny: Optional[int]=None, z0_const=0.01, **kwargs) -> domain.Domain:
    lon, lat, x, y, z0 = None, None, None, None, None
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        grid_type = int(numpy.reshape(nc['grid_type'], ()))
        if grid_type == 1:
            # cartesian
            pass
        elif grid_type == 2:
            # spherical
            assert nlev is not None
            nclon = nc['lon']
            nclat = nc['lat']
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

            H = domain.read_centers_to_supergrid(nc['bathymetry'], ioffset, joffset, nx, ny)
            z0 = z0_const if 'z0' not in nc.variables else domain.read_centers_to_supergrid(nc['z0'], ioffset, joffset, nx, ny)
            return domain.Domain(nx, ny, nlev, lon=lon, lat=lat, H=numpy.ma.filled(H), z0=numpy.ma.filled(z0), spherical=True, mask=~numpy.ma.getmaskarray(H), **kwargs)
        elif grid_type == 3:
            # planar curvilinear
            pass
        elif grid_type == 4:
            # spherical curvilinear
            pass

def load_bdyinfo(dom: domain.Domain, path: str):
    with open(path) as f:
        def get_line() -> str:
            while True:
                l = f.readline()
                assert l != '', 'End-of-file reached in %s while trying to read next line.' % path
                l = l.rstrip()
                if l and not (l.startswith('!') or l.startswith('#')):
                    return l

        for side in (domain.WEST, domain.NORTH, domain.EAST, domain.SOUTH):
            n = int(get_line())
            for _ in range(n):
                # Note: for Western and Eastern boundaries, l and m are indices in x and y dimensions, respectively, 
                # but that is the other way around (y and , respectively) for Northern and Southern boundaries.
                l, mstart, mstop, type_2d, type_3d = map(int, get_line().split())
                dom.add_open_boundary(side, l, mstart, mstop, type_2d, type_3d)