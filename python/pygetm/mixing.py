from typing import Optional

import numpy

from . import core
from . import domain
from . import _pygotm
from .constants import *

class Turbulence:
    def __init__(self, domain: domain.Domain):
        self.nuh = domain.T.array(z=INTERFACES, name='nuh', units='m2 s-1', long_name='turbulent diffusivity of heat', fill_value=FILL_VALUE)
        self.num = domain.T.array(z=INTERFACES, name='num', units='m2 s-1', long_name='turbulent diffusivity of momentum', fill_value=FILL_VALUE)
        self.nuh.fill(0.)
        self.num.fill(0.)

class GOTM(Turbulence):
    def __init__(self, domain: domain.Domain, nml_path: Optional[str]=None):
        super().__init__(domain)
        self.mix = _pygotm.Mixing(domain.T.nz, b'' if nml_path is None else nml_path.encode('ascii'))
        self.domain = domain
        self.nuh.fill(self.mix.nuh[:, numpy.newaxis, numpy.newaxis])
        self.num.fill(self.mix.num[:, numpy.newaxis, numpy.newaxis])
        self.tke = domain.T.array(fill=self.mix.tke[:, numpy.newaxis, numpy.newaxis], z=INTERFACES, name='tke', units='m2 s-2', long_name='turbulent kinetic energy')
        self.tkeo = domain.T.array(fill=self.mix.tkeo[:, numpy.newaxis, numpy.newaxis], z=INTERFACES, name='tkeo', units='m2 s-2', long_name='turbulent kinetic energy at previous timestep')
        self.eps = domain.T.array(fill=self.mix.eps[:, numpy.newaxis, numpy.newaxis], z=INTERFACES, name='eps', units='m2 s-3', long_name='energy dissipation rate')
        self.L = domain.T.array(fill=self.mix.L[:, numpy.newaxis, numpy.newaxis], z=INTERFACES, name='L', units='m', long_name='turbulence length scale')

    def __call__(self, timestep: float, u_taus: core.Array, u_taub: core.Array, z0s: core.Array, z0b: core.Array, NN: core.Array, SS: core.Array):
        assert u_taus.grid is self.domain.T and u_taus.ndim == 2
        assert u_taub.grid is self.domain.T and u_taub.ndim == 2
        assert z0s.grid is self.domain.T and z0s.ndim == 2
        assert z0b.grid is self.domain.T and z0b.ndim == 2
        assert NN.grid is self.domain.T and NN.z == INTERFACES
        assert SS.grid is self.domain.T and SS.z == INTERFACES

        nz, ny, nx = self.domain.T.hn.all_values.shape
        self.mix.turbulence_3d(nx, ny, nz, 2, nx - 2, 2, ny - 2, timestep, self.domain.T.mask.all_values,
            self.domain.T.hn.all_values, self.domain.T.D.all_values, u_taus.all_values, u_taub.all_values,
            z0s.all_values, z0b.all_values, NN.all_values, SS.all_values,
            self.tke.all_values, self.tkeo.all_values, self.eps.all_values, self.L.all_values, self.num.all_values, self.nuh.all_values)

        # Take viscosity at open boundary from nearest interior point
        # Viscosity (at T points) needs to be valid at open boundary points to so it can be interpolated to inward-adjacent U/V points.
        # However, it cannot be computed as the shear frequency SS is not available at the boundary.
        self.num.update_boundary(ZERO_GRADIENT)