from typing import Optional
from . import core
from . import domain

import numpy

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
        loc_h = numpy.empty((self.domain.T.hn.shape[0] + 1,), dtype=self.domain.T.hn.dtype)
        loc_NN = numpy.empty((NN.shape[0],), dtype=NN.dtype)
        loc_SS = numpy.empty((SS.shape[0],), dtype=SS.dtype)
        loc_tke = self.mix.tke
        loc_tkeo = self.mix.tkeo
        loc_eps = self.mix.eps
        loc_L = self.mix.L
        for j in range(self.domain.T.ny):
            for i in range(self.domain.T.nx):
                if self.domain.T.mask[j, i] == 1:   # skip boundary points (mask=2) as SS is not defined there
                    loc_tke[:] = self.tke[:, j, i]
                    #loc_tkeo[:] = self.tkeo[:, j, i]
                    loc_eps[:] = self.eps[:, j, i]
                    loc_L[:] = self.L[:, j, i]
                    loc_h[1:] = self.domain.T.hn[:, j, i]
                    loc_NN[:] = NN[:, j, i]
                    loc_SS[:] = SS[:, j, i]
                    self.mix.nuh[:] = self.nuh[:, j, i]
                    self.mix.num[:] = self.num[:, j, i]
                    self.mix.turbulence(timestep, loc_h, self.domain.T.D[j, i], u_taus[j, i], u_taub[j, i], z0s[j, i], z0b[j, i], loc_NN, loc_SS)
                    self.tke[:, j, i] = loc_tke
                    #self.tkeo[:, j, i] = loc_tkeo
                    self.eps[:, j, i] = loc_eps
                    self.L[:, j, i] = loc_L
                    self.nuh[:, j, i] = self.mix.nuh
                    self.num[:, j, i] = self.mix.num

        # Take viscosity at open boundary from nearest interior point
        # Viscosity (at T points) needs to be valid at open boundary points to so it can be interpolated to inward-adjacent U/V points.
        # However, it cannot be computed as the shear frequency SS is not available at the boundary.
        self.num.update_boundary(ZERO_GRADIENT)