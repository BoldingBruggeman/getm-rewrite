from typing import Optional
from . import core
from . import domain

import numpy

import pygsw
from . import _pygotm

def get_buoyancy_frequency(SA: core.Array, ct: core.Array, p: core.Array=None, out: core.Array=None):
    assert SA.grid is ct.grid
    assert SA.ndim == 3
    if out is None:
        out = SA.grid.wgrid.array(is_3d=True)
    assert out.grid is SA.grid.wgrid
    if p is None:
        p = -SA.grid.zc
    assert p.grid is SA.grid
    pygsw.nsquared(SA.grid.hn.all_values, SA.all_values, ct.all_values, p.all_values, SA.grid.lat.all_values, out.all_values[1:-1, :, :])
    return out

class Turbulence:
    def __init__(self, domain: domain.Domain):
        self.mix = _pygotm.Mixing(domain.nz)
        self.domain = domain
        self.tke = domain.W.array(fill=0., is_3d=True, name='tke', units='m2 s-2', long_name='turbulent kinetic energy')
        self.eps = domain.W.array(fill=0., is_3d=True, name='eps', units='m2 s-3', long_name='energy dissipation rate')
        self.L = domain.W.array(fill=0., is_3d=True, name='L', units='m', long_name='turbulence length scale')
        self.nuh = domain.W.array(fill=0., is_3d=True, name='nuh', units='m2 s-1', long_name='turbulent diffusivity of heat')
        self.num = domain.W.array(fill=0., is_3d=True, name='num', units='m2 s-1', long_name='turbulent diffusivity of momentum')

    def __call__(self, timestep: float, u_taus: core.Array, u_taub: core.Array, z0s: core.Array, z0b: core.Array, NN: core.Array, SS: core.Array):
        assert u_taus.grid is self.domain.T and u_taus.ndim == 2
        assert u_taub.grid is self.domain.T and u_taub.ndim == 2
        assert z0s.grid is self.domain.T and z0s.ndim == 2
        assert z0b.grid is self.domain.T and z0b.ndim == 2
        assert NN.grid is self.domain.W and NN.ndim == 3
        assert SS.grid is self.domain.W and SS.ndim == 3
        loc_h = numpy.empty((self.domain.T.hn.shape[0] + 1,), dtype=self.domain.T.hn.dtype)
        loc_NN = numpy.empty((NN.shape[0],), dtype=NN.dtype)
        loc_SS = numpy.empty((SS.shape[0],), dtype=SS.dtype)
        for j in range(self.domain.ny):
            for i in range(self.domain.nx):
                if self.domain.T.mask[j, i] == 1:
                    self.mix.tke[:] = self.tke[:, j, i]
                    self.mix.eps[:] = self.eps[:, j, i]
                    self.mix.L[:] = self.L[:, j, i]
                    loc_h[1:] = self.domain.T.hn[:, j, i]
                    loc_NN[:] = NN[:, j, i]
                    loc_SS[:] = SS[:, j, i]
                    self.mix.turbulence(timestep, loc_h, self.domain.T.D[j, i], u_taus[j, i], u_taub[j, i], z0s[j, i], z0b[j, i], loc_NN, loc_SS)
                    self.tke[:, j, i] = self.mix.tke
                    self.eps[:, j, i] = self.mix.eps
                    self.L[:, j, i] = self.mix.L
                    self.nuh[:, j, i] = self.mix.nuh
                    self.num[:, j, i] = self.mix.num
