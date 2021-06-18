import numpy

from . import _pygetm
from . import core
from . import domain
from . import output

Advection = _pygetm.Advection

class Simulation(_pygetm.Simulation):
    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: int=4, apply_bottom_friction: bool=True):
        self.output_manager = output.OutputManager(rank=dom.tiling.rank)
        dom.field_manager = self.output_manager

        assert not dom.initialized
        _pygetm.Simulation.__init__(self, dom, runtype, apply_bottom_friction)

        self.update_depth()

        self.uadv = _pygetm.Advection(dom.U, scheme=advection_scheme)
        self.vadv = _pygetm.Advection(dom.V, scheme=advection_scheme)

        self.uua = dom.UU.array(fill=numpy.nan)
        self.uva = dom.UV.array(fill=numpy.nan)
        self.vua = dom.VU.array(fill=numpy.nan)
        self.vva = dom.VV.array(fill=numpy.nan)

        for name in ('U', 'V', 'fU', 'fV', 'advU', 'advV', 'u1', 'v1'):
            setattr(self, name, self.wrap(core.Array(name=name), name.encode('ascii'), source=1))
        for name in ('dpdx', 'dpdy'):
            setattr(self, name, self.wrap(core.Array(name=name), name.encode('ascii'), source=2))

    def uv_momentum_2d(self, timestep: float, tausx: core.Array, tausy: core.Array, dpdx: core.Array, dpdy: core.Array):
        self.u1.all_values[:, :] = self.U.all_values / self.U.grid.D.all_values
        self.v1.all_values[:, :] = self.V.all_values / self.V.grid.D.all_values

        itimestep = 1. / timestep

        # Update column depth on advection grids
        self.domain.UU.D.all_values[:, :-1] = self.domain.T.D.all_values[:, 1:]
        self.domain.VV.D.all_values[:-1, :] = self.domain.T.D.all_values[1:, :]
        self.domain.UV.D.all_values[:, :] = self.domain.VU.D.all_values[:, :] = self.domain.X.D.all_values[1:, 1:]

        # Advect U using velocities interpolated to its own advection grids
        self.U.interp(self.domain.UU, out=self.uua)
        self.V.interp(self.domain.UV, out=self.uva)
        self.uua.all_values[...] /= self.domain.UU.D.all_values
        self.uva.all_values[...] /= self.domain.UV.D.all_values
        self.uadv.calculate(self.uua, self.uva, timestep, self.u1)
        self.advU.all_values[...] = (self.u1.all_values * self.uadv.D - self.U.all_values) * itimestep

        # Advect V using velocities interpolated to its own advection grids
        self.U.interp(self.domain.VU, out=self.vua)
        self.V.interp(self.domain.VV, out=self.vva)
        self.vua.all_values[...] /= self.domain.VU.D.all_values
        self.vva.all_values[...] /= self.domain.VV.D.all_values
        self.vadv.calculate(self.vua, self.vva, timestep, self.v1)
        self.advV.all_values[...] = (self.v1.all_values * self.vadv.D - self.V.all_values) * itimestep

        _pygetm.Simulation.uv_momentum_2d(self, timestep, tausx, tausy, dpdx, dpdy)

    def update_depth(self):
        # Halo exchange for sea level on T grid
        self.domain.T.z.update_halos()

        # Compute sea level on U, V, X grids
        self.update_sealevel_uvx()

        # Halo exchange for sea level on U, V, X grids
        self.domain.U.z.update_halos()
        self.domain.V.z.update_halos()
        self.domain.X.z.update_halos()

        # Update total water depth D on T, U, V, X grids
        # This also processes the halos; no further halo exchange needed.
        self.domain.update_depths()

    @property
    def Ekin(self, rho0: float=1025.):
        dom = self.domain
        U = self.U.interp(dom.T)
        V = self.V.interp(dom.T)
        vel2_D2 = U**2 + V**2
        return 0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D
