import numpy

from . import _pygetm, core, domain

Advection = _pygetm.Advection

class Simulation(_pygetm.Simulation):
    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: int=4, apply_bottom_friction: bool=True):
        assert not dom.initialized
        _pygetm.Simulation.__init__(self, dom, runtype, advection_scheme, apply_bottom_friction)

        self.update_depth()

        for name in ('U', 'V', 'fU', 'fV', 'advU', 'advV'):
            setattr(self, name, self.wrap(core.Array(), name.encode('ascii'), source=1))
        for name in ('dpdx', 'dpdy'):
            setattr(self, name, self.wrap(core.Array(), name.encode('ascii'), source=2))

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
