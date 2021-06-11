import numpy

from . import _pygetm, core, domain

Advection = _pygetm.Advection

class Simulation(_pygetm.Simulation, core.FortranObject):
    def __init__(self, dom: domain.Domain, runtype: int, advection_scheme: int=4, apply_bottom_friction: bool=True):
        assert not dom.initialized
        self.domain = dom
        _pygetm.Simulation.__init__(self, dom, runtype, advection_scheme, apply_bottom_friction)

        self.dist_zT = None
        if self.domain.tiling:
            self.dist_zT = self.domain.distribute(self.domain.T.z_)
            self.dist_zU = self.domain.distribute(self.domain.U.z_)
            self.dist_zV = self.domain.distribute(self.domain.V.z_)
            self.dist_zX = self.domain.distribute(self.domain.X.z_)
            self.update_depth()

        self.get_arrays(('U', 'V', 'fU', 'fV', 'advU', 'advV'), source=0, halo=self.domain.halo)
        self.get_arrays(('dpdx', 'dpdy'), source=1, halo=self.domain.halo)

    def update_depth(self):
        # Halo exchange for sea level on T grid
        if self.dist_zT:
            self.dist_zT.update_halos()

        # Compute sea level on U, V, X grids
        self.update_sealevel_uvx()

        if self.dist_zT:
            # Halo exchange for sea level on U, V, X grids
            self.dist_zU.update_halos()
            self.dist_zV.update_halos()
            self.dist_zX.update_halos()

        # Update total water depth D on T, U, V, X grids
        # This also processes the halos; no further halo exchange needed.
        self.domain.update_depths()

    @property
    def Ekin(self, rho0: float=1025.):
        dom = self.domain
        U = dom.T.map(self.U_, dom.U)[2:-2, 2:-2]
        V = dom.T.map(self.V_, dom.V)[2:-2, 2:-2]
        vel2_D2 = U**2 + V**2
        return numpy.ma.array(0.5 * rho0 * dom.T.area * vel2_D2 / dom.T.D, mask=dom.T.mask==0)
