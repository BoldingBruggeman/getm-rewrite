import enum
from typing import Optional, Iterable

from . import _pygetm
from . import domain
from . import core
from . import parallel
from .constants import INTERFACES

class AdvectionScheme(enum.IntEnum):
    HSIMT = 1
    MUSCL = 2
    P2_PDM = 3
    SPLMAX13 = 4
    SUPERBEE = 5
    UPSTREAM = 6 

class Advection(_pygetm.Advection):
    __slots__ = ()

    def __init__(self, grid: domain.Grid, scheme: AdvectionScheme=AdvectionScheme.HSIMT):
        super().__init__(grid, scheme)
    
    def __call__(self, u: core.Array, v: core.Array, timestep: float, var: core.Array, Ah: float=0., skip_initial_halo_exchange: bool=False):
        assert u.grid is self.ugrid
        assert v.grid is self.vgrid
        assert var.grid is self.grid
        self.D[...] = self.grid.D.all_values
        if not skip_initial_halo_exchange:
            var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_2d(v, Ah, 0.5 * timestep, var)
        var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_2d(u, Ah, timestep, var)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_2d(v, Ah, 0.5 * timestep, var)

    def apply_3d(self, u: core.Array, v: core.Array, w: core.Array, timestep: float, var: core.Array, Ah: float=0., new_h: bool=False, skip_initial_halo_exchange: bool=False, w_var: Optional[core.Array]=None):
        if w_var is None:
            w_var = w
        assert u.grid is self.ugrid, 'grid mismatch for u: expected %s, got %s' % (self.ugrid.postfix, u.grid.postfix)
        assert v.grid is self.vgrid, 'grid mismatch for v: expected %s, got %s' % (self.vgrid.postfix, v.grid.postfix)
        assert w.grid is self.grid, 'grid mismatch for w: expected %s, got %s' % (self.grid.postfix, w.grid.postfix)
        assert w.z == INTERFACES, 'grid mismatch for w: expected values at layer interfaces'
        assert w_var.grid is self.grid, 'grid mismatch for w_var: expected %s, got %s' % (self.grid.postfix, w_var.grid.postfix)
        assert w_var.z == INTERFACES, 'grid mismatch for w_var: expected values at layer interfaces'
        assert var.grid is self.grid, 'grid mismatch for advected quantity: expected %s, got %s' % (self.grid.postfix, var.grid.postfix)
        assert var.ndim == 3 and var.z != INTERFACES, 'grid mismatch for advected quantity: expected 3D variable defined at layer centers'
        self.h[...] = (self.grid.hn if new_h else self.grid.ho).all_values
        if not skip_initial_halo_exchange:
            var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_3d(u, Ah, 0.5 * timestep, var)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_3d(v, Ah, 0.5 * timestep, var)
        self.w_3d(w, w_var, timestep, var)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_3d(v, Ah, 0.5 * timestep, var)
        var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_3d(u, Ah, 0.5 * timestep, var)

    def apply_3d_batch(self, u: core.Array, v: core.Array, w: core.Array, timestep: float, vars: Iterable[core.Array], Ah: float=0., new_h: bool=False, skip_initial_halo_exchange: bool=False, w_vars: Optional[Iterable[core.Array]]=None):
        if w_vars is None:
            w_vars = [w] * len(vars)
        assert u.grid is self.ugrid, 'grid mismatch for u: expected %s, got %s' % (self.ugrid.postfix, u.grid.postfix)
        assert v.grid is self.vgrid, 'grid mismatch for v: expected %s, got %s' % (self.vgrid.postfix, v.grid.postfix)
        assert w.grid is self.grid, 'grid mismatch for w: expected %s, got %s' % (self.grid.postfix, w.grid.postfix)
        assert w.z == INTERFACES, 'grid mismatch for w: expected values at layer interfaces'
        for w_var in w_vars:
            assert w_var.grid is self.grid, 'grid mismatch for w_var: expected %s, got %s' % (self.grid.postfix, w_var.grid.postfix)
            assert w_var.z == INTERFACES, 'grid mismatch for w_var: expected values at layer interfaces'
        current_h = (self.grid.hn if new_h else self.grid.ho).all_values.copy()
        for var in vars:
            self.h[...] = current_h
            if not skip_initial_halo_exchange:
                var.update_halos_finish(parallel.Neighbor.LEFT_AND_RIGHT)
            self.u_3d(u, Ah, 0.5 * timestep, var)
            var.update_halos_start(parallel.Neighbor.TOP_AND_BOTTOM)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.TOP_AND_BOTTOM)
            self.v_3d(v, Ah, 0.5 * timestep, var)
        current_h[...] = self.h
        for var, w_var in zip(vars, w_vars):
            self.h[...] = current_h
            self.w_3d(w, w_var, timestep, var)
            var.update_halos_start(parallel.Neighbor.TOP_AND_BOTTOM)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.TOP_AND_BOTTOM)
            self.v_3d(v, Ah, 0.5 * timestep, var)
            var.update_halos_start(parallel.Neighbor.LEFT_AND_RIGHT)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.LEFT_AND_RIGHT)
            self.u_3d(u, Ah, 0.5 * timestep, var)

VerticalDiffusion = _pygetm.VerticalDiffusion
