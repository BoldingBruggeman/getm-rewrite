import enum
from typing import Optional, Iterable

from . import _pygetm
from . import domain
from . import core
from . import parallel
from .constants import CENTERS, INTERFACES


class AdvectionScheme(enum.IntEnum):
    HSIMT = 1
    MUSCL = 2
    P2_PDM = 3
    SPLMAX13 = 4
    SUPERBEE = 5
    UPSTREAM = 6


class Advection(_pygetm.Advection):
    __slots__ = ()

    def __init__(
        self, grid: domain.Grid, scheme: AdvectionScheme = AdvectionScheme.HSIMT,
    ):
        super().__init__(grid, scheme)

    def __call__(
        self,
        u: core.Array,
        v: core.Array,
        timestep: float,
        var: core.Array,
        Ah_u: Optional[core.Array] = None,
        Ah_v: Optional[core.Array] = None,
        skip_initial_halo_exchange: bool = False,
    ):
        assert u.grid is self.ugrid and not u.z
        assert v.grid is self.vgrid and not v.z
        assert var.grid is self.grid and not var.z
        self.D[...] = self.grid.D.all_values
        if not skip_initial_halo_exchange:
            var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_2d(v, 0.5 * timestep, var, Ah_v)
        var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_2d(u, timestep, var, Ah_u)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_2d(v, 0.5 * timestep, var, Ah_v)

    def apply_3d(
        self,
        u: core.Array,
        v: core.Array,
        w: core.Array,
        timestep: float,
        var: core.Array,
        Ah_u: Optional[core.Array] = None,
        Ah_v: Optional[core.Array] = None,
        new_h: bool = False,
        skip_initial_halo_exchange: bool = False,
        w_var: Optional[core.Array] = None,
    ):
        if w_var is None:
            w_var = w
        assert u.grid is self.ugrid and u.z == CENTERS
        assert v.grid is self.vgrid and v.z == CENTERS
        assert w.grid is self.grid and w.z == INTERFACES
        assert w_var.grid is self.grid and w_var.z == INTERFACES
        assert var.grid is self.grid and var.z == CENTERS
        self.h[...] = (self.grid.hn if new_h else self.grid.ho).all_values
        if not skip_initial_halo_exchange:
            var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_3d(u, 0.5 * timestep, var, Ah_u)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_3d(v, 0.5 * timestep, var, Ah_v)
        self.w_3d(w, w_var, timestep, var)
        var.update_halos(parallel.Neighbor.TOP_AND_BOTTOM)
        self.v_3d(v, 0.5 * timestep, var, Ah_v)
        var.update_halos(parallel.Neighbor.LEFT_AND_RIGHT)
        self.u_3d(u, 0.5 * timestep, var, Ah_u)

    def apply_3d_batch(
        self,
        u: core.Array,
        v: core.Array,
        w: core.Array,
        timestep: float,
        vars: Iterable[core.Array],
        Ah: Optional[core.Array] = None,
        new_h: bool = False,
        skip_initial_halo_exchange: bool = False,
        w_vars: Optional[Iterable[core.Array]] = None,
    ):
        if w_vars is None:
            w_vars = [w] * len(vars)
        assert u.grid is self.ugrid and u.z == CENTERS
        assert v.grid is self.vgrid and v.z == CENTERS
        assert w.grid is self.grid and w.z == INTERFACES
        for var in vars:
            assert var.grid is self.grid and var.z == CENTERS
        for w_var in w_vars:
            assert w_var.grid is self.grid and w_var.z == INTERFACES
        current_h = (self.grid.hn if new_h else self.grid.ho).all_values.copy()
        for var in vars:
            self.h[...] = current_h
            if not skip_initial_halo_exchange:
                var.update_halos_finish(parallel.Neighbor.LEFT_AND_RIGHT)
            self.u_3d(u, 0.5 * timestep, var, Ah)
            var.update_halos_start(parallel.Neighbor.TOP_AND_BOTTOM)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.TOP_AND_BOTTOM)
            self.v_3d(v, 0.5 * timestep, var, Ah)
        current_h[...] = self.h
        for var, w_var in zip(vars, w_vars):
            self.h[...] = current_h
            self.w_3d(w, w_var, timestep, var)
            var.update_halos_start(parallel.Neighbor.TOP_AND_BOTTOM)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.TOP_AND_BOTTOM)
            self.v_3d(v, 0.5 * timestep, var, Ah)
            var.update_halos_start(parallel.Neighbor.LEFT_AND_RIGHT)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(parallel.Neighbor.LEFT_AND_RIGHT)
            self.u_3d(u, 0.5 * timestep, var, Ah)


VerticalDiffusion = _pygetm.VerticalDiffusion
