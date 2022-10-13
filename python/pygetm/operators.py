import enum
from typing import Optional, Iterable
import functools

from . import _pygetm
from . import domain
from . import core
from . import parallel
from .constants import CENTERS, INTERFACES


class AdvectionScheme(enum.IntEnum):
    HSIMT = 1  #: `Wu & Zhu (2010) <https://doi.org/10.1016/j.ocemod.2009.12.001>`_
    MUSCL = 2
    P2_PDM = 3
    SPLMAX13 = 4
    SUPERBEE = 5
    UPSTREAM = 6
    DEFAULT = SUPERBEE


class AdvectionSplit(enum.Enum):
    FULL = 1  #: full splitting (first order in time): u-v in 2D, u-v-w in 3D
    HALF = 2  #: Strang splitting (second order in time): u/2-v-u/2 in 2D, u/2-v/2-w-v/2-u/2 in 3D
    HALF_ALWAYS = 3  #: Strang splitting (second order in time): u/2-v/2-v/2-u/2 in 2D, u/2-v/2-w-v/2-u/2 in 3D


class Advection(_pygetm.Advection):
    __slots__ = ("_ufirst", "halo1", "halo2", "split_2d")

    def __init__(
        self,
        grid: domain.Grid,
        scheme: AdvectionScheme = AdvectionScheme.DEFAULT,
        split_2d: AdvectionSplit = AdvectionSplit.HALF,
    ):
        super().__init__(grid, scheme)
        self.ufirst = False
        self.split_2d = split_2d

    @property
    def ufirst(self) -> bool:
        return self._ufirst

    @ufirst.setter
    def ufirst(self, value: bool):
        self._ufirst = value
        self.halo1 = parallel.Neighbor.LEFT_AND_RIGHT
        self.halo2 = parallel.Neighbor.TOP_AND_BOTTOM
        if not self._ufirst:
            self.halo1, self.halo2 = self.halo2, self.halo1

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
        assert Ah_u is None or (Ah_u.grid is self.ugrid and not Ah_u.z)
        assert Ah_v is None or (Ah_v.grid is self.vgrid and not Ah_v.z)

        adv1 = functools.partial(self.u_2d, u, Ah=Ah_u)
        adv2 = functools.partial(self.v_2d, v, Ah=Ah_v)
        if not self._ufirst:
            adv1, adv2 = adv2, adv1

        self.D[...] = self.grid.D.all_values
        if not skip_initial_halo_exchange:
            var.update_halos(self.halo1)
        if self.split_2d == AdvectionSplit.FULL:
            adv1(timestep, var)
            var.update_halos(self.halo2)
            adv2(timestep, var)
        elif self.split_2d == AdvectionSplit.HALF:
            adv1(0.5 * timestep, var)
            var.update_halos(self.halo2)
            adv2(timestep, var)
            var.update_halos(self.halo1)
            adv1(0.5 * timestep, var)
        else:
            adv1(0.5 * timestep, var)
            var.update_halos(self.halo2)
            adv2(0.5 * timestep, var)
            var.update_halos(self.halo2)
            adv2(0.5 * timestep, var)
            var.update_halos(self.halo1)
            adv1(0.5 * timestep, var)

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
        assert Ah_u is None or (Ah_u.grid is self.ugrid and not Ah_u.z)
        assert Ah_v is None or (Ah_v.grid is self.vgrid and not Ah_v.z)

        adv1 = functools.partial(self.u_3d, u, 0.5 * timestep, var, Ah_u)
        adv2 = functools.partial(self.v_3d, v, 0.5 * timestep, var, Ah_v)
        if not self._ufirst:
            adv1, adv2 = adv2, adv1

        self.h[...] = (self.grid.hn if new_h else self.grid.ho).all_values
        if not skip_initial_halo_exchange:
            var.update_halos(self.halo1)
        adv1()
        var.update_halos(self.halo2)
        adv2()
        self.w_3d(w, w_var, timestep, var)
        var.update_halos(self.halo2)
        adv2()
        var.update_halos(self.halo1)
        adv1()

    def apply_3d_batch(
        self,
        u: core.Array,
        v: core.Array,
        w: core.Array,
        timestep: float,
        vars: Iterable[core.Array],
        Ah_u: Optional[core.Array] = None,
        Ah_v: Optional[core.Array] = None,
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

        adv1 = functools.partial(self.u_3d, u, 0.5 * timestep, Ah=Ah_u)
        adv2 = functools.partial(self.v_3d, v, 0.5 * timestep, Ah=Ah_v)
        if not self._ufirst:
            adv1, adv2 = adv2, adv1

        current_h = (self.grid.hn if new_h else self.grid.ho).all_values.copy()
        for var in vars:
            self.h[...] = current_h
            if not skip_initial_halo_exchange:
                var.update_halos_finish(self.halo1)
            adv1(var)
            var.update_halos_start(self.halo2)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(self.halo2)
            adv2(var)
        current_h[...] = self.h
        for var, w_var in zip(vars, w_vars):
            self.h[...] = current_h
            self.w_3d(w, w_var, timestep, var)
            var.update_halos_start(self.halo2)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(self.halo2)
            adv2(var)
            var.update_halos_start(self.halo1)
        current_h[...] = self.h
        for var in vars:
            self.h[...] = current_h
            var.update_halos_finish(self.halo1)
            adv1(var)


VerticalDiffusion = _pygetm.VerticalDiffusion
