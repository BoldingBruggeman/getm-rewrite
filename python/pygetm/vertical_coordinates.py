from typing import Optional, Union
import logging

import numpy as np
import numpy.typing as npt

from .constants import CENTERS, INTERFACES, TimeVarying, ZERO_GRADIENT, FILL_VALUE
from . import core
from . import _pygetm
from .open_boundaries import ArrayOpenBoundaries


class Base:
    """Base class. It is responsible for updating layer thicknesses
    ``hn`` for all grids provided to initialize"""

    logger: logging.Logger

    def __init__(self, nz: int):
        if nz <= 0:
            raise Exception("Number of layers nz must be a positive number")
        self.nz = nz

    def initialize(
        self, ref_grid: core.Grid, *other_grids: core.Grid, logger: logging.Logger
    ):
        self.logger = logger

    def update(self, timestep: float):
        """Update layer thicknesses hn for all grids"""
        raise NotImplementedError


class PerGrid(Base):
    """Base class for vertical coordinate types that apply the same operation
    to every grid."""

    def initialize(self, *grids: core.Grid, logger: logging.Logger):
        super().initialize(*grids, logger=logger)
        self.grid_info = [self.prepare_update_args(grid) for grid in grids]

    def prepare_update_args(self, grid: core.Grid):
        """Prepare grid-specific information that will be passed as
        arguments to __call__"""
        return (grid.Dclip.all_values, grid.hn.all_values, grid.mask.all_values)

    def update(self, timestep: float):
        """Update all grids"""
        for gi in self.grid_info:
            self(*gi)

    def __call__(
        self,
        D: np.ndarray,
        out: Optional[np.ndarray] = None,
        where: Union[bool, npt.ArrayLike] = True,
    ):
        """Calculate layer thicknesses

        Args:
            D: water depths (m)
            out: array to hold layer thicknesses. It must have shape ``(nz,) + D.shape``
            where: locations where to compute thicknesses (typically: water points).
                It must be broadcastable to the shape of ``D``
        """
        raise NotImplementedError(
            "Classes that inherit from PerGrid must implement __call__"
        )


def calculate_sigma(nz: int, ddl: float = 0.0, ddu: float = 0.0) -> np.ndarray:
    """Return sigma thicknesses (fraction of column depth) of all layers,
    using a formulation that allows zooming towards the surface and bottom.

    Args:
        nz: number of layers
        ddl: zoom factor at bottom (0: no zooming, 2: strong zooming)
        ddu: zoom factor at surface (0: no zooming, 2: strong zooming)
    """
    if ddl <= 0.0 and ddu <= 0.0:
        return np.broadcast_to(1.0 / nz, (nz,))
    ddl, ddu = max(ddl, 0.0), max(ddu, 0.0)

    # This zooming routine is from Antoine Garapon, ICCH, DK
    ga = np.linspace(0.0, 1.0, nz + 1)
    ga[1:-1] = np.tanh((ddl + ddu) * ga[1:-1] - ddl) + np.tanh(ddl)
    ga[1:-1] /= np.tanh(ddl) + np.tanh(ddu)
    dga = ga[1:] - ga[:-1]
    assert (dga > 0.0).all(), f"ga not monotonically increasing: {ga}"
    return dga


class Sigma(PerGrid):
    """Sigma coordinates with optional zooming towards bottom and surface"""

    def __init__(self, nz: int, *, ddl: float = 0.0, ddu: float = 0.0):
        """
        Args:
            nz: number of layers
            ddl: zoom factor at bottom (0: no zooming, 2: strong zooming)
            ddu: zoom factor at surface (0: no zooming, 2: strong zooming)
        """
        super().__init__(nz)
        self.dga = calculate_sigma(nz, ddl, ddu)[:, np.newaxis, np.newaxis]

    def __call__(
        self,
        D: np.ndarray,
        out: Optional[np.ndarray] = None,
        where: Union[bool, npt.ArrayLike] = True,
    ) -> np.ndarray:
        # From sigma thicknesses as fraction [dga] to layer thicknesses in m [hn]
        return np.multiply(self.dga, D, out=out, where=np.asarray(where, dtype=bool))


class GVC(PerGrid):
    """Generalized Vertical Coordinates

    This blends equidistant and surface/bottom-zoomed coordinates as described in
    `Burchard & Petersen (1997)
    <https://doi.org/10.1002/(SICI)1097-0363(19971115)25%3A9%3C1003%3A%3AAID-FLD600%3E3.0.CO%3B2-E>`_.
    It is designed to keep the thickness of either the surface or bottom layer at a constant
    value, except in shallow water where all layers are assigned equal thickness.
    """

    def __init__(
        self,
        nz: int,
        *,
        ddl: float = 0.0,
        ddu: float = 0.0,
        gamma_surf: bool = True,
        Dgamma: float = 0.0,
    ):
        """
        Args:
            nz: number of layers
            ddl: zoom factor at bottom (0: no zooming, 2: strong zooming)
            ddu: zoom factor at surface (0: no zooming, 2: strong zooming)
            gamma_surf: use layers of constant thickness ``Dgamma/nz`` at surface
                (otherwise, at bottom)
            Dgamma: water depth below which to use equal layer thicknesses
        """
        if ddl <= 0.0 and ddu <= 0.0:
            raise Exception("ddl and/or ddu must be a positive number")
        if Dgamma <= 0.0:
            raise Exception("Dgamma must be a positive number")

        super().__init__(nz)

        self.dsigma = 1.0 / nz
        self.dbeta = calculate_sigma(nz, ddl, ddu)
        self.k_ref = -1 if gamma_surf else 0

        if self.dbeta[self.k_ref] >= self.dsigma:
            raise Exception(
                "This GVC parameterization would always result in equidistant layers."
                " If this is desired, use Sigma instead."
            )

        self.Dgamma = Dgamma

        # Calculate valid limit for alpha, and from that, the maximum water depth
        # NB the max alpha would be alpha_lim[alpha_lim > 0.0].min()
        # However, where alpha_lim > 0, we know dsigma - dbeta < 0.
        # Since dsigma - dbeta > -dbeta (while both are < 0), the limit must
        # then be >= 1. This limit is irrelevant as alpha <= 1
        # (unless dbeta[k_ref] > dsigma, but that case was already eliminated above)
        alpha_lim = -self.dbeta / (self.dsigma - self.dbeta)
        alpha_min = alpha_lim[alpha_lim < 0.0].max()
        denom = alpha_min * self.dsigma + (1.0 - alpha_min) * self.dbeta[self.k_ref]
        self.D_max = np.inf if abs(denom) < 1e-15 else (Dgamma * self.dsigma) / denom

    def initialize(self, *grids: core.Grid, logger: logging.Logger):
        super().initialize(*grids, logger=logger)
        self.logger.info(
            f"This GVC parameterization supports water depths up to {self.D_max:.3f} m"
        )

    def __call__(
        self,
        D: np.ndarray,
        out: Optional[np.ndarray] = None,
        where: Optional[np.ndarray] = None,
    ):
        if out is None:
            out = np.empty(self.dbeta.shape + D.shape)
        if where is None:
            where = np.full(D.shape, 1, dtype=np.intc)
        _pygetm.update_gvc(
            self.dsigma, self.dbeta, self.Dgamma, self.k_ref, D, where, out
        )
        return out


class Adaptive(Base):
    """
    Adaptive vertical coordinates based on `Hofmeister et al. (2010)
    <https://doi.org/10.1016/j.ocemod.2009.12.003>`_ and `their GETM
    implementation <https://sourceforge.net/p/getm/code/ci/iow/tree/src/3d/adaptive_coordinates_6.F90>`_
    """

    def __init__(
        self,
        nz: int,
        *,
        ddu: float = 0.0,
        ddl: float = 0.0,
        gamma_surf: bool = True,
        Dgamma: float = 0.0,
        csigma: float = 0.01,
        cgvc: float = 0.0,
        hpow: int = 3,
        chsurf: float = 0.5,
        hsurf: float = 0.5,
        chmidd: float = 0.2,
        hmidd: float = -4.0,
        chbott: float = 0.3,
        hbott: float = -0.25,
        decay: float = 2.0 / 3.0,
        cneigh: float = 0.1,
        rneigh: float = 0.25,
        cNN: float = -1.0,
        drho: float = 0.3,
        cSS: float = -1.0,
        dvel: float = 0.1,
        chmin: float = -0.1,
        hmin: float = 0.3,
        nvfilter: int = 1,
        vfilter: float = 0.2,
        nhfilter: int = 1,
        hfilter: float = 0.1,
        split: int = 1,
        timescale: float = 14400.0,
    ):
        """
        Args:
            nz: number of layers
            ddu: zoom factor at surface (0: no zooming, 2: strong zooming)
            ddl: zoom factor at bottom (0: no zooming, 2: strong zooming)
            gamma_surf: use layers of constant thickness `Dgamma/nz` at surface (otherwise, at bottom)
            Dgamma: water depth below which to use equal layer thicknesses
            hpow: exponent for growth of Dgrid (ramp between 0 and c* tendencies)
            csigma: tendency to uniform sigma
            cgvc: tendency to "standard" gvc (w/ddu,ddl)
            chsurf: tendency to keep surface layer bounded
            hsurf: reference thickness for surface layer
                (absolute thickness in m if >0, relative to average thickness D/nz if <0)
            chmidd: tendency to keep all layers bounded
            hmidd: reference thickness for other layers
                (absolute thickness in m if >0, relative to average thickness D/nz if <0)
            chbott: tendency to keep bottom layer bounded
            hbott: reference thickness for bottom layer
                (absolute thickness in m if >0, relative to average thickness D/nz if <0)
            decay: fraction of surface/bottom tendencies (controlled by chsurf, chbott)
                to preserve for each additional layer away from surface/bottom
                (0: only apply tendency in targeted layer, 1: apply same tendency in all layers)
            cneigh: tendency to keep neighbors of similar size
            rneigh: reference relative growth between neighbors
            cNN: dependence on NN (density zooming)
            drho: reference value for NN density between neighbor cells
            cSS: dependence on SS (shear zooming)
            dvel: reference value for SS absolute shear between neighbor cells
            chmin: internal nug coeff for shallow-water regions
            hmin: minimum depth
            nvfilter: number of vertical filter iterations
            vfilter: strength of vertical filter-of-Dgrid [0:~0.5]
            nhfilter: number of horizontal filter iterations
            hfilter: strength of horizontal filter-of-Dgrid [0:~0.5]
            split: Take this many partial-steps for for vertical filtering == 1 for now
            timescale: time scale of grid adaptation (s)
        """

        if timescale <= 0.0:
            raise Exception("timescale must be a positive value")
        if decay < 0.0 or decay > 1.0:
            raise Exception("decay must be between 0 and 1")
        if chsurf < 0.0:
            raise Exception("chsurf must be non-negative")
        if chbott < 0.0:
            raise Exception("chbott must be non-negative")
        if chmidd < 0.0:
            raise Exception("chmidd must be non-negative")
        if chsurf > 0.0 and hsurf == 0.0:
            raise Exception("hsurf must be non-zero when chsurf is positive")
        if chbott > 0.0 and hbott == 0.0:
            raise Exception("hbott must be non-zero when chbott is positive")
        if chmidd > 0.0 and hmidd == 0.0:
            raise Exception("hmidd must be non-zero when chmidd is positive")
        if nvfilter < 0:
            raise Exception("nvfilter must be non-negative")
        if nhfilter < 0:
            raise Exception("nhfilter must be non-negative")
        if vfilter < 0:
            raise Exception("vfilter must be non-negative")
        if hfilter < 0:
            raise Exception("hfilter must be non-negative")

        if vfilter == 0.0:
            nvfilter = 0
        if hfilter == 0.0:
            nhfilter = 0

        super().__init__(nz)
        self.decay = decay
        self.hpow = hpow
        self.csigma = max(0.0, csigma)
        self.cgvc = max(0.0, cgvc)
        self.chsurf = chsurf
        self.chbott = chbott
        self.chmidd = chmidd
        self.hsurf = hsurf
        self.hbott = hbott
        self.hmidd = hmidd
        self.cneigh = cneigh
        self.rneigh = rneigh
        self.cNN = cNN
        self.drho = drho
        self.cSS = cSS
        self.dvel = dvel
        self.chmin = chmin
        self.hmin = hmin
        self.nvfilter = nvfilter
        self.vfilter = vfilter
        self.nhfilter = nhfilter
        self.hfilter = hfilter
        # self.split = split
        self.split = 1
        self.timescale = timescale

        if (ddl <= 0.0 and ddu <= 0.0) or Dgamma <= 0.0:
            self._gvc = Sigma(nz, ddl=ddl, ddu=ddu)
        else:
            self._gvc = GVC(nz, ddu=ddu, ddl=ddl, gamma_surf=gamma_surf, Dgamma=Dgamma)

    def initialize(
        self,
        tgrid: core.Grid,
        *other_grids: core.Grid,
        logger: logging.Logger,
    ):
        super().initialize(tgrid, *other_grids, logger=logger)
        logger.warning("Support for adaptive vertical coordinates is experimental.")

        if self.csigma > 0.0 and self.cgvc > 0.0:
            logger.warning(f"Relaxing to both sigma and gvc coordinates.")
        elif self.csigma == 0.0 and self.cgvc == 0.0:
            logger.warning(
                f"Not relaxing to any background layer distribution (sigma or gvc)."
            )

        self._gvc.initialize(tgrid, logger=logger)

        self.other_grids = other_grids
        self.dga_t = tgrid.array(
            name="dga",
            z=CENTERS,
            attrs=dict(_time_varying=TimeVarying.MACRO, _mask_output=True),
        )
        self.dga_other = tuple(grid.array(z=CENTERS) for grid in other_grids)

        self.tgrid = tgrid
        self.nug = tgrid.array(
            name="nug",
            units="s-1",
            long_name="vertical grid diffusivity",
            z=CENTERS,
            attrs=dict(_time_varying=TimeVarying.MACRO, _mask_output=True),
            fill_value=FILL_VALUE,
        )

        self.ga = tgrid.array(
            name="ga",
            long_name="gamma coordinate",
            z=INTERFACES,
            attrs=dict(_time_varying=TimeVarying.MACRO, _mask_output=True),
        )

        if hasattr(self.tgrid, "open_boundaries"):
            self.nug.open_boundaries = ArrayOpenBoundaries(self.nug, type=ZERO_GRADIENT)
            self.dga_t.open_boundaries = ArrayOpenBoundaries(
                self.dga_t, type=ZERO_GRADIENT
            )

        # Obtain additional fields used by adaptive coordinates
        # NN and SS should maybe be interpolated to centers
        self.NN = tgrid.fields["NN"]
        self.SS = tgrid.fields["SS"]

    def update(self, timestep: float = 0.0):

        if timestep == 0.0:
            # Simulation is initializing
            self._gvc(self.tgrid.Dclip.all_values, self.tgrid.hn.all_values)
            self.dga_t.all_values = (
                self.tgrid.hn.all_values / self.tgrid.Dclip.all_values
            )

            self.update_other_grids()

            return

        # Reconstruct old sigma positions (from -1 at bottom to 0 at surface)
        # NB the values of ga will usually not be identical to those produced
        # by the previous call to "update" due to river inflow and precipitation.
        # Therefore we compute it anew from old layer thicknesses,
        # which already incorporate these freshwater inputs (unlike, for
        # example, the grid's interface coordinates zf!)
        _pygetm.thickness2interface_depth(self.tgrid.mask, self.tgrid.ho, self.ga)
        self.ga.all_values *= -1.0 / self.ga.all_values[0]

        # Construct the grid diffusivity, which is defined per layer, and equal
        # to the relative rate (tendency, in s-1) at which the layer "gives" its
        # thickness to each neighbor. That is, interior layers (non-surface/bottom)
        # will lose thickness delta_sigma over both interfaces, at a combined rate
        # 2*nug*delta_sigma. Simultaneously they will gain thickness from both of
        # their neighbors at rates nug(k-1)*delta_sigma(k-1) and
        # nug(k+1)*delta_sigma(k+1).
        # While constructing the tendencies, we use dimensionless rates (csigma etc.)
        # These are later divided by the adaptation timescale to obtain
        # the final nug in s-1.

        # Tendency towards pure (non-zoomed) sigma coordinates
        # This is constant for all layers in order to ensure that an equilibrium
        # distribution (equal thickness fluxes in and out) is reached when all
        # layers have equal thickness.
        self.nug.all_values = self.csigma

        # Tendency towards Generalized Vertical Coordinates (incl. zoomed sigma)
        # NB the stationary solution of AVC has thicknesses proportional to
        # the inverse of diffusivity (nug). The division by hn below is
        # therefore necessary; its scaling with surface thicknesses only
        # ensures the tendency equals exactly cgvc at the surface.
        if self.cgvc > 0.0:
            self._gvc(self.tgrid.Dclip.all_values, self.tgrid.hn.all_values)
            self.nug.all_values += np.divide(
                self.cgvc * self.tgrid.hn.all_values[-1], self.tgrid.hn.all_values
            )

        # then add contributions handled by Fortran
        _pygetm.update_adaptive(
            self.nug,
            self.ga,
            self.NN.all_values,
            self.SS.all_values,
            self.decay,
            self.hpow,
            self.chsurf,
            self.hsurf,
            self.chmidd,
            self.hmidd,
            self.chbott,
            self.hbott,
            self.cneigh,
            self.rneigh,
            self.cNN,
            self.drho,
            self.cSS,
            self.dvel,
            self.chmin,
            self.hmin,
        )

        # apply diffusion timescale
        self.nug.all_values *= 1.0 / self.timescale

        # apply vertical filtering from ~/python/src/filters.F90
        if self.nvfilter > 0:
            _pygetm.vertical_filter(self.nvfilter, self.nug, self.vfilter)

        # apply horizontal filtering from ~/python/src/filters.F90
        # requires halo updates
        for _ in range(self.nhfilter):
            if self.nug.open_boundaries is not None:
                self.nug.open_boundaries.update()
            self.nug.update_halos()
            _pygetm.horizontal_filter(self.nug, self.hfilter)

        # now the grid diffusion field is ready to be applied
        _pygetm.tridiagonal(self.nug, self.ga, timestep)

        # assure consistent ga on boundaries
        self.ga.open_boundaries.update()

        # Calculate sigma thicknesses from interface positions
        np.subtract(
            self.ga.all_values[1:], self.ga.all_values[:-1], out=self.dga_t.all_values
        )

        # assure consistent dga on boundaries and in halo zones
        if self.dga_t.open_boundaries is not None:
            self.dga_t.open_boundaries.update()
        self.dga_t.update_halos()

        self.update_other_grids()

    def update_other_grids(self):
        """Update layer thicknesses hn for all other grids"""
        # Interpolate dga from T grid to other grids
        for dga in self.dga_other:
            self.dga_t.interp(dga)
            dga.mirror()
            dga.update_halos()

        # From dga to layer thicknesses in m [hn]
        all_grids = (self.tgrid,) + self.other_grids
        all_dga = (self.dga_t,) + self.dga_other
        for grid, dga in zip(all_grids, all_dga):
            np.multiply(
                dga, grid.Dclip.all_values, where=grid._water, out=grid.hn.all_values
            )
