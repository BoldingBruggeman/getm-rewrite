import enum
import operator
import logging
from typing import Optional

import numpy as np

from . import core
from . import parallel
from . import operators
import pygetm.domain
import pygetm._pygetm
from .constants import FILL_VALUE, BAROTROPIC_2D, CENTERS


class CoriolisScheme(enum.IntEnum):
    OFF = 0
    DEFAULT = 1
    ESPELID = 2  #: `Espelid et al. (2000) <https://doi.org/10.1002/1097-0207(20001230)49:12<1521::AID-NME9>3.0.CO;2-F>`_


MASK_ZERO_2D = ("U", "V", "u1", "v1")
MASK_ZERO_3D = ("pk", "qk", "uk", "vk", "Ui", "Vi")


class Momentum(pygetm._pygetm.Momentum):
    _arrays = (
        "U",
        "V",
        "fU",
        "fV",
        "advU",
        "advV",
        "diffU",
        "diffV",
        "dampU",
        "dampV",
        "u1",
        "v1",
        "uk",
        "vk",
        "ru",
        "rru",
        "rv",
        "rrv",
        "pk",
        "qk",
        "ww",
        "advpk",
        "advqk",
        "diffpk",
        "diffqk",
        "Ui",
        "Vi",
        "SS",
        "fpk",
        "fqk",
        "ustar2_s",
        "ustar2_b",
        "SxB",
        "SyB",
        "SxA",
        "SyA",
        "SxD",
        "SyD",
        "SxF",
        "SyF",
    )
    _all_fortran_arrays = tuple(["_%s" % name for name in _arrays]) + (
        "uadv",
        "vadv",
        "uua",
        "uva",
        "vua",
        "vva",
        "uua3d",
        "uva3d",
        "vua3d",
        "vva3d",
    )
    __slots__ = _all_fortran_arrays + (
        "runtype",
        "logger",
        "Ui_tmp",
        "Vi_tmp",
        "_ufirst",
        "_u3dfirst",
        "diffuse_momentum",
        "apply_bottom_friction",
        "_Am_const",
        "_An_const",
        "cnpar",
        "advection_scheme",
        "advection_split_2d",
        "coriolis_scheme",
        "An",
        "An_uu",
        "An_uv",
        "An_vu",
        "An_vv",
        "u_V",
        "v_U",
        "hU_bot",
        "hV_bot",
        "u_bot",
        "v_bot",
        "udev",
        "vdev",
    )

    _array_args = {
        "U": dict(
            units="m2 s-1",
            long_name="depth-integrated transport in x-direction",
            attrs=dict(_part_of_state=True, _mask_output=True),
        ),
        "V": dict(
            units="m2 s-1",
            long_name="depth-integrated transport in y-direction",
            attrs=dict(_part_of_state=True, _mask_output=True),
        ),
        "Ui": dict(units="m2 s-1", attrs=dict(_part_of_state=True, _mask_output=True),),
        "Vi": dict(units="m2 s-1", attrs=dict(_part_of_state=True, _mask_output=True),),
        "u1": dict(
            units="m s-1",
            long_name="depth-averaged velocity in x-direction",
            attrs=dict(_mask_output=True),
        ),
        "v1": dict(
            units="m s-1",
            long_name="depth-averaged velocity in y-direction",
            attrs=dict(_mask_output=True),
        ),
        "pk": dict(
            units="m2 s-1",
            long_name="layer-integrated transport in x-direction",
            attrs=dict(_part_of_state=True, _mask_output=True),
        ),
        "qk": dict(
            units="m2 s-1",
            long_name="layer-integrated transport in y-direction",
            attrs=dict(_part_of_state=True, _mask_output=True),
        ),
        "uk": dict(
            units="m s-1",
            long_name="velocity in x-direction",
            attrs=dict(_mask_output=True, standard_name="sea_water_x_velocity"),
        ),
        "vk": dict(
            units="m s-1",
            long_name="velocity in y-direction",
            attrs=dict(_mask_output=True, standard_name="sea_water_y_velocity"),
        ),
        "ww": dict(
            units="m s-1",
            long_name="vertical velocity",
            attrs=dict(standard_name="upward_sea_water_velocity"),
        ),
        "SS": dict(units="s-2", long_name="shear frequency squared"),
        "SxA": dict(
            units="m-2 s-2",
            long_name="slow advection in x-direction",
            attrs=dict(_mask_output=True),
        ),
        "SyA": dict(
            units="m-2 s-2",
            long_name="slow advection in y-direction",
            attrs=dict(_mask_output=True),
        ),
        "SxB": dict(
            units="m-2 s-2",
            long_name="depth-integrated internal pressure in x-direction",
            attrs=dict(_mask_output=True),
        ),
        "SyB": dict(
            units="m-2 s-2",
            long_name="depth-integrated internal pressure in y-direction",
            attrs=dict(_mask_output=True),
        ),
        "SxD": dict(
            units="m-2 s-2",
            long_name="slow diffusion in x-direction",
            attrs=dict(_mask_output=True),
        ),
        "SyD": dict(
            units="m-2 s-2",
            long_name="slow diffusion in y-direction",
            attrs=dict(_mask_output=True),
        ),
        "SxF": dict(
            units="m-2 s-2",
            long_name="slow bottom friction in x-direction",
            attrs=dict(_mask_output=True),
        ),
        "SyF": dict(
            units="m-2 s-2",
            long_name="slow bottom friction in y-direction",
            attrs=dict(_mask_output=True),
        ),
        "fU": dict(attrs=dict(_mask_output=True),),
        "fV": dict(attrs=dict(_mask_output=True),),
        "fpk": dict(attrs=dict(_mask_output=True),),
        "fqk": dict(attrs=dict(_mask_output=True),),
        "advU": dict(attrs=dict(_mask_output=True),),
        "advV": dict(attrs=dict(_mask_output=True),),
        "advpk": dict(attrs=dict(_mask_output=True),),
        "advqk": dict(attrs=dict(_mask_output=True),),
    }

    def __init__(
        self,
        Am: float = 0.0,
        An: float = 0.0,
        cnpar: float = 1.0,
        advection_scheme: Optional[operators.AdvectionScheme] = None,
        advection_split_2d: operators.AdvectionSplit = operators.AdvectionSplit.FULL,
        coriolis_scheme: CoriolisScheme = CoriolisScheme.DEFAULT,
    ):
        """Create momentum handler

        Args:
            Am: horizontal viscosity (m2 s-1)
                If provided, this must be a constant.
            An: horizontal diffusivity (m2 s-1)
                If provided, this must be a constant.
                It can subsequently be changed through attribute :attr:`An`, which also
                allows spatially varying diffusivities to be set.
            cnpar: parameter for the Crank–Nicolson vertical diffusion solver
            advection_scheme: advection scheme
            advection_split_2d: directional splitting for advection solver
            coriolis_scheme: interpolation method to use to recontruct velocities for
                Coriolis terms
        """
        self._Am_const = Am
        self._An_const = An
        self.cnpar = cnpar
        self.advection_scheme = advection_scheme
        self.advection_split_2d = advection_split_2d
        self.coriolis_scheme = coriolis_scheme

    def initialize(
        self,
        logger: logging.Logger,
        domain: pygetm.domain.Domain,
        runtype: int,
        default_advection_scheme: operators.AdvectionScheme,
    ):
        super().__init__(
            domain, runtype, self._Am_const, self.cnpar, self.coriolis_scheme
        )

        for name in self._arrays:
            kwargs = dict(fill_value=FILL_VALUE)
            kwargs.update(self._array_args.get(name, {}))
            setattr(
                self,
                "_%s" % name,
                self.wrap(core.Array(name=name, **kwargs), name.encode("ascii"),),
            )

        self.logger = logger
        self.runtype = runtype

        # Disable bottom friction if physical bottom roughness is 0 everywhere
        z0b = domain.z0b_min[domain.mask > 0]
        self.apply_bottom_friction = not (z0b == 0.0).any()
        if not self.apply_bottom_friction:
            self.logger.warning(
                "Disabling bottom friction because bottom roughness is 0 in"
                " %i of %i unmasked supergrid cells." % ((z0b == 0.0).sum(), z0b.size)
            )

        self.diffuse_momentum = self._Am_const > 0.0
        if not self.diffuse_momentum:
            self.logger.info(
                "Disabling horizontal diffusion because horizontal viscosity Am is 0"
            )

        ZERO_EVERYWHERE = MASK_ZERO_2D
        ZERO_UNMASKED = (
            "SxA",
            "SyA",
            "SxB",
            "SyB",
            "SxD",
            "SyD",
            "SxF",
            "SyF",
            "dampU",
            "dampV",
            "advU",
            "advV",
            "diffU",
            "diffV",
            "ru",
            "rv",
            "rru",
            "rrv",
            "fU",
            "fV",
        )
        if runtype > BAROTROPIC_2D:
            ZERO_EVERYWHERE = ZERO_EVERYWHERE + MASK_ZERO_3D
            ZERO_UNMASKED = ZERO_UNMASKED + (
                "SS",
                "advpk",
                "advqk",
                "diffpk",
                "diffqk",
                "fpk",
                "fqk",
                "ww",
            )
        for v in ZERO_EVERYWHERE:
            array = getattr(self, v)
            array.all_values[...] = 0.0
        for v in ZERO_UNMASKED:
            getattr(self, v).fill(0.0)

        # Zero elements of Coriolis terms that wil never be set, but still
        # multiplied by f. This prevents overflow warnings
        # We could zero the whole arrays, but by doing this selectively, we
        # can verify in tests that values on other masked points are not used
        self.fU.all_values[:, 0] = 0.0
        self.fU.all_values[-1, :] = 0.0
        self.fV.all_values[:, -1] = 0.0
        self.fV.all_values[0, :] = 0.0
        if runtype > BAROTROPIC_2D:
            self.fpk.all_values[:, :, 0] = 0.0
            self.fpk.all_values[:, -1, :] = 0.0
            self.fqk.all_values[:, :, -1] = 0.0
            self.fqk.all_values[:, 0, :] = 0.0

        self.Ui_tmp = np.zeros_like(self.Ui.all_values)
        self.Vi_tmp = np.zeros_like(self.Vi.all_values)

        if self.advection_scheme is None:
            self.advection_scheme = default_advection_scheme
        self.logger.info("Advection scheme: %s" % self.advection_scheme.name)
        self.uadv = operators.Advection(
            domain.U, scheme=self.advection_scheme, split_2d=self.advection_split_2d
        )
        self.vadv = operators.Advection(
            domain.V, scheme=self.advection_scheme, split_2d=self.advection_split_2d
        )

        self.logger.info("Crank-Nicolson parameter: %s" % self.cnpar)
        self.logger.info("Coriolis interpolation: %s" % self.coriolis_scheme.name)

        self.uua = domain.UU.array(fill=np.nan)
        self.uva = domain.UV.array(fill=np.nan)
        self.vua = domain.VU.array(fill=np.nan)
        self.vva = domain.VV.array(fill=np.nan)

        self.u_V = domain.V.array()
        self.v_U = domain.U.array()

        if self.runtype > BAROTROPIC_2D:
            self.uua3d = domain.UU.array(fill=np.nan, z=CENTERS)
            self.uva3d = domain.UV.array(fill=np.nan, z=CENTERS)
            self.vua3d = domain.VU.array(fill=np.nan, z=CENTERS)
            self.vva3d = domain.VV.array(fill=np.nan, z=CENTERS)

            self.hU_bot = domain.U.hn.isel(z=0)
            self.hV_bot = domain.V.hn.isel(z=0)
            self.u_bot = self.uk.isel(z=0)
            self.v_bot = self.vk.isel(z=0)

        self.An = domain.T.array(
            name="An",
            units="m2 s-1",
            long_name="horizontal diffusivity of momentum",
            fill_value=FILL_VALUE,
            attrs=dict(_require_halos=True, _time_varying=False),
        )
        self.An.fill(self._An_const)
        self.An_uu = self.An_uv = self.An_vu = self.An_vv = None

        #: Whether to start the depth-integrated (2D) momentum update with u
        # (as opposed to v)
        self._ufirst = False

        #: Whether to start the depth-explicit (3D) momentum update with u
        # (as opposed to v)
        self._u3dfirst = False

        self.udev = domain.U.array(
            name="udev",
            units="m s-1",
            attrs=dict(_mask_output=True),
            fill_value=FILL_VALUE,
        )
        self.vdev = domain.V.array(
            name="vdev",
            units="m s-1",
            attrs=dict(_mask_output=True),
            fill_value=FILL_VALUE,
        )

    def start(self):
        # Ensure transports and velocities are 0 in masked points
        # NB velocities will be computed from transports, but only in unmasked points,
        # so zeroing them here is needed.
        ZERO = MASK_ZERO_2D
        if self.runtype > BAROTROPIC_2D:
            ZERO = ZERO + MASK_ZERO_3D
        for v in ZERO:
            array = getattr(self, v)
            edges = array.grid._land & array.grid._water_contact
            array.all_values[...] = 0.0

        self.An.update_halos()
        self.An.all_values[self.domain.T._land] = self.An.fill_value
        if (self.An.all_values[self.domain.T._water] == 0.0).all():
            self.logger.info("Disabling numerical damping because An is 0 everywhere")
        else:
            self.logger.info(
                "Horizontal diffusivity An used for numerical damping ranges between"
                " %s and %s m2 s-1" % (self.An.ma.min(), self.An.ma.max())
            )
            self.An_uu = self.domain.UU.array(fill=np.nan)
            self.An_uv = self.domain.UV.array(fill=np.nan)
            self.An_vu = self.domain.VU.array(fill=np.nan)
            self.An_vv = self.domain.VV.array(fill=np.nan)
            An_x = self.An.interp(self.domain.X)
            self.An_uu.all_values[:, :-1] = self.An.all_values[:, 1:]
            self.An_vv.all_values[:-1, :] = self.An.all_values[1:, :]
            self.An_uv.all_values[:, :] = self.An_vu.all_values[:, :] = An_x.all_values[
                1:, 1:
            ]

    def advance_depth_integrated(
        self,
        timestep: float,
        tausx: core.Array,
        tausy: core.Array,
        dpdx: core.Array,
        dpdy: core.Array,
    ):
        """Update depth-integrated transports (:attr:`U`, :attr:`V`) and depth-averaged
        velocities (:attr:`u1`, :attr:`v1`). This will also update their halos.

        Args:
            timestep: time step (s)
            tausx: surface stress (Pa) in x-direction
            tausy: surface stress (Pa) in y-direction
            dpdx: surface pressure gradient (dimensionless) in x-direction
            dpdx: surface pressure gradient (dimensionless) in y-direction
        """
        # Update 2D transports from t-1/2 to t+1/2.
        # This uses advection, diffusion and bottom friction terms
        # (advU, advV, diffU, diffV, ru, rv) that were calculated by the call
        # to update_2d_momentum_diagnostics
        if self._ufirst:
            self.u_2d(timestep, tausx, dpdx)
            self.U.mirror()
            self.U.update_halos()
            self.coriolis(self.U, self.fU)
            self.v_2d(timestep, tausy, dpdy)
            self.V.mirror()
            self.V.update_halos()
            self.coriolis(self.V, self.fV)
        else:
            self.v_2d(timestep, tausy, dpdy)
            self.V.mirror()
            self.V.update_halos()
            self.coriolis(self.V, self.fV)
            self.u_2d(timestep, tausx, dpdx)
            self.U.mirror()
            self.U.update_halos()
            self.coriolis(self.U, self.fU)
        self._ufirst = not self._ufirst

        self.Ui_tmp += self.U.all_values
        self.Vi_tmp += self.V.all_values

    def update_depth_integrated_diagnostics(
        self, timestep: float, skip_coriolis: bool = False, update_z0b: bool = False
    ):
        """Update 2D momentum diagnostics, including the Coriolis terms that will drive
        the next 2D update.

        Args:
            timestep: time step (s) to calculate advection of momentum over
            skip_coriolis: flag to indicate that Coriolis terms are already up-to-date
                and do not need recomputing, for instance, after a recent call to
                :meth:`advance_depth_integrated`
            update_z0b: whether to iteratively update hydrodynamic bottom roughness
        """
        if not skip_coriolis:
            self.coriolis(self.U, self.fU)
            self.coriolis(self.V, self.fV)

        # Calculate sources of transports U and V due to advection (advU, advV)
        # and diffusion (diffU, diffV)
        # Transports generally come in at time=-1/2 and are then advanced to time+1/2
        self.transport_2d_momentum(
            self.U,
            self.V,
            timestep,
            self.advU,
            self.advV,
            self.diffU,
            self.diffV,
            update_z0b,
        )

        if self.An_uu is not None:
            # Numerical damping of transports
            pygetm._pygetm.horizontal_diffusion(
                self.U, self.An_uu, self.An_uv, timestep, self.dampU
            )
            pygetm._pygetm.horizontal_diffusion(
                self.V, self.An_vu, self.An_vv, timestep, self.dampV
            )

    def advance(
        self,
        timestep: float,
        split_factor: int,
        tausx: core.Array,
        tausy: core.Array,
        dpdx: core.Array,
        dpdy: core.Array,
        idpdx: core.Array,
        idpdy: core.Array,
        viscosity: core.Array,
    ):
        """Update depth-explicit transports (:attr:`pk`, :attr:`qk`) and velocities
        (:attr:`uk`, :attr:`vk`). This will also update their halos.

        Args:
            timestep: (macro) time step (s)
            split_factor: number of microtimesteps per macrotimestep
            tausx: surface stress (Pa) in x-direction
            tausy: surface stress (Pa) in y-direction
            dpdx: surface pressure gradient (dimensionless) in x-direction
            dpdy: surface pressure gradient (dimensionless) in y-direction
            idpdx: internal pressure gradient (m2 s-2) in x-direction
            idpdy: internal pressure gradient (m2 s-2) in y-direction
            viscosity: turbulent viscosity (m2 s-1)
        """

        # Depth-integrated transports have been summed over all microtimesteps.
        # Average them, then reset depth-integrated transports that will be incremented
        # over the next macrotimestep.
        np.multiply(self.Ui_tmp, 1.0 / split_factor, out=self.Ui.all_values)
        np.multiply(self.Vi_tmp, 1.0 / split_factor, out=self.Vi.all_values)
        self.Ui_tmp.fill(0.0)
        self.Vi_tmp.fill(0.0)

        # Do the halo exchange for viscosity, as this needs to be interpolated
        # to the U and V grids. For that, information from the halos is used.
        viscosity.update_halos(parallel.Neighbor.TOP_AND_RIGHT)

        def adv(advancer, taus, dp, idp, tp3d, tp2d, cor, dev):
            grid = tp3d.grid
            advancer(timestep, taus, dp, idp, viscosity.interp(grid))
            tp3d.all_values.sum(axis=0, out=dev.all_values)
            dev.all_values -= tp2d.all_values
            dev.all_values /= grid.D.all_values
            tp3d.all_values -= dev.all_values * grid.hn.all_values
            tp3d.mirror()
            tp3d.update_halos()
            self.coriolis(tp3d, cor)

        # Update horizontal transports. Also update the halos so that transports
        # (and more importantly, the velocities derived subsequently) are valid there.
        # Information from these halos is needed for many reasons:
        # - the Coriolis update requires horizontal velocities at the four points
        #   surrounding each U/V point
        # - to advect the horizontal velocities themselves, for which they need to be
        #   valid in the halos in the direction of transport
        # - to advect quantities defined on the T grid, as this requires horizontal
        #   velocities at the boundaries of every T cell of the subdomain interior;
        #   this includes cells at the very Western and Southern boundary,
        #   which for U and V grids lie within the halo
        # - to allow interpolation of horizontal velocities to the advection grids for
        #   momentum (UU, UV, VU, VV), which again requires halos values
        # - to calculate vertical velocities, which requires horizontal transports at
        #   the four interfaces around every T point
        if self._u3dfirst:
            adv(self.pk_3d, tausx, dpdx, idpdx, self.pk, self.Ui, self.fpk, self.udev)
            adv(self.qk_3d, tausy, dpdy, idpdy, self.qk, self.Vi, self.fqk, self.vdev)
        else:
            adv(self.qk_3d, tausy, dpdy, idpdy, self.qk, self.Vi, self.fqk, self.vdev)
            adv(self.pk_3d, tausx, dpdx, idpdx, self.pk, self.Ui, self.fpk, self.udev)
        self._u3dfirst = not self._u3dfirst

        self.update_diagnostics(timestep, viscosity, skip_coriolis=True)

    def update_diagnostics(
        self, timestep: float, viscosity: core.Array, skip_coriolis: bool = False
    ):
        """Update 3D momentum diagnostics, including the vertical velocity :attr:`ww`,
        the slow terms that will drive the 2D updates over the next macrotimestep,
        and the bottom friction and Coriolis terms that will drive the next 3D update.
        NB the Coriolis update is already done as part of the momentum update itself,
        so needed only when starting from a restart.

        Args:
            timestep: time step (s)
            viscosity: turbulent viscosity (T grid, layer interfaces, m2 s-1)
            skip_coriolis: flag to indicate that Coriolis terms are already up-to-date
                and do not need recomputing
        """
        if not skip_coriolis:
            self.coriolis(self.pk, self.fpk)
            self.coriolis(self.qk, self.fqk)

        # Infer vertical velocity from horizontal transports and desired layer height
        # change (ho -> hn). This is done at all points surrounding U and V points, so
        # no further halo exchange of w is needed to support interpolation to U and V
        # grids later on. This does require that transports are up to date in halos.
        self.w_3d(timestep)

        itimestep = 1.0 / timestep

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by
        # layer heights Both velocities and U/V thicknesses are now at time 1/2
        np.divide(
            self.pk.all_values, self.pk.grid.hn.all_values, out=self.uk.all_values
        )
        np.divide(
            self.qk.all_values, self.qk.grid.hn.all_values, out=self.vk.all_values
        )

        # Use updated velocities (uk, vk) to compute shear frequency (SS) at T points
        # (interior only, not in halos)
        self.update_shear_frequency(viscosity)

        # Calculate bottom friction from updated velocities (and syncronized layer
        # thicknesses hn). This needs to be done before derived quantities such as
        # bottom stress are calculated
        if self.apply_bottom_friction:
            self.bottom_friction(
                self.u_bot, self.v_bot, self.hU_bot, self.hV_bot, self.rru, self.rrv
            )

        # Interpolate 3D velocities to advection grids.
        # This needs to be done before uk/vk are changed by the advection operator.
        self.uk.interp(self.uua3d)
        self.vk.interp(self.uva3d)
        self.uk.interp(self.vua3d)
        self.vk.interp(self.vva3d)

        # Advect 3D u and v velocity from time=1/2 to 1 1/2 using velocities
        # interpolated to its own advection grids. Store the resulting trend, which
        # will be applied as part of the momentum update in the next timestep.
        # They will also be used to calculate the slow advection contribution to
        # depth-integrated momentum equations.
        # JB the alternative would be to interpolate transports and then divide by
        # (colocated) layer heights, like we do for 2D
        self.uadv.apply_3d(
            self.uua3d,
            self.uva3d,
            self.ww.interp(self.uk.grid),
            timestep,
            self.uk,
            new_h=True,
            skip_initial_halo_exchange=True,
        )
        self.advpk.all_values[...] = (
            self.uk.all_values * self.uadv.h - self.pk.all_values
        ) * itimestep
        self.vadv.apply_3d(
            self.vua3d,
            self.vva3d,
            self.ww.interp(self.vk.grid),
            timestep,
            self.vk,
            new_h=True,
            skip_initial_halo_exchange=True,
        )
        self.advqk.all_values[...] = (
            self.vk.all_values * self.vadv.h - self.qk.all_values
        ) * itimestep

        # Restore velocity at time=1/2
        # (the final value at the end of the current timestep)
        np.divide(
            self.pk.all_values, self.pk.grid.hn.all_values, out=self.uk.all_values
        )
        np.divide(
            self.qk.all_values, self.qk.grid.hn.all_values, out=self.vk.all_values
        )

        if self.diffuse_momentum:
            # Calculate the momentum trends (diffpk, diffqk) associated with diffusion
            # of 3D u and v velocity between time=1/2 to 1 1/2. Note that thicknesses
            # should be in sync with velocities uk and vk. This means they should lag
            # 1/2 a timestep behind the T grid (already the case for X, but for T we
            # use 1/2(ho+hn))
            self.momentum_diffusion_driver(
                self.domain.h_T_half,
                self.domain.X.hn,
                self.uk,
                self.vk,
                self.diffpk,
                self.diffqk,
            )

        # Compute slow (3D) advection and diffusion contribution to to the
        # depth-integrated momentum equations. This is done by comparing the
        # depth-integrated 3D transport calculated above (between centers of the
        # current and next macrotime step) with the newly calculated depth-integrated
        # transport based on accumulated 2D transports (accumulated over the
        # current macrotimestep, and thus representative for its center).
        self.transport_2d_momentum(
            self.Ui, self.Vi, timestep, self.SxA, self.SyA, self.SxD, self.SyD, False
        )
        self.SxA.all_values[...] = (
            self.advpk.all_values.sum(axis=0) - self.SxA.all_values
        )
        self.SyA.all_values[...] = (
            self.advqk.all_values.sum(axis=0) - self.SyA.all_values
        )
        self.SxD.all_values[...] = (
            self.diffpk.all_values.sum(axis=0) - self.SxD.all_values
        )
        self.SyD.all_values[...] = (
            self.diffqk.all_values.sum(axis=0) - self.SyD.all_values
        )

        if self.apply_bottom_friction:
            # Note: ru and rv have been updated by transport_2d_momentum, using
            # accumulated transports Ui and Vi (representative for t=1/2, just like uk,
            # vk, rru, rrv). Slow bottom friction (stress/density) is derived by taking
            # the difference between 3D bottom friction and the inferred ru and rv.
            self.SxF.all_values[...] = (
                -self.rru.all_values * self.uk.all_values[0, ...]
                + self.ru.all_values * self.u1.all_values
            )
            self.SyF.all_values[...] = (
                -self.rrv.all_values * self.vk.all_values[0, ...]
                + self.rv.all_values * self.v1.all_values
            )

    def transport_2d_momentum(
        self,
        U: core.Array,
        V: core.Array,
        timestep: float,
        advU: core.Array,
        advV: core.Array,
        diffU: core.Array,
        diffV: core.Array,
        update_z0b: bool,
    ):
        """Advect and optionally diffuse depth-integrated transports in x and y
        direction (arguments ``U`` and ``V``). From these, first the depth-averaged
        velocities are calculated and stored in :attr:`u1` and :attr:`v1`.
        This routine also updates bottom friction :attr:`ru` and  :attr:`rv`.

        Args:
            U: depth-integrated velocity (m2 s-1) in x-direction
            V: depth-integrated velocity (m2 s-1) in y-direction
            timestep: time step (s) to calculate advection over
            advU: array for storing the change in transport ``U`` (m2 s-2)
                due to advection
            advV: array for storing the change in transport ``V`` (m2 s-2)
                due to advection
            diffU: array for storing the change in transport ``U`` (m2 s-2)
                due to diffusion
            diffV: array for storing the change in transport ``V`` (m2 s-2)
                due to diffusion
            update_z0b: whether to iteratively update hydrodynamic bottom roughness
        """
        np.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
        np.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)

        if self.diffuse_momentum:
            # Compute velocity diffusion contribution to transport sources.
            # This uses depth-averaged velocities u1 and v1, which therefore have to be
            # up to date. Water depths should be in sync with velocities, which means
            # they should lag 1/2 a timestep behind the tracer/T grid
            self.momentum_diffusion_driver(
                self.domain.D_T_half, self.domain.X.D, self.u1, self.v1, diffU, diffV
            )

        # Calculate bottom friction (ru and rv) using updated depth-averaged velocities
        # u1 and v1. Warning: this uses velocities u1 and v1 at masked points, which
        # therefore need to be kept at 0
        if self.apply_bottom_friction:
            self.bottom_friction(
                self.u1,
                self.v1,
                self.domain.U.D,
                self.domain.V.D,
                self.ru,
                self.rv,
                update_z0b,
            )

        itimestep = 1.0 / timestep

        # Advection of u velocity (u1)
        U.interp(self.uua)
        V.interp(self.uva)
        self.uua.all_values /= self.domain.UU.D.all_values
        self.uva.all_values /= self.domain.UV.D.all_values
        self.uadv(
            self.uua, self.uva, timestep, self.u1, skip_initial_halo_exchange=True,
        )
        advU.all_values[...] = (
            self.u1.all_values * self.uadv.D - U.all_values
        ) * itimestep

        # Advection of v velocity (v1)
        U.interp(self.vua)
        V.interp(self.vva)
        self.vua.all_values /= self.domain.VU.D.all_values
        self.vva.all_values /= self.domain.VV.D.all_values
        self.vadv(
            self.vua, self.vva, timestep, self.v1, skip_initial_halo_exchange=True,
        )
        advV.all_values[...] = (
            self.v1.all_values * self.vadv.D - V.all_values
        ) * itimestep

        # Restore depth-averaged velocities as they need to be valid on exit
        np.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
        np.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)

    def bottom_friction(
        self,
        u: core.Array,
        v: core.Array,
        DU: core.Array,
        DV: core.Array,
        ru: core.Array,
        rv: core.Array,
        update_z0b: bool = False,
    ):
        """Calculate bottom friction
        
        Args:
            u: velocity in x-direction (m s-1) in bottom layer [U grid]
            v: velocity in y-direction (m s-1) in bottom layer [V grid]
            DU: thickness (m) of bottom layer on the U grid
            DV: thickness (m) of bottom layer on the V grid
            ru: the square of friction velocity (= stress in Pa divided by density),
                per actual velocity at the layer center on the U grid
            rv: the square of friction velocity (= stress in Pa divided by density),
                per actual velocity at the layer center on the V grid
            update_z0b: whether to iteratively update the hydrodynamic bottom roughness
        """
        u.interp(self.u_V)
        v.interp(self.v_U)
        pygetm._pygetm.bottom_friction(u, self.v_U, DU, ru, update_z0b)
        pygetm._pygetm.bottom_friction(self.u_V, v, DV, rv, update_z0b)

    def coriolis(self, U: core.Array, out: core.Array):
        """Calculate change in transport due to Coriolis force

        Args:
            U: 2d or 3d transport in x- or y-direction (m2 s-1)
            out: array to store change in complementary transport (y or x-direction)
                due to Coriolis force (m2 s-2)
        """
        if self.coriolis_scheme == CoriolisScheme.ESPELID:
            Din = U.grid.D if U.ndim == 2 else U.grid.hn
            Dout = out.grid.D if out.ndim == 2 else out.grid.hn
            sqrtDin = np.sqrt(Din, where=U.grid._water)
            sqrtDin.all_values[..., U.grid._land] = 1.0
            (U / sqrtDin).interp(out)
            out.all_values *= np.sqrt(Dout.all_values, where=out.grid._water)
        else:
            U.interp(out)
        # todo: cord_curv
        out.all_values *= out.grid.cor.all_values


# Expose all Fortran arrays that are a member of Momentum as read-only properties
# The originals are members with and underscore as prefix, and therefore not visible to
# the user. This ensures the user will not accidentally disconnect the Python variable
# from the underlying Fortran libraries/data
for membername in Momentum._all_fortran_arrays:
    attrs = Momentum._array_args.get(membername[1:], {})
    long_name = attrs.get("long_name")
    units = attrs.get("units")
    doc = ""
    if long_name is not None:
        doc = long_name
        if units:
            doc += " (%s)" % units
    prop = property(operator.attrgetter(membername), doc=doc)
    setattr(Momentum, membername[1:], prop)
