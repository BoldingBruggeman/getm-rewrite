import enum
import operator
import logging
from typing import Union

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
    ESPELID = 2  #: Espelid et al. [2000], IJNME 49, 1521-1545


class Momentum(pygetm._pygetm.Momentum):
    _arrays = (
        "U",
        "V",
        "fU",
        "fV",
        "advU",
        "advV",
        "diffu1",
        "diffv1",
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
        "An",
    )

    _array_args = {
        "U": dict(
            units="m2 s-1",
            long_name="depth-integrated transport in x direction",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "V": dict(
            units="m2 s-1",
            long_name="depth-integrated transport in y direction",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "Ui": dict(
            units="m2 s-1",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "Vi": dict(
            units="m2 s-1",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "u1": dict(
            units="m s-1",
            long_name="depth-averaged velocity in x direction",
            fill_value=FILL_VALUE,
            attrs={"_mask_output": True},
        ),
        "v1": dict(
            units="m s-1",
            long_name="depth-averaged velocity in y direction",
            fill_value=FILL_VALUE,
            attrs={"_mask_output": True},
        ),
        "pk": dict(
            units="m2 s-1",
            long_name="layer-integrated transport in x direction",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "qk": dict(
            units="m2 s-1",
            long_name="layer-integrated transport in y direction",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True, "_mask_output": True},
        ),
        "uk": dict(
            units="m s-1",
            long_name="velocity in x direction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, standard_name="sea_water_x_velocity"),
        ),
        "vk": dict(
            units="m s-1",
            long_name="velocity in y direction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, standard_name="sea_water_y_velocity"),
        ),
        "ww": dict(
            units="m s-1",
            long_name="vertical velocity",
            fill_value=FILL_VALUE,
            attrs=dict(standard_name="upward_sea_water_velocity"),
        ),
        "SS": dict(
            units="s-2", long_name="shear frequency squared", fill_value=FILL_VALUE
        ),
    }

    def __init__(
        self,
        logger: logging.Logger,
        domain: pygetm.domain.Domain,
        runtype: int,
        apply_bottom_friction: bool = True,
        Am: float = 0.0,
        An: Union[float, core.Array] = 0.0,
        cnpar: float = 1.0,
        advection_scheme: operators.AdvectionScheme = operators.AdvectionScheme.HSIMT,
        coriolis_scheme: CoriolisScheme = CoriolisScheme.DEFAULT,
    ):
        super().__init__(domain, runtype, Am, cnpar, coriolis_scheme)

        for name in self._arrays:
            setattr(
                self,
                "_%s" % name,
                self.wrap(
                    core.Array(name=name, **self._array_args.get(name, {})),
                    name.encode("ascii"),
                ),
            )

        self.logger = logger
        self.runtype = runtype

        # Disable bottom friction if physical bottom roughness is 0 everywhere
        if (
            apply_bottom_friction
            and (np.ma.array(domain.z0b_min, mask=domain.mask == 0) == 0.0).any()
        ):
            self.logger.warning(
                "Disabling bottom friction because bottom roughness is 0"
                " in one or more points."
            )
            apply_bottom_friction = False
        self.apply_bottom_friction = apply_bottom_friction
        self.diffuse_momentum = Am > 0.0
        if not self.diffuse_momentum:
            self.logger.info("Diffusion of momentum is off because Am is 0")

        self.U.all_values.fill(0.0)
        self.V.all_values.fill(0.0)
        self.u1.all_values.fill(0.0)
        self.v1.all_values.fill(0.0)
        self.Ui.all_values.fill(0.0)
        self.Vi.all_values.fill(0.0)
        self.Ui_tmp = np.zeros_like(self.Ui.all_values)
        self.Vi_tmp = np.zeros_like(self.Vi.all_values)
        if runtype > BAROTROPIC_2D:
            self.pk.all_values.fill(0.0)
            self.qk.all_values.fill(0.0)
            self.uk.all_values.fill(0.0)
            self.vk.all_values.fill(0.0)
            self.ww.all_values.fill(0.0)
            self.SS.fill(0.0)  # for surface/bottom interfaces, which are not updated

        self.uadv = operators.Advection(domain.U, scheme=advection_scheme)
        self.vadv = operators.Advection(domain.V, scheme=advection_scheme)

        self.uua = domain.UU.array(fill=np.nan)
        self.uva = domain.UV.array(fill=np.nan)
        self.vua = domain.VU.array(fill=np.nan)
        self.vva = domain.VV.array(fill=np.nan)

        self.uua3d = domain.UU.array(fill=np.nan, z=CENTERS)
        self.uva3d = domain.UV.array(fill=np.nan, z=CENTERS)
        self.vua3d = domain.VU.array(fill=np.nan, z=CENTERS)
        self.vva3d = domain.VV.array(fill=np.nan, z=CENTERS)

        self.An = domain.T.array(
            name="An",
            units="m2 s-1",
            long_name="horizontal diffusivity of momentum",
            fill_value=FILL_VALUE,
            attrs=dict(_require_halos=True),
        )
        self.An.fill(An)

        #: Whether to start the depth-integrated (2D) momentum update with u
        # (as opposed to v)
        self._ufirst = False

        #: Whether to start the depth-explicit (3D) momentum update with u
        # (as opposed to v)
        self._u3dfirst = False

    def start(self):
        # Ensure transports and velocities are 0 in masked points
        # NB velocities will be computed from transports, but only in unmasked points,
        # so zeroing them here is needed.
        zero_masked = [self.U, self.V, self.u1, self.v1]
        if self.runtype > BAROTROPIC_2D:
            zero_masked += [
                self.Ui,
                self.Vi,
                self.pk,
                self.qk,
                self.uk,
                self.vk,
                self.ww,
            ]
        for array in zero_masked:
            array.all_values[..., array.grid.mask.all_values == 0] = 0.0
        if (self.An.ma == 0.0).all():
            self.logger.info("Disabling An because it is 0 everywhere")
            self.An = None

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
            tausx: surface stress (Pa) in x direction
            tausy: surface stress (Pa) in y direction
            dpdx: surface pressure gradient (dimensionless) in x direction
            dpdx: surface pressure gradient (dimensionless) in y direction
        """
        # Update 2D transports from t-1/2 to t+1/2.
        # This uses advection, diffusion and bottom friction terms
        # (advU, advV, diffu1, diffv1, ru, rv) that were calculated by the call
        # to update_2d_momentum_diagnostics
        if self._ufirst:
            self.u_2d(timestep, tausx, dpdx)
            self.U.update_halos()
            self.coriolis_fu()
            self.v_2d(timestep, tausy, dpdy)
            self.V.update_halos()
            self.coriolis_fv()
        else:
            self.v_2d(timestep, tausy, dpdy)
            self.V.update_halos()
            self.coriolis_fv()
            self.u_2d(timestep, tausx, dpdx)
            self.U.update_halos()
            self.coriolis_fu()
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
        """
        if not skip_coriolis:
            self.coriolis_fu()
            self.coriolis_fv()

        # Calculate sources of transports U and V due to advection (advU, advV)
        # and diffusion (diffu1, diffv1)
        # Transports generally come in at time=-1/2 and are then advanced to time+1/2
        self.transport_2d_momentum(
            self.U,
            self.V,
            timestep,
            self.advU,
            self.advV,
            self.diffu1,
            self.diffv1,
            update_z0b,
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
            tausx: surface stress (Pa) in x direction
            tausy: surface stress (Pa) in y direction
            dpdx: surface pressure gradient (dimensionless) in x direction
            dpdy: surface pressure gradient (dimensionless) in y direction
            idpdx: internal pressure gradient (m2 s-2) in x direction
            idpdy: internal pressure gradient (m2 s-2) in y direction
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
            self.pk_3d(timestep, tausx, dpdx, idpdx, viscosity.interp(self.domain.U))
            self.pk.update_halos()
            self.coriolis_fpk()
            self.qk_3d(timestep, tausy, dpdy, idpdy, viscosity.interp(self.domain.V))
            self.qk.update_halos()
            self.coriolis_fqk()
        else:
            self.qk_3d(timestep, tausy, dpdy, idpdy, viscosity.interp(self.domain.V))
            self.qk.update_halos()
            self.coriolis_fqk()
            self.pk_3d(timestep, tausx, dpdx, idpdx, viscosity.interp(self.domain.U))
            self.pk.update_halos()
            self.coriolis_fpk()
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
            self.coriolis_fpk()
            self.coriolis_fqk()

        # Infer vertical velocity from horizontal transports and desired layer height
        # change (ho -> hn). This is done at all points surrounding U and V points, so
        # no further halo exchange of w is needed to support interpolation to U and V
        # grids later on. This does require that transports are up to date in halos.
        self.w_3d(timestep)

        itimestep = 1.0 / timestep

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by
        # layer heights Both velocities and U/V thicknesses are now at time 1/2
        np.divide(
            self.pk.all_values,
            self.pk.grid.hn.all_values,
            where=self.pk.grid.mask.all_values != 0,
            out=self.uk.all_values,
        )
        np.divide(
            self.qk.all_values,
            self.qk.grid.hn.all_values,
            where=self.qk.grid.mask.all_values != 0,
            out=self.vk.all_values,
        )

        # Use updated velocities (uk, vk) to compute shear frequency (SS) at T points
        # (interior only, not in halos)
        self.update_shear_frequency(viscosity)

        # Calculate bottom friction from updated velocities (and syncronized layer
        # thicknesses hn). This needs to be done before derived quantities such as
        # bottom stress are calculated
        if self.apply_bottom_friction:
            self.bottom_friction_3d()

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
            Ah=self.An,
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
            Ah=self.An,
            new_h=True,
            skip_initial_halo_exchange=True,
        )
        self.advqk.all_values[...] = (
            self.vk.all_values * self.vadv.h - self.qk.all_values
        ) * itimestep

        # Restore velocity at time=1/2
        # (the final value at the end of the current timestep)
        np.divide(
            self.pk.all_values,
            self.pk.grid.hn.all_values,
            where=self.pk.grid.mask.all_values != 0,
            out=self.uk.all_values,
        )
        np.divide(
            self.qk.all_values,
            self.qk.grid.hn.all_values,
            where=self.qk.grid.mask.all_values != 0,
            out=self.vk.all_values,
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
            U: depth-integrated velocity (m2 s-1) in x direction
            V: depth-integrated velocity (m2 s-1) in y direction
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
            self.bottom_friction_2d(update_z0b)

        itimestep = 1.0 / timestep

        # Advection of u velocity (u1)
        U.interp(self.uua)
        V.interp(self.uva)
        self.uua.all_values /= self.domain.UU.D.all_values
        self.uva.all_values /= self.domain.UV.D.all_values
        self.uadv(
            self.uua,
            self.uva,
            timestep,
            self.u1,
            Ah=self.An,
            skip_initial_halo_exchange=True,
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
            self.vua,
            self.vva,
            timestep,
            self.v1,
            Ah=self.An,
            skip_initial_halo_exchange=True,
        )
        advV.all_values[...] = (
            self.v1.all_values * self.vadv.D - V.all_values
        ) * itimestep

        # Restore depth-averaged velocities as they need to be valid on exit
        np.divide(U.all_values, U.grid.D.all_values, out=self.u1.all_values)
        np.divide(V.all_values, V.grid.D.all_values, out=self.v1.all_values)


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
