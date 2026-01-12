import enum
import logging
from typing import Optional

import numpy as np

from . import core
from . import parallel
from . import operators
import pygetm._pygetm
from .constants import (
    FILL_VALUE,
    RunType,
    CENTERS,
    RHO0,
    INTERFACES,
    TimeVarying,
    CellType,
)


class CoriolisScheme(enum.IntEnum):
    OFF = 0
    DEFAULT = 1
    ESPELID = 2  #: `Espelid et al. (2000) <https://doi.org/10.1002/1097-0207(20001230)49%3A12%3C1521%3A%3AAID-NME9%3E3.0.CO%3B2-F>`_


MASK_ZERO_2D = ("U", "V", "u1", "v1")
MASK_ZERO_3D = ("pk", "qk", "uk", "vk", "Ui", "Vi", "rru", "rrv")


class Momentum:
    _arrays = (
        "U",
        "V",
        "corU",
        "corV",
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
        "corpk",
        "corqk",
        "SxB",
        "SyB",
        "SxA",
        "SyA",
        "SxD",
        "SyD",
        "SxF",
        "SyF",
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
    __slots__ = _arrays + (
        "runtype",
        "logger",
        "_U_cum",
        "_V_cum",
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
        "_vertical_diffusion",
        "ea2",
        "ea4",
        "avmmol",
        "ustar2_bx",
        "ustar2_by",
        "ustar_b",
        "taub",
    )

    def __init__(
        self,
        Am: float = 0.0,
        An: float = 0.0,
        cnpar: float = 1.0,
        advection_scheme: Optional[operators.AdvectionScheme] = None,
        advection_split_2d: operators.AdvectionSplit = operators.AdvectionSplit.FULL,
        coriolis_scheme: CoriolisScheme = CoriolisScheme.DEFAULT,
        avmmol: float = 1.8e-6,
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
            avmmol: molecular viscosity (m2 s-1)
        """
        self._Am_const = Am
        self._An_const = An
        self.cnpar = cnpar
        self.advection_scheme = advection_scheme
        self.advection_split_2d = advection_split_2d
        self.coriolis_scheme = coriolis_scheme
        self.avmmol = avmmol

    def _create_depth_integrated_arrays(self, tgrid: core.Grid):
        ugrid = tgrid.ugrid
        vgrid = tgrid.vgrid

        self.U = ugrid.array(
            name="U",
            units="m2 s-1",
            long_name="depth-integrated velocity in x-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_Y),
            ),
        )

        self.V = vgrid.array(
            name="V",
            units="m2 s-1",
            long_name="depth-integrated velocity in y-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_X),
            ),
        )
        self.u1 = ugrid.array(
            name="u1",
            units="m s-1",
            long_name="depth-averaged velocity in x-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_Y),
            ),
        )
        self.v1 = vgrid.array(
            name="v1",
            units="m s-1",
            long_name="depth-averaged velocity in y-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_X),
            ),
        )
        self.SxA = ugrid.array(
            name="SxA",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to slow advection",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SyA = vgrid.array(
            name="SyA",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to slow advection",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SxB = ugrid.array(
            name="SxB",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to internal pressure",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SyB = vgrid.array(
            name="SyB",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to internal pressure",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SxD = ugrid.array(
            name="SxD",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to slow diffusion",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SyD = vgrid.array(
            name="SyD",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to slow diffusion",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SxF = ugrid.array(
            name="SxF",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to slow bottom friction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.SyF = vgrid.array(
            name="SyF",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to slow bottom friction",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.corU = ugrid.array(
            name="corU",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to Coriolis effect",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.corV = vgrid.array(
            name="corV",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to Coriolis effect",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.advU = ugrid.array(
            name="advU",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to advection",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.advV = vgrid.array(
            name="advV",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to advection",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.diffU = ugrid.array(
            name="diffU",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to horizontal diffusion",
            fill_value=FILL_VALUE,
        )
        self.diffV = vgrid.array(
            name="diffV",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to horizontal diffusion",
            fill_value=FILL_VALUE,
        )
        self.dampU = ugrid.array(
            name="dampU",
            units="m-2 s-2",
            long_name="tendency of depth-integrated x-velocity due to numerical damping",
            fill_value=FILL_VALUE,
        )
        self.dampV = vgrid.array(
            name="dampV",
            units="m-2 s-2",
            long_name="tendency of depth-integrated y-velocity due to numerical damping",
            fill_value=FILL_VALUE,
        )
        self.ru = ugrid.array(
            name="ru",
            units="m s-1",
            long_name="bottom drag coefficient multiplied by depth-averaged velocity",
            fill_value=FILL_VALUE,
        )
        self.rv = vgrid.array(
            name="rv",
            units="m s-1",
            long_name="bottom drag coefficient multiplied by depth-averaged velocity",
            fill_value=FILL_VALUE,
        )

        self.uua = ugrid.ugrid.array(fill=np.nan)
        self.uva = ugrid.vgrid.array(fill=np.nan)
        self.vua = vgrid.ugrid.array(fill=np.nan)
        self.vva = vgrid.vgrid.array(fill=np.nan)

    def _create_depth_explicit_arrays(self, tgrid: core.Grid):
        ugrid = tgrid.ugrid
        vgrid = tgrid.vgrid

        self.pk = ugrid.array(
            name="pk",
            z=CENTERS,
            units="m2 s-1",
            long_name="layer-integrated velocity in x-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_Y),
            ),
        )
        self.qk = vgrid.array(
            name="qk",
            z=CENTERS,
            units="m2 s-1",
            long_name="layer-integrated velocity in y-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_X),
            ),
        )
        self.uk = ugrid.array(
            name="uk",
            z=CENTERS,
            units="m s-1",
            long_name="velocity in x-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                standard_name="sea_water_x_velocity",
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_Y),
            ),
        )
        self.vk = vgrid.array(
            name="vk",
            z=CENTERS,
            units="m s-1",
            long_name="velocity in y-direction",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                standard_name="sea_water_y_velocity",
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_X),
            ),
        )
        self.ww = tgrid.array(
            name="ww",
            z=INTERFACES,
            units="m s-1",
            long_name="vertical velocity",
            fill_value=FILL_VALUE,
            attrs=dict(
                standard_name="upward_sea_water_velocity",
                _valid_at=(CellType.BOUNDARY,),
            ),
        )
        self.SS = tgrid.array(
            name="SS",
            z=INTERFACES,
            units="s-2",
            long_name="shear frequency squared",
            fill_value=FILL_VALUE,
        )
        self.corpk = ugrid.array(
            name="corpk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated x-velocity due to Coriolis effect",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.corqk = vgrid.array(
            name="corqk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated y-velocity due to Coriolis effect",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.advpk = ugrid.array(
            name="advpk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated x-velocity due to advection",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.advqk = vgrid.array(
            name="advqk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated y-velocity due to advection",
            z=CENTERS,
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True),
        )
        self.diffpk = ugrid.array(
            name="diffpk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated x-velocity due to horizontal diffusion",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.diffqk = vgrid.array(
            name="diffqk",
            units="m-2 s-2",
            long_name="tendency of layer-integrated y-velocity due to horizontal diffusion",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.rru = ugrid.array(
            name="rru",
            units="m s-1",
            long_name="bottom drag coefficient multiplied by near-bottom velocity",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                _time_varying=TimeVarying.MACRO,
                _valid_at=(
                    CellType.EDGE_Y,
                ),  # must be finite to calculate stresses at T points
            ),
        )
        self.rrv = vgrid.array(
            name="rrv",
            units="m s-1",
            long_name="bottom drag coefficient multiplied by near-bottom velocity",
            fill_value=FILL_VALUE,
            attrs=dict(
                _mask_output=True,
                _time_varying=TimeVarying.MACRO,
                _valid_at=(
                    CellType.EDGE_X,
                ),  # must be finite to calculate stresses at T points
            ),
        )

        self.Ui = ugrid.array(
            name="Ui",
            units="m2 s-1",
            long_name="depth-integrated velocity in x-direction averaged over previous macrotimestep",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _time_varying=TimeVarying.MACRO,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_Y),
            ),
        )
        self.Vi = vgrid.array(
            name="Vi",
            units="m2 s-1",
            long_name="depth-integrated velocity in y-direction averaged over previous macrotimestep",
            fill_value=FILL_VALUE,
            attrs=dict(
                _part_of_state=True,
                _mask_output=True,
                _time_varying=TimeVarying.MACRO,
                _valid_at=(CellType.MIRROR_INT, CellType.MIRROR_EXT, CellType.EDGE_X),
            ),
        )

        self._U_cum = np.zeros_like(self.Ui.all_values)
        self._V_cum = np.zeros_like(self.Vi.all_values)

        self.udev = ugrid.array(
            name="udev",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, _time_varying=TimeVarying.MACRO),
        )
        self.vdev = vgrid.array(
            name="vdev",
            units="m s-1",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, _time_varying=TimeVarying.MACRO),
        )

        self.uua3d = ugrid.ugrid.array(fill=np.nan, z=CENTERS)
        self.uva3d = ugrid.vgrid.array(fill=np.nan, z=CENTERS)
        self.vua3d = vgrid.ugrid.array(fill=np.nan, z=CENTERS)
        self.vva3d = vgrid.vgrid.array(fill=np.nan, z=CENTERS)

        self.ustar2_bx = ugrid.array(
            name="ustar2_bx",
            long_name="tendency of depth-integrated x-velocity due to bottom friction",
            units="m2 s-2",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, _time_varying=TimeVarying.MACRO),
        )
        self.ustar2_by = vgrid.array(
            name="ustar2_by",
            long_name="tendency of depth-integrated y-velocity due to bottom friction",
            units="m2 s-2",
            fill_value=FILL_VALUE,
            attrs=dict(_mask_output=True, _time_varying=TimeVarying.MACRO),
        )
        self.ustar_b = tgrid.array(
            name="ustar_b",
            units="m s-1",
            long_name="bottom shear velocity",
            fill_value=FILL_VALUE,
            attrs=dict(_time_varying=TimeVarying.MACRO),
        )
        self.taub = tgrid.array(
            fill=0.0,
            name="taub",
            units="Pa",
            long_name="bottom shear stress",
            fill_value=FILL_VALUE,
            fabm_standard_name="bottom_stress",
            attrs=dict(_mask_output=True, _time_varying=TimeVarying.MACRO),
        )

    def initialize(
        self,
        logger: logging.Logger,
        tgrid: core.Grid,
        runtype: RunType,
        default_advection_scheme: operators.AdvectionScheme,
    ):
        self._create_depth_integrated_arrays(tgrid)
        if runtype > RunType.BAROTROPIC_2D:
            self._create_depth_explicit_arrays(tgrid)

        ugrid = tgrid.ugrid
        vgrid = tgrid.vgrid

        self.logger = logger
        self.runtype = runtype

        # Disable bottom friction if minimum hydrodynamic bottom roughness is 0 everywhere
        z0b_u = ugrid.z0b_min.all_values[~ugrid.z0b_min.all_mask] == 0.0
        z0b_v = vgrid.z0b_min.all_values[~vgrid.z0b_min.all_mask] == 0.0
        self.apply_bottom_friction = not (z0b_u.any() or z0b_v.any())
        if not self.apply_bottom_friction:
            self.logger.warning(
                "Disabling bottom friction because minimum bottom roughness is 0"
                f" in {z0b_u.sum()} U points, {z0b_v.sum()} V points"
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
            "corU",
            "corV",
        )
        if runtype > RunType.BAROTROPIC_2D:
            ZERO_EVERYWHERE = ZERO_EVERYWHERE + MASK_ZERO_3D
            ZERO_UNMASKED = ZERO_UNMASKED + (
                "SS",
                "advpk",
                "advqk",
                "diffpk",
                "diffqk",
                "corpk",
                "corqk",
                "ww",
            )
        for v in ZERO_EVERYWHERE:
            array = getattr(self, v)
            array.all_values = 0.0
        for v in ZERO_UNMASKED:
            getattr(self, v).fill(0.0)

        # Zero elements of Coriolis terms that will never be set, but still
        # multiplied by f. This prevents overflow warnings
        # We could zero the whole arrays, but by doing this selectively, we
        # can verify in tests that values on other masked points are not used
        self.corU.all_values[:, -1] = 0.0
        self.corU.all_values[0, :] = 0.0
        self.corV.all_values[:, 0] = 0.0
        self.corV.all_values[-1, :] = 0.0
        if runtype > RunType.BAROTROPIC_2D:
            self.corpk.all_values[:, :, -1] = 0.0
            self.corpk.all_values[:, 0, :] = 0.0
            self.corqk.all_values[:, :, 0] = 0.0
            self.corqk.all_values[:, -1, :] = 0.0

        # Advection operators
        if self.advection_scheme is None:
            self.advection_scheme = default_advection_scheme
        self.logger.info(f"Advection scheme: {self.advection_scheme.name}")
        self.uadv = operators.Advection(
            ugrid, scheme=self.advection_scheme, split_2d=self.advection_split_2d
        )
        self.vadv = operators.Advection(
            vgrid, scheme=self.advection_scheme, split_2d=self.advection_split_2d
        )

        if self.runtype > RunType.BAROTROPIC_2D:
            # Layer thicknesses and velocities in bottom layer (for bottom friction)
            self.hU_bot = ugrid.hn.isel(z=0)
            self.hV_bot = vgrid.hn.isel(z=0)
            self.u_bot = self.uk.isel(z=0)
            self.v_bot = self.vk.isel(z=0)

            # Vertical diffusion operator and helper arrays for source terms
            self.logger.info(
                f"Crank-Nicolson parameter for vertical diffusion: {self.cnpar}"
            )
            self._vertical_diffusion = operators.VerticalDiffusion(
                tgrid, cnpar=self.cnpar
            )
            self.ea2 = tgrid.array(fill=0.0, z=CENTERS)
            self.ea4 = tgrid.array(fill=0.0, z=CENTERS)

        self.logger.info(f"Coriolis interpolation: {self.coriolis_scheme.name}")

        self.An = tgrid.array(
            name="An",
            units="m2 s-1",
            long_name="horizontal diffusivity of momentum",
            fill_value=FILL_VALUE,
            attrs=dict(
                _require_halos=True,
                _time_varying=False,
                _valid_at=(CellType.BOUNDARY,),
            ),
        )
        self.An.fill(self._An_const)

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
        ZERO = MASK_ZERO_2D
        if self.runtype > RunType.BAROTROPIC_2D:
            ZERO = ZERO + MASK_ZERO_3D
        for v in ZERO:
            array = getattr(self, v)
            array.all_values[..., array.grid._land] = 0.0

        # Mirror transports (restarts may not include mirrored points)
        # This is done twice because some external mirrored points (MIRROR_EXT)
        # are mirrored from internal mirrored points (MIRROR_INT)
        for _ in range(2):
            self.U.mirror()
            self.V.mirror()
            if self.runtype > RunType.BAROTROPIC_2D:
                self.pk.mirror()
                self.qk.mirror()
                self.Ui.mirror()
                self.Vi.mirror()

        # If horizontal diffusivity [damping] is specified, interpolate it to
        # U and V grids
        self.An.update_halos()
        self.An.all_values[self.An.all_mask] = self.An.fill_value
        if self.An.all_values.any(where=~self.An.all_mask):
            self.logger.info(
                "Horizontal diffusivity An used for numerical damping ranges between"
                f" {self.An.ma.min()} and {self.An.ma.max()} m2 s-1"
            )
            self.An_uu = self.An.interp(self.U.grid.ugrid.array(fill=np.nan))
            self.An_uv = self.An.interp(self.U.grid.vgrid.array(fill=np.nan))
            self.An_vu = self.An.interp(self.V.grid.ugrid.array(fill=np.nan))
            self.An_vv = self.An.interp(self.V.grid.vgrid.array(fill=np.nan))
        else:
            self.logger.info("Disabling numerical damping because An is 0 everywhere")
            self.An_uu = self.An_uv = self.An_vu = self.An_vv = None

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
            dpdy: surface pressure gradient (dimensionless) in y-direction
        """

        def u():
            pygetm._pygetm.advance_2d_transport(
                self.U,
                dpdx,
                tausx,
                self.corU,
                self.advU,
                self.diffU,
                self.dampU,
                self.SxA,
                self.SxB,
                self.SxD,
                self.SxF,
                self.ru,
                timestep,
            )
            self.U.mirror()
            self.U.update_halos()
            self.coriolis(self.U, self.corV, True)

        def v():
            pygetm._pygetm.advance_2d_transport(
                self.V,
                dpdy,
                tausy,
                self.corV,
                self.advV,
                self.diffV,
                self.dampV,
                self.SyA,
                self.SyB,
                self.SyD,
                self.SyF,
                self.rv,
                timestep,
            )
            self.V.mirror()
            self.V.update_halos()
            self.coriolis(self.V, self.corU, False)

        # Update 2D transports from t-1/2 to t+1/2.
        # This uses advection, diffusion, damping and bottom friction terms
        # (advU, advV, diffU, diffV, dampU, dampV, ru, rv) that were calculated
        # in update_depth_integrated_diagnostics, and slow terms (SxA, SyA, SxB,
        # SyB, SxD, SyD, SxF, SyF) that were calculated in update_diagnostics
        if self._ufirst:
            u()
            v()
        else:
            v()
            u()
        self._ufirst = not self._ufirst

        if self.runtype > RunType.BAROTROPIC_2D:
            self._U_cum += self.U.all_values
            self._V_cum += self.V.all_values

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
            self.coriolis(self.U, self.corV, True)
            self.coriolis(self.V, self.corU, False)

        # Calculate sources of transports U and V due to advection (advU, advV)
        # and diffusion (diffU, diffV)
        # Transports generally come in at time=-1/2 and are then advanced to time+1/2
        self._transport_2d_momentum(
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
        np.multiply(self._U_cum, 1.0 / split_factor, out=self.Ui.all_values)
        np.multiply(self._V_cum, 1.0 / split_factor, out=self.Vi.all_values)
        self._U_cum.fill(0.0)
        self._V_cum.fill(0.0)

        # Do the halo exchange for viscosity, as this needs to be interpolated
        # to the U and V grids. For that, information from the halos is used.
        viscosity.update_halos(parallel.Neighbor.TOP_AND_RIGHT)

        def u():
            self.advance_3d_transport(
                self.pk,
                self.Ui,
                self.uk,
                dpdx,
                self.corpk,
                self.advpk,
                self.diffpk,
                idpdx,
                tausx,
                self.rru,
                viscosity,
                self.udev,
                timestep,
            )
            self.coriolis(self.pk, self.corqk, True)

        def v():
            self.advance_3d_transport(
                self.qk,
                self.Vi,
                self.vk,
                dpdy,
                self.corqk,
                self.advqk,
                self.diffqk,
                idpdy,
                tausy,
                self.rrv,
                viscosity,
                self.vdev,
                timestep,
            )
            self.coriolis(self.qk, self.corpk, False)

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
            u()
            v()
        else:
            v()
            u()
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
            self.coriolis(self.pk, self.corqk, True)
            self.coriolis(self.qk, self.corpk, False)

        # Infer vertical velocity from horizontal transports and desired layer height
        # change (ho -> hn). This is done at all points surrounding U and V points, so
        # no further halo exchange of w is needed to support interpolation to U and V
        # grids later on. This does require that transports are up to date in halos.
        pygetm._pygetm.w_momentum_3d(self.pk, self.qk, timestep, self.ww)

        # Compute 3D velocities (m s-1) from 3D transports (m2 s-1) by dividing by
        # layer heights. Both velocities and U/V thicknesses are now at time 1/2.
        np.divide(
            self.pk.all_values, self.pk.grid.hn.all_values, out=self.uk.all_values
        )
        np.divide(
            self.qk.all_values, self.qk.grid.hn.all_values, out=self.vk.all_values
        )

        # Use updated velocities (uk, vk) to compute shear frequency (SS) at T points
        # (interior only, not in halos)
        self.shear_frequency(self.uk, self.vk, viscosity, self.SS)

        # Calculate bottom friction from updated velocities (and synchronized layer
        # thicknesses hn). This needs to be done before derived quantities such as
        # bottom stress are calculated.
        if self.apply_bottom_friction:
            self.bottom_friction(
                self.u_bot, self.v_bot, self.hU_bot, self.hV_bot, self.rru, self.rrv
            )

        pygetm._pygetm.bottom_shear_velocity(
            self.u_bot,
            self.v_bot,
            self.rru,
            self.rrv,
            self.ustar2_bx,
            self.ustar2_by,
            self.ustar_b,
        )
        if self.taub.saved:
            # compute total bottom stress in Pa for e.g. FABM
            self.taub.all_values = self.ustar_b.all_values**2 * RHO0

        # Advect 3D u and v velocity from time=1/2 to 1 1/2 using velocities
        # interpolated to its own advection grids. Store the resulting trend, which
        # will be applied as part of the momentum update in the next timestep.
        # They will also be used to calculate the slow advection contribution to
        # depth-integrated momentum equations.
        # JB the alternative would be to interpolate transports and then divide by
        # (colocated) layer heights, like we do for 2D

        # Advection of 3D u velocity (uk)
        self.advpk.all_values = self.uk.all_values
        self.uadv.apply_3d(
            self.uk.interp(self.uua3d),
            self.vk.interp(self.uva3d),
            self.ww.interp(self.uk.grid),
            timestep,
            self.advpk,  # the "tracer" (velocity u) being advected; updated in-place
            new_h=True,
            skip_initial_halo_exchange=True,
        )
        # Reconstruct the change in transport pk from updated velocity (in advpk),
        # updated layer thickness (uadv.h), and old transport pk.
        # The result (m2 s-2) goes into advpk.
        pygetm._pygetm.reconstruct_transport_change(
            self.advpk, self.uadv.h, self.pk, timestep
        )

        # Advection of 3D v velocity (vk)
        self.advqk.all_values = self.vk.all_values
        self.vadv.apply_3d(
            self.uk.interp(self.vua3d),
            self.vk.interp(self.vva3d),
            self.ww.interp(self.vk.grid),
            timestep,
            self.advqk,  # the "tracer" (velocity v) being advected; updated in-place
            new_h=True,
            skip_initial_halo_exchange=True,
        )
        # Reconstruct the change in transport qk from updated velocity (in advqk),
        # updated layer thickness (vadv.h), and old transport qk.
        # The result (m2 s-2) goes into advqk.
        pygetm._pygetm.reconstruct_transport_change(
            self.advqk, self.vadv.h, self.qk, timestep
        )

        if self.diffuse_momentum:
            # Calculate the momentum trends (diffpk, diffqk) associated with diffusion
            # of 3D u and v velocity between time=1/2 to 1 1/2. Note that thicknesses
            # should be in sync with velocities uk and vk. This means they should lag
            # 1/2 a timestep behind the T grid (already the case for X, but for T we
            # use 1/2(ho+hn))
            pygetm._pygetm.momentum_diffusion(
                self.uua.grid.hn,
                self.uva.grid.hn,
                self.vua.grid.hn,
                self.vva.grid.hn,
                self.uk,
                self.vk,
                self._Am_const,
                self.diffpk,
                self.diffqk,
            )

        # Compute slow (3D) advection and diffusion contribution to to the
        # depth-integrated momentum equations. This is done by comparing the
        # depth-integrated 3D transport calculated above (between centers of the
        # current and next macrotime step) with the newly calculated depth-integrated
        # transport based on accumulated 2D transports (accumulated over the
        # current macrotimestep, and thus representative for its center).
        self._transport_2d_momentum(
            self.Ui, self.Vi, timestep, self.SxA, self.SyA, self.SxD, self.SyD, False
        )
        self.SxA.all_values = self.advpk.all_values.sum(axis=0) - self.SxA.all_values
        self.SyA.all_values = self.advqk.all_values.sum(axis=0) - self.SyA.all_values
        self.SxD.all_values = self.diffpk.all_values.sum(axis=0) - self.SxD.all_values
        self.SyD.all_values = self.diffqk.all_values.sum(axis=0) - self.SyD.all_values

        if self.apply_bottom_friction:
            # Note: ru and rv have been updated by _transport_2d_momentum, using
            # accumulated transports Ui and Vi (representative for t=1/2, just like uk,
            # vk, rru, rrv). The associated tendencies of depth-integrated transport
            # are -ru*u1 and -rv*v1. Slow bottom friction (stress in Pa divided by
            # density) is derived by taking the difference between the tendency of 3D
            # (bottom) transport and the inferred depth-integrated tendencies.
            self.SxF.all_values = (
                self.ustar2_bx.all_values + self.ru.all_values * self.u1.all_values
            )
            self.SyF.all_values = (
                self.ustar2_by.all_values + self.rv.all_values * self.v1.all_values
            )

    def _transport_2d_momentum(
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
            pygetm._pygetm.momentum_diffusion(
                self.uua.grid.D,
                self.uva.grid.D,
                self.vua.grid.D,
                self.vva.grid.D,
                self.u1,
                self.v1,
                self._Am_const,
                diffU,
                diffV,
            )

        # Calculate bottom friction (ru and rv) using updated depth-averaged velocities
        # u1 and v1. Warning: this uses velocities u1 and v1 at masked points, which
        # therefore need to be kept at 0
        if self.apply_bottom_friction:
            self.bottom_friction(
                self.u1,
                self.v1,
                self.u1.grid.D,
                self.v1.grid.D,
                self.ru,
                self.rv,
                update_z0b,
            )

        # Advection of u velocity (u1)
        U.interp(self.uua)
        V.interp(self.uva)
        self.uua.all_values /= self.uua.grid.D.all_values
        self.uva.all_values /= self.uva.grid.D.all_values
        advU.all_values = self.u1.all_values
        self.uadv(self.uua, self.uva, timestep, advU, skip_initial_halo_exchange=True)
        pygetm._pygetm.reconstruct_transport_change(advU, self.uadv.D, U, timestep)

        # Advection of v velocity (v1)
        U.interp(self.vua)
        V.interp(self.vva)
        self.vua.all_values /= self.vua.grid.D.all_values
        self.vva.all_values /= self.vva.grid.D.all_values
        advV.all_values = self.v1.all_values
        self.vadv(self.vua, self.vva, timestep, advV, skip_initial_halo_exchange=True)
        pygetm._pygetm.reconstruct_transport_change(advV, self.vadv.D, V, timestep)

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
        u_V = u.interp(v.grid._work)
        v_U = v.interp(u.grid._work)
        pygetm._pygetm.bottom_friction(u, v_U, DU, self.avmmol, ru, update_z0b)
        pygetm._pygetm.bottom_friction(u_V, v, DV, self.avmmol, rv, update_z0b)

    def coriolis(self, U: core.Array, out: core.Array, x_direction: bool):
        """Calculate change in transport due to Coriolis force

        Args:
            U: 2d or 3d transport in x- or y-direction (m2 s-1)
            out: array to store change in complementary transport (y or x-direction)
                due to Coriolis force (m2 s-2)
        """
        assert U.grid is not out.grid
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
        if x_direction:
            np.negative(out.all_values, out=out.all_values)

    def shear_frequency(
        self, uk: core.Array, vk: core.Array, viscosity: core.Array, out: core.Array
    ):
        """Calculate squared shear frequency

        Args:
            uk: depth-explicit velocity in x-direction (m s-1)
            vk: depth-explicit velocity in y-direction (m s-1)
            viscosity: vertical viscosity (m2 s-1)
            out: array to store squared shear frequency (s-2)
        """
        pygetm._pygetm.shear_frequency(uk, vk, out)
        # pygetm._pygetm.shear_frequency2(uk, vk, viscosity, out)

    def advance_3d_transport(
        self,
        tp3d: core.Array,
        tp2d: core.Array,
        vel3d: core.Array,
        dp: core.Array,
        cor: core.Array,
        adv: core.Array,
        diff: core.Array,
        idp: core.Array,
        taus: core.Array,
        rr: core.Array,
        viscosity: core.Array,
        dev: core.Array,
        timestep: float,
    ):
        """Advance 3D transports (layer-integrated velocities in m2 s-1) in x-
        or y-direction, using all relevant tendencies supplied as arguments.
        This also updates the halos.

        Note that 3D velocities (transports divided by thicknesses) are not
        updated here and therefore will be out of sync with transports after
        this function returns.
        """
        grid = tp3d.grid
        pygetm._pygetm.collect_3d_momentum_sources(
            dp, cor, adv, diff, idp, taus, rr, timestep, self.ea2, self.ea4
        )
        self._vertical_diffusion(
            viscosity.interp(grid),
            timestep,
            vel3d,
            molecular=self.avmmol,
            ea2=self.ea2,
            ea4=self.ea4,
            use_ho=True,
        )
        np.multiply(vel3d.all_values, grid.hn.all_values, out=tp3d.all_values)

        # Shift 3D transports to ensure their depth integral equals the depth
        # integrated velocities calculated by the 2D update. This is done by
        # adding the same offset in velocity throughout the water column.
        tp3d.all_values.sum(axis=0, out=dev.all_values)
        dev.all_values -= tp2d.all_values
        dev.all_values /= grid.D.all_values
        tp3d.all_values -= dev.all_values * grid.hn.all_values

        # Mirror into/across the open boundaries
        tp3d.mirror()

        tp3d.update_halos()
