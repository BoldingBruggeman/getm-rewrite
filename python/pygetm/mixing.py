from typing import Optional
import itertools
import os.path

import numpy as np

from . import core
from . import domain
from .open_boundaries import ArrayOpenBoundaries
from . import _pygotm
from .constants import INTERFACES, FILL_VALUE, ZERO_GRADIENT


class Turbulence:
    """Base class that provides the turbulent viscosity :attr:`num` and diffusivity
    :attr:`nuh`. When using this class directly, viscosity and diffusivity are
    prescribed, not calculated. In this case, both default to zero; assign to
    :attr:`num`/:attr:`nuh` or call ``num.set``/``nuh.set`` to change this.
    """

    def initialize(self, grid: domain.Grid):
        self.grid = grid
        self.logger = grid.domain.root_logger.getChild(self.__class__.__name__)
        self.nuh = grid.array(
            z=INTERFACES,
            name="nuh",
            units="m2 s-1",
            long_name="turbulent diffusivity of heat",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True},
        )
        self.num = grid.array(
            z=INTERFACES,
            name="num",
            units="m2 s-1",
            long_name="turbulent diffusivity of momentum",
            fill_value=FILL_VALUE,
            attrs={"_part_of_state": True},
        )
        self.nuh.fill(0.0)
        self.num.fill(0.0)

    def advance(
        self,
        timestep: float,
        ustar_s: core.Array,
        ustar_b: core.Array,
        z0s: core.Array,
        z0b: core.Array,
        NN: core.Array,
        SS: core.Array,
    ):
        pass


class GOTM(Turbulence):
    """Calculate the turbulent viscosity :attr:`num` and diffusivity :attr:`nuh`
    using the `General Ocean Turbulence Model (GOTM) <https://gotm.net>`_.
    """

    def __init__(self, path: Optional[str] = None):
        super().__init__()
        if path and not os.path.isfile(path):
            raise Exception(f"Configuration file {path} does not exist")
        self.path = path

    def initialize(self, grid: domain.Grid):
        super().initialize(grid)

        has_yaml = self.path and self.path.endswith(".yaml")
        nml_path = b"" if not self.path or has_yaml else self.path.encode("ascii")
        yaml_path = b"" if not has_yaml else self.path.encode("ascii")
        self.mix = _pygotm.Mixing(grid.nz, nml_path, yaml_path)
        self.tke = grid.array(
            fill_value=FILL_VALUE,
            z=INTERFACES,
            name="tke",
            units="m2 s-2",
            long_name="turbulent kinetic energy",
            attrs=dict(
                _part_of_state=True,
                standard_name="specific_turbulent_kinetic_energy_of_sea_water",
            ),
        )
        self.eps = grid.array(
            fill_value=FILL_VALUE,
            z=INTERFACES,
            name="eps",
            units="m2 s-3",
            long_name="energy dissipation rate",
            attrs=dict(
                _part_of_state=True,
                standard_name="specific_turbulent_kinetic_energy_dissipation_in_sea_water",
            ),
        )
        self.L = grid.array(
            fill_value=FILL_VALUE,
            z=INTERFACES,
            name="L",
            units="m",
            long_name="turbulence length scale",
            attrs=dict(
                _part_of_state=True,
                standard_name="turbulent_mixing_length_of_sea_water",
            ),
        )
        self.tke.fill(self.mix.tke[:, np.newaxis, np.newaxis])
        self.eps.fill(self.mix.eps[:, np.newaxis, np.newaxis])
        self.L.fill(self.mix.L[:, np.newaxis, np.newaxis])
        self.num.fill(self.mix.num[:, np.newaxis, np.newaxis])
        self.nuh.fill(self.mix.nuh[:, np.newaxis, np.newaxis])
        self._log()

        self.num.open_boundaries = ArrayOpenBoundaries(self.num, type=ZERO_GRADIENT)

    def _log(self):
        """Copy lines written by GOTM to stdout/stderr to the log
        """
        for line in itertools.chain(_pygotm.stdout, _pygotm.stderr):
            line = line.rstrip()
            if line:
                self.logger.info(line)

    def advance(
        self,
        timestep: float,
        ustar_s: core.Array,
        ustar_b: core.Array,
        z0s: core.Array,
        z0b: core.Array,
        NN: core.Array,
        SS: core.Array,
    ):
        """Update turbulent quantities and calculate turbulent diffusivity ``nuh`` and
        turbulent viscosity ``num``

        Args:
            timestep: time step (s)
            ustar_s: surface friction velocity (m s-1)
            ustar_b: bottom friction velocity (m s-1)
            z0s: hydrodynamic surface roughness (m)
            z0b: hydrodynamic bottom roughness (m)
            NN: squared buoyancy frequency (s-2)
            SS: squared shear frequency (s-2)
        """
        assert ustar_s.grid is self.grid and ustar_s.ndim == 2
        assert ustar_b.grid is self.grid and ustar_b.ndim == 2
        assert z0s.grid is self.grid and z0s.ndim == 2
        assert z0b.grid is self.grid and z0b.ndim == 2
        assert NN.grid is self.grid and NN.z == INTERFACES
        assert SS.grid is self.grid and SS.z == INTERFACES

        nz, ny, nx = self.grid.hn.all_values.shape
        self.mix.turbulence_3d(
            nx,
            ny,
            nz,
            2,
            nx - 2,
            2,
            ny - 2,
            timestep,
            self.grid.mask.all_values,
            self.grid.hn.all_values,
            self.grid.D.all_values,
            ustar_s.all_values,
            ustar_b.all_values,
            z0s.all_values,
            z0b.all_values,
            NN.all_values,
            SS.all_values,
            self.tke.all_values,
            self.eps.all_values,
            self.L.all_values,
            self.num.all_values,
            self.nuh.all_values,
        )

        # Take viscosity at open boundary from nearest interior point
        # Viscosity (at T points) needs to be valid at open boundary points to so it
        # can be interpolated to inward-adjacent U/V points. However, it cannot be
        # computed as the shear frequency SS is not available at the boundary.
        self.num.open_boundaries.update()

        self._log()
