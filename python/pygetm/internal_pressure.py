from . import core
from . import _pygetm
from .constants import FILL_VALUE, CENTERS


class Base:
    idpdx: core.Array
    idpdy: core.Array

    def initialize(self, ugrid: core.Grid, vgrid: core.Grid):
        self.idpdx = ugrid.array(
            name="idpdx",
            units="m2 s-2",
            long_name="tendency of layer-integrated x-velocity due to internal pressure",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.idpdy = vgrid.array(
            name="idpdy",
            units="m2 s-2",
            long_name="tendency of layer-integrated y-velocity due to internal pressure",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )

    def __call__(self, buoy: core.Array):
        """Update tendencies of layer-integrated velocity due to internal pressure.

        These are defined as -hu/rho0 * dp/dx and -hv/rho0 * dp/dy,
        with all quantities defined at the current (tracer) time level.

        Args:
            buoy: buoyancy (m s-2)
        """
        raise NotImplementedError


class Prescribed(Base):
    """Prescribed internal pressure gradient.

    To set the gradients, provide arguments `idpdx` and `idpdy` when creating
    an instance, or call :meth:`pygetm.core.Array.set` on the `idpdx` and `idpdy`
    attributes of the instance after initialization."""

    def __init__(self, idpdx: float = 0.0, idpdy: float = 0.0):
        """
        Args:
            idpdx: initial tendency of layer-integrated x-velocity (m2 s-2)
            idpdy: initial tendency of layer-integrated y-velocity (m2 s-2)
        """
        super().__init__()
        self._initial_idpdx = idpdx
        self._initial_idpdy = idpdy

    def initialize(self, ugrid: core.Grid, vgrid: core.Grid):
        super().initialize(ugrid, vgrid)
        self.idpdx.fill(self._initial_idpdx)
        self.idpdy.fill(self._initial_idpdy)

    def __call__(self, buoy: core.Array):
        return


class BlumbergMellor(Base):
    """Internal pressure gradient calculated as described by `Blumberg and Mellor (1987)
    <https://doi.org/10.1029/CO004p0001>`_."""

    def __call__(self, buoy: core.Array):
        _pygetm.blumberg_mellor(buoy, self.idpdx, self.idpdy)


class ShchepetkinMcwilliams(Base):
    """Internal pressure gradient calculated as described by `Shchepetkin and McWilliams
    (2003) <https://doi.org/10.1029/2001JC001047>`_."""

    def __call__(self, buoy: core.Array):
        _pygetm.shchepetkin_mcwilliams(buoy, self.idpdx, self.idpdy)
