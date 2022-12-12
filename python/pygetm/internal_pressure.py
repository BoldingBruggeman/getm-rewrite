from . import core
from . import domain
from . import _pygetm
from .constants import FILL_VALUE, CENTERS


class Base:
    idpdx: core.Array
    idpdy: core.Array

    def __init__(self, idpdx_fill: float = 0.0, idpdy_fill: float = 0.0):
        self.idpdx_fill = idpdx_fill
        self.idpdy_fill = idpdy_fill

    def initialize(self, domain: domain.Domain):
        self.idpdx = domain.U.array(
            name="idpdx",
            units="m2 s-2",
            long_name="internal pressure gradient in x-direction",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.idpdy = domain.V.array(
            name="idpdy",
            units="m2 s-2",
            long_name="internal pressure gradient in y-direction",
            z=CENTERS,
            fill_value=FILL_VALUE,
        )
        self.idpdx.fill(self.idpdx_fill)
        self.idpdy.fill(self.idpdy_fill)

    def __call__(self, buoy: core.Array):
        raise NotImplementedError


class Constant(Base):
    def __call__(self, buoy: core.Array):
        return


class BlumbergMellor(Base):
    def __call__(self, buoy: core.Array):
        _pygetm.blumberg_mellor(buoy, self.idpdx, self.idpdy)


class ShchepetkinMcwilliams(Base):
    def __call__(self, buoy: core.Array):
        _pygetm.shchepetkin_mcwilliams(buoy, self.idpdx, self.idpdy)
