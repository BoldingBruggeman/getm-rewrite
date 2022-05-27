import enum

from ._pyairsea import *


class HumidityMeasure(enum.IntEnum):
    """Measure used to specify air humidity"""

    RELATIVE_HUMIDITY = 1  #: relative humidity in %
    WET_BULB_TEMPERATURE = 2  #: wet bulb temperature in in degrees Celsius
    DEW_POINT_TEMPERATURE = 3  #: dewpoint temperature in in degrees Celsius
    SPECIFIC_HUMIDITY = 4  #: specific humidity in kg kg-1


class LongwaveMethod(enum.IntEnum):
    """Method used to calculate longwave radiation"""

    CLARK = 1  #: Clark et al. (1974)
    HASTENRATH_LAMB = 2  #: Hastenrath and Lamb (1978)
    BIGNAMI = 3  #: Bignami et al. (1995)
    BERLIAND_BERLIAND = 4  #: Berliand & Berliand (1952)
    JOSEY1 = 5  #: Josey et.al. 2003 - (J1,9)
    JOSEY2 = 6  #: Josey et.al. 2003 - (J2,14)


class AlbedoMethod(enum.IntEnum):
    """Method used to calculate albedo"""

    PAYNE = 1  #: Payne (1972)
    COGLEY = 2  #: Cogley (1979)
