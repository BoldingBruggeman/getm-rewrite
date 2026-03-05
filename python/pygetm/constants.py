import enum

RHO0 = 1025.0
GRAVITY = 9.81  #: gravitational acceleration

CENTERS = 1
INTERFACES = 2

ZERO_GRADIENT = 1
SOMMERFELD = 2
CLAMPED = 3
FLATHER_ELEV = 4
FLATHER_TRANSPORT = -4
SPONGE = 5

FILL_VALUE = -2e20


class RunType(enum.IntEnum):
    BAROTROPIC_2D = 1
    BAROTROPIC_3D = 2
    FROZEN_DENSITY = 3
    BAROCLINIC = 4


BAROTROPIC = BAROTROPIC_2D = RunType.BAROTROPIC_2D
BAROTROPIC_3D = RunType.BAROTROPIC_3D
FROZEN_DENSITY = RunType.FROZEN_DENSITY
BAROCLINIC = RunType.BAROCLINIC


class TimeVarying(enum.Enum):
    MACRO = enum.auto()
    MICRO = enum.auto()


class CoordinateType(enum.Enum):
    XY = enum.auto()
    LONLAT = enum.auto()
    IJ = enum.auto()


class CellType(enum.IntEnum):
    UNRESOLVED = 0  #: land point
    ACTIVE = 1  #: actively modeled wet point
    BOUNDARY = 2  #: wet point that is part of an open boundary
    MIRROR_INT = 3  #: wet point within boundary, set by mirroring the interior
    MIRROR_EXT = 4  #: wet point just outside boundary, set by mirroring the interior
    EDGE_X = (
        -1
    )  #: edge between actively modeled and unresolved points (e.g., coastlines)
    EDGE_Y = (
        -2
    )  #: edge between actively modeled and unresolved points (e.g., coastlines)


class EdgeTreatment(enum.Enum):
    MISSING = enum.auto()
    CLAMP = enum.auto()
    PERIODIC = enum.auto()
    EXTRAPOLATE = enum.auto()
    EXTRAPOLATE_PERIODIC = enum.auto()
