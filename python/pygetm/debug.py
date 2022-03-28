import logging

from . import core
import numpy
import numpy.typing

def check_zero(name: str, values: numpy.typing.ArrayLike, tol: float=1e-15, logger: logging.Logger=None) -> bool:
    logger = logger or logging.getLogger()
    max_val = numpy.abs(values).max()
    success = max_val < tol
    logger.log(logging.INFO if success else logging.ERROR, '%s: %.6e' % (name, max_val))
    return success

def check_equal(name: str, values: numpy.typing.ArrayLike, values_ref: numpy.typing.ArrayLike, rtol=1e-15, atol=1e-15, logger: logging.Logger=None) -> bool:
    logger = logger or logging.getLogger()
    values = numpy.asarray(values)
    values_ref = numpy.asarray(values_ref)
    abs_dif = abs(values - values_ref)
    valid = abs_dif < rtol * abs(values_ref) + atol
    success = valid.all()
    logger.log(logging.INFO if success else logging.ERROR, '%s: maximum absolute difference %.6e' % (name, abs_dif.max()))
    return success

def check_range(name: str, values: numpy.typing.ArrayLike, rtol=1e-12, atol=1e-12, target_value=None, logger: logging.Logger=None, plot: bool=False) -> bool:
    logger = logger or logging.getLogger()
    values = numpy.asarray(values)
    absrange = values.max() - values.min()
    mean = values.mean()
    crit = rtol * abs(mean) + atol
    valid = absrange < crit
    result = 'OK' if valid else 'FAILED because field is not homogeneous'
    relrange = 0 if absrange == 0 else absrange / abs(mean)
    if valid and target_value is not None:
        valid = abs(mean - target_value) < rtol * abs(target_value) + atol
        if not valid:
            result = 'FAILED to match target value %.6g' % target_value
    logger.log(logging.INFO if valid else logging.ERROR, '%s %s (mean %.6g, spread %.6g, relative spread %.6g)' % (name, result, mean, absrange, relrange))
    if not valid and plot:
        plot_2d('%s.png' % name, values)
    return valid

def check_finite(arr: core.Array, logger: logging.Logger=None) -> bool:
    finite = numpy.isfinite(arr.values)
    if not finite.all():
        logger = logger or logging.getLogger()
        logger.error('%s not finite in subdomain %i,%i' % (arr.name,), arr.grid.domain.tiling.irow, arr.grid.domain.tiling.icol)
        return False
    return True

def log_range(arr: core.Array, logger: logging.Logger=None):
    logger = logger or logging.getLogger()
    minval = arr.ma.min()
    maxval = arr.ma.max()
    logger.info('%s: %s - %s in subdomain %i,%i' % (arr.name, minval, maxval, arr.grid.domain.tiling.irow, arr.grid.domain.tiling.icol))

def plot_2d(path: str, values: numpy.ndarray):
    import matplotlib.pyplot
    fig, ax = matplotlib.pyplot.subplots()
    pc = ax.pcolormesh(values)
    fig.colorbar(pc)
    fig.savefig(path, dpi=300)
