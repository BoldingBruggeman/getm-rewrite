import logging

import numpy
import numpy.typing

def check_zero(name: str, values: numpy.typing.ArrayLike, tol: float=1e-15, logger: logging.Logger=None):
    logger = logger or logging.getLogger()
    max_val = numpy.abs(values).max()
    success = max_val < tol
    logger.log(logging.INFO if success else logging.ERROR, '%s: %.6e' % (name, max_val))
    return success