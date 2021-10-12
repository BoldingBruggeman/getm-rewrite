import sys
import numpy

class Linear2DGridInterpolator:
    def __init__(self, x, y, xp, yp):
        xp = numpy.array(xp)
        yp = numpy.array(yp)
        x = numpy.array(x)
        y = numpy.array(y)
        assert xp.ndim == 1, 'source x coordinate must be 1D but has shape %s' % (xp.shape,)
        assert yp.ndim == 1, 'source y coordinate must be 1D but has shape %s' % (yp.shape,)
        self.nxp, self.nyp = xp.size, yp.size
        x, y = numpy.broadcast_arrays(x, y)
        dxp = numpy.diff(xp)
        dyp = numpy.diff(yp)
        assert (dxp > 0).all() or (dxp < 0).all(), 'source x coordinate must be monotonically increasing or decreasing'
        assert (dyp > 0).all() or (dyp < 0).all(), 'source y coordinate must be monotonically increasing or decreasing'
        if dxp[0] < 0:
            # reversed source x
            xp = xp[::-1]
        if dyp[0] < 0:
            # reversed source y
            yp = yp[::-1]
        assert (x >= xp[0]).all() and (x <= xp[-1]).all(), 'One or more target x coordinates (%s - %s) fall outside of source range (%s - %s)' % (x.min(), x.max(), xp[0], xp[-1])
        assert (y >= yp[0]).all() and (y <= yp[-1]).all(), 'One or more target y coordinates (%s - %s) fall outside of source range (%s - %s)' % (y.min(), y.max(), yp[0], yp[-1])
        self.ix_right = numpy.minimum(xp.searchsorted(x, side='right'), xp.size - 1)
        self.ix_left = self.ix_right - 1
        self.iy_right = numpy.minimum(yp.searchsorted(y, side='right'), yp.size - 1)
        self.iy_left = self.iy_right - 1
        wx_left = (xp[self.ix_right] - x) / (xp[self.ix_right] - xp[self.ix_left])
        wy_left = (yp[self.iy_right] - y) / (yp[self.iy_right] - yp[self.iy_left])
        self.w11 = wx_left * wy_left
        self.w12 = wx_left * (1. - wy_left)
        self.w21 = (1. - wx_left) * wy_left
        self.w22 = (1. - wx_left) * (1. - wy_left)
        if dxp[0] < 0:
            self.ix_left, self.ix_right = xp.size - self.ix_left - 1, xp.size - self.ix_right - 1
        if dyp[0] < 0:
            self.iy_left, self.iy_right = yp.size - self.iy_left - 1, yp.size - self.iy_right - 1

    def __call__(self, fp):
        assert fp.ndim == 2
        assert fp.shape[0] == self.nxp
        assert fp.shape[1] == self.nyp
        f11 = fp[self.ix_left, self.iy_left]
        f12 = fp[self.ix_left, self.iy_right]
        f21 = fp[self.ix_right, self.iy_left]
        f22 = fp[self.ix_right, self.iy_right]
        return self.w11 * f11 + self.w12 * f12 + self.w21 * f21 + self.w22 * f22

def test():
    # Generate random x and y source grid (1D) and a corresponding random f (2D)
    # Draw random x and y from the source space, do linear interpolation, and compare with authoratitive scipy result.
    import numpy.random
    import scipy.interpolate
    eps = 5 * numpy.finfo(float).eps
    dx = numpy.random.random_sample((99,))
    dy = numpy.random.random_sample((100,))
    xp = numpy.cumsum(dx)
    yp = numpy.cumsum(dy)
    shape = (50, 51)
    x = numpy.random.uniform(xp[0], xp[-1], shape)
    y = numpy.random.uniform(yp[0], yp[-1], shape)

    def compare(a, b):
        assert a == b, '%s vs %s' % (a, b)

    fp = numpy.random.random_sample((dx.size, dy.size,))
    compare(Linear2DGridInterpolator(xp[0], yp[0], xp, yp)(fp), fp[0, 0])
    compare(Linear2DGridInterpolator(xp[0], yp[-1], xp, yp)(fp), fp[0, -1])
    compare(Linear2DGridInterpolator(xp[-1], yp[0], xp, yp)(fp), fp[-1, 0])
    compare(Linear2DGridInterpolator(xp[-1], yp[-1], xp, yp)(fp), fp[-1, -1])

    ip = Linear2DGridInterpolator(x, y, xp, yp)
    f = ip(fp)
    f_check = scipy.interpolate.interpn((xp, yp), fp, (x, y))
    print('Maximum relative error: %s' % numpy.abs(f / f_check - 1).max())
    return numpy.isclose(f, f_check, rtol=eps, atol=eps).all()

if __name__ == '__main__':
    for _ in range(100):
        if not test():
            print('Large error found')
            sys.exit(1)
