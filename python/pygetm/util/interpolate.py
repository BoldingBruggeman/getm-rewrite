import sys
import numpy
import argparse

class Linear2DGridInterpolator:
    def __init__(self, x, y, xp, yp, preslice=(Ellipsis,), ndim_trailing: int=0, mask=None):
        assert ndim_trailing >= 0
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
        ix_right = numpy.minimum(xp.searchsorted(x, side='right'), xp.size - 1)
        ix_left = ix_right - 1
        iy_right = numpy.minimum(yp.searchsorted(y, side='right'), yp.size - 1)
        iy_left = iy_right - 1
        wx_left = (xp[ix_right] - x) / (xp[ix_right] - xp[ix_left])
        wy_left = (yp[iy_right] - y) / (yp[iy_right] - yp[iy_left])
        self.w11 = wx_left * wy_left
        self.w12 = wx_left * (1. - wy_left)
        self.w21 = (1. - wx_left) * wy_left
        self.w22 = (1. - wx_left) * (1. - wy_left)

        # Ensure weights are broadcastable to shape of data array
        wshape = x.shape + (1,) * ndim_trailing
        self.w11 = numpy.reshape(self.w11, wshape)
        self.w12 = numpy.reshape(self.w12, wshape)
        self.w21 = numpy.reshape(self.w21, wshape)
        self.w22 = numpy.reshape(self.w22, wshape)

        # If we reversed source coordinates, compute the correct indices and swap left and right
        if dxp[0] < 0:
            ix_left, ix_right = xp.size - ix_left - 1, xp.size - ix_right - 1
        if dyp[0] < 0:
            iy_left, iy_right = yp.size - iy_left - 1, yp.size - iy_right - 1

        # Store slices into data array
        self.slice11 = (Ellipsis, ix_left, iy_left) + tuple([slice(None)] * ndim_trailing)
        self.slice12 = (Ellipsis, ix_left, iy_right) + tuple([slice(None)] * ndim_trailing)
        self.slice21 = (Ellipsis, ix_right, iy_left) + tuple([slice(None)] * ndim_trailing)
        self.slice22 = (Ellipsis, ix_right, iy_right) + tuple([slice(None)] * ndim_trailing)

        if mask is not None:
            # Force weights to zero for masked points and renormalize the weights so their sum is 1
            shape = mask.shape[:-ndim_trailing - 2] + x.shape + mask.shape[-ndim_trailing:]
            w11_full = numpy.empty(shape, dtype=self.w11.dtype)
            w12_full = numpy.empty(shape, dtype=self.w12.dtype)
            w21_full = numpy.empty(shape, dtype=self.w21.dtype)
            w22_full = numpy.empty(shape, dtype=self.w22.dtype)
            w11_full[...] = self.w11
            w12_full[...] = self.w12
            w21_full[...] = self.w21
            w22_full[...] = self.w22
            w11_full[mask[self.slice11]] = 0.
            w12_full[mask[self.slice12]] = 0.
            w21_full[mask[self.slice21]] = 0.
            w22_full[mask[self.slice22]] = 0.
            norm = 1. / (w11_full + w12_full + w21_full + w22_full)
            self.w11 = w11_full * norm
            self.w12 = w12_full * norm
            self.w21 = w21_full * norm
            self.w22 = w22_full * norm

        self.idim1 = -2 - ndim_trailing
        self.idim2 = -1 - ndim_trailing

    def __call__(self, fp):
        assert fp.shape[self.idim1] == self.nxp
        assert fp.shape[self.idim2] == self.nyp
        f11 = fp[self.slice11]
        f12 = fp[self.slice12]
        f21 = fp[self.slice21]
        f22 = fp[self.slice22]
        return self.w11 * f11 + self.w12 * f12 + self.w21 * f21 + self.w22 * f22

def test():
    success = True

    # Generate random x and y source grid (1D) and a corresponding random f (2D)
    # Draw random x and y from the source space, do linear interpolation, and compare with authoratitive scipy result.
    import numpy.random
    import scipy.interpolate

    eps = 5 * numpy.finfo(float).eps

    # Random source grid (1D x and y) - but still monotonically increasing because of the use of cumsum
    dx = numpy.random.random_sample((99,))
    dy = numpy.random.random_sample((100,))
    xp = numpy.cumsum(dx)
    yp = numpy.cumsum(dy)

    # Random source field
    fp = numpy.random.random_sample((xp.size, yp.size,))

    # Interpolate to source corner positions and verify results are identical to source values
    def compare(name: str, a, b) -> bool:
        if a != b:
            print('Interpolated value does not match source (%s vs. %s) at %s.' % (a, b, name))
        return a == b
    success = compare('min x, min y', Linear2DGridInterpolator(xp[0], yp[0], xp, yp)(fp), fp[0, 0]) and success
    success = compare('min x, max y', Linear2DGridInterpolator(xp[0], yp[-1], xp, yp)(fp), fp[0, -1]) and success
    success = compare('max x, min y', Linear2DGridInterpolator(xp[-1], yp[0], xp, yp)(fp), fp[-1, 0]) and success
    success = compare('max x, max y', Linear2DGridInterpolator(xp[-1], yp[-1], xp, yp)(fp), fp[-1, -1]) and success

    f_check = scipy.interpolate.interpn((xp, yp), fp, (xp[:, numpy.newaxis], yp[numpy.newaxis, :]))
    print('  Interpolate to original grid: Maximum relative error: %s' % (numpy.abs(fp / f_check - 1).max(),))
    success = (fp - f_check == 0).all() and success
    
    # Interpolate to 2D target grid and comapre resuts with that of scipy.interpolate.interpn
    shape = (50, 51)
    x = numpy.random.uniform(xp[0], xp[-1], shape)
    y = numpy.random.uniform(yp[0], yp[-1], shape)
    ip = Linear2DGridInterpolator(x, y, xp, yp)
    f = ip(fp)
    f_check = scipy.interpolate.interpn((xp, yp), fp, (x, y))
    print('  2D: Maximum relative error: %s' % (numpy.abs(f / f_check - 1).max(),))
    success = numpy.isclose(f, f_check, rtol=eps, atol=eps).all() and success

    # Interpolate to 1D target grid and comapre resuts with that of scipy.interpolate.interpn
    shape = (100,)
    x = numpy.random.uniform(xp[0], xp[-1], shape)
    y = numpy.random.uniform(yp[0], yp[-1], shape)
    ip = Linear2DGridInterpolator(x, y, xp, yp)
    f = ip(fp)
    f_check = scipy.interpolate.interpn((xp, yp), fp, (x, y))
    print('  1D: Maximum relative error: %s' % (numpy.abs(f / f_check - 1).max(),))
    success = numpy.isclose(f, f_check, rtol=eps, atol=eps).all() and success

    return success

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('N', type=int, help='number of tests')
    args = parser.parse_args()
    for itest in range(args.N):
        print('Test %i: ' % itest)
        if not test():
            print('Large error found')
            sys.exit(1)
