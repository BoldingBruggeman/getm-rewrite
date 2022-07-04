import sys
import numpy as np
from numpy.typing import ArrayLike
import argparse


class Linear2DGridInterpolator:
    def __init__(
        self,
        x: ArrayLike,
        y: ArrayLike,
        xp: ArrayLike,
        yp: ArrayLike,
        preslice=(Ellipsis,),
        ndim_trailing: int = 0,
        mask=None,
    ):
        assert ndim_trailing >= 0
        xp = np.array(xp)
        yp = np.array(yp)
        x = np.array(x)
        y = np.array(y)
        assert xp.ndim == 1, "source x coordinate must be 1D but has shape %s" % (
            xp.shape,
        )
        assert yp.ndim == 1, "source y coordinate must be 1D but has shape %s" % (
            yp.shape,
        )
        self.nxp, self.nyp = xp.size, yp.size
        assert (
            self.nxp > 1
        ), "source x coordinate must have length > 1, but has length %i" % (self.nxp,)
        assert (
            self.nyp > 1
        ), "source y coordinate must have length > 1, but has length %i" % (self.nyp,)
        x, y = np.broadcast_arrays(x, y)
        dxp = np.diff(xp)
        dyp = np.diff(yp)
        assert (dxp > 0).all() or (
            dxp < 0
        ).all(), "source x coordinate must be monotonically increasing or decreasing"
        assert (dyp > 0).all() or (
            dyp < 0
        ).all(), "source y coordinate must be monotonically increasing or decreasing"
        if dxp[0] < 0:
            # reversed source x
            xp = xp[::-1]
        if dyp[0] < 0:
            # reversed source y
            yp = yp[::-1]
        assert (x >= xp[0]).all() and (x <= xp[-1]).all(), (
            "One or more target x coordinates (%s - %s) fall outside of source range (%s - %s)"
            % (x.min(), x.max(), xp[0], xp[-1])
        )
        assert (y >= yp[0]).all() and (y <= yp[-1]).all(), (
            "One or more target y coordinates (%s - %s) fall outside of source range (%s - %s)"
            % (y.min(), y.max(), yp[0], yp[-1])
        )
        ix_right = np.minimum(xp.searchsorted(x, side="right"), xp.size - 1)
        ix_left = ix_right - 1
        iy_right = np.minimum(yp.searchsorted(y, side="right"), yp.size - 1)
        iy_left = iy_right - 1
        wx_left = (xp[ix_right] - x) / (xp[ix_right] - xp[ix_left])
        wy_left = (yp[iy_right] - y) / (yp[iy_right] - yp[iy_left])
        self.w11 = wx_left * wy_left
        self.w12 = wx_left * (1.0 - wy_left)
        self.w21 = (1.0 - wx_left) * wy_left
        self.w22 = (1.0 - wx_left) * (1.0 - wy_left)

        # Ensure weights are broadcastable to shape of data array
        wshape = x.shape + (1,) * ndim_trailing
        self.w11 = np.reshape(self.w11, wshape)
        self.w12 = np.reshape(self.w12, wshape)
        self.w21 = np.reshape(self.w21, wshape)
        self.w22 = np.reshape(self.w22, wshape)

        # If we reversed source coordinates, compute the correct indices
        if dxp[0] < 0:
            ix_left, ix_right = xp.size - ix_left - 1, xp.size - ix_right - 1
        if dyp[0] < 0:
            iy_left, iy_right = yp.size - iy_left - 1, yp.size - iy_right - 1

        # Store slices into data array
        self.slice11 = (Ellipsis, ix_left, iy_left) + (slice(None),) * ndim_trailing
        self.slice12 = (Ellipsis, ix_left, iy_right) + (slice(None),) * ndim_trailing
        self.slice21 = (Ellipsis, ix_right, iy_left) + (slice(None),) * ndim_trailing
        self.slice22 = (Ellipsis, ix_right, iy_right) + (slice(None),) * ndim_trailing

        if mask is not None:
            # Force weights to zero for masked points and renormalize the weights
            # so their sum is 1
            shape = (
                mask.shape[: -ndim_trailing - 2] + x.shape + mask.shape[-ndim_trailing:]
            )
            w11_full = np.empty(shape, dtype=self.w11.dtype)
            w12_full = np.empty(shape, dtype=self.w12.dtype)
            w21_full = np.empty(shape, dtype=self.w21.dtype)
            w22_full = np.empty(shape, dtype=self.w22.dtype)
            w11_full[...] = self.w11
            w12_full[...] = self.w12
            w21_full[...] = self.w21
            w22_full[...] = self.w22
            w11_full[mask[self.slice11]] = 0.0
            w12_full[mask[self.slice12]] = 0.0
            w21_full[mask[self.slice21]] = 0.0
            w22_full[mask[self.slice22]] = 0.0
            norm = 1.0 / (w11_full + w12_full + w21_full + w22_full)
            self.w11 = w11_full * norm
            self.w12 = w12_full * norm
            self.w21 = w21_full * norm
            self.w22 = w22_full * norm

        self.idim1 = -2 - ndim_trailing
        self.idim2 = -1 - ndim_trailing

    def __call__(self, fp: ArrayLike) -> np.ndarray:
        assert fp.shape[self.idim1] == self.nxp
        assert fp.shape[self.idim2] == self.nyp
        result = self.w11 * fp[self.slice11]
        result += self.w12 * fp[self.slice12]
        result += self.w21 * fp[self.slice21]
        result += self.w22 * fp[self.slice22]
        return result


class LinearVectorized1D:
    def __init__(self, x, xp, axis=0, fill_value=np.nan):
        x = np.asarray(x)
        xp = np.asarray(xp)
        assert x.ndim == 1
        assert axis >= -xp.ndim and axis < xp.ndim
        xp_slice = [slice(None)] * xp.ndim
        final_shape = list(xp.shape)
        final_shape[axis] = x.size
        ix_right = np.empty(final_shape, dtype=int)
        for ix, cur_x in enumerate(x):
            xp_slice[axis] = ix
            ix_right[tuple(xp_slice)] = (xp <= cur_x).sum(axis=axis)
        ix_left = np.maximum(ix_right - 1, 0)
        ix_right = np.minimum(ix_right, xp.shape[axis] - 1)
        valid = ix_left != ix_right
        xp_right = np.take_along_axis(xp, ix_right, axis=axis)
        xp_left = np.take_along_axis(xp, ix_left, axis=axis)
        dxp = xp_right - xp_left
        x_shape = [1] * xp.ndim
        x_shape[axis] = x.size
        x_bc = x.reshape(x_shape)
        w_left = np.true_divide(xp_right - x_bc, dxp, where=valid)
        np.putmask(w_left, ~valid, 1.0)
        self.ix_left = ix_left
        self.ix_right = ix_right
        self.w_left = w_left
        self.axis = axis
        self.valid = valid
        self.fill_value = fill_value

    def __call__(self, yp):
        yp = np.asarray(yp)
        yp_left = np.take_along_axis(yp, self.ix_left, axis=self.axis)
        yp_right = np.take_along_axis(yp, self.ix_right, axis=self.axis)
        y = self.w_left * yp_left + (1.0 - self.w_left) * yp_right
        y = np.where(self.valid, y, self.fill_value)
        return y


def interp_1d(x, xp, fp, axis: int = 0):
    x = np.asarray(x)
    xp = np.asarray(xp)
    fp = np.ma.filled(fp, np.nan)
    assert fp.ndim == x.ndim, (
        "Number of dimensions %i of source values does not match %i"
        " of target coordinate." % (fp.ndim, x.ndim)
    )
    assert xp.ndim == 1, "Source coordinate must be 1D but has shape %s." % (xp.shape,)
    assert (
        fp.shape[:axis] == x.shape[:axis]
        and fp.shape[axis + 1 :] == x.shape[axis + 1 :]
    ), (
        "Shapes of source values %s and target coordinate %s should match everywhere"
        " except the depth dimension (%i)" % (fp.shape, x.shape, axis)
    )
    assert fp.shape[axis] == xp.shape[0]

    dxp = np.diff(xp)
    assert (dxp > 0).all() or (
        dxp < 0
    ).all(), "source coordinate must be monotonically increasing or decreasing"
    if dxp[0] < 0:
        # reversed source coordinate
        xp = xp[::-1]

    # Deal with masked sections at beginning or end of source values
    invalid = np.isnan(fp)
    if invalid.any():
        valid = ~invalid
        ind = np.arange(fp.shape[axis])
        ind.shape = [1 if i != axis else -1 for i in range(fp.ndim)]
        ind = np.broadcast_to(ind, fp.shape)
        first = ind.min(axis=axis, where=valid, initial=fp.shape[axis])
        last = ind.max(axis=axis, where=valid, initial=0)
        first = np.minimum(first, last)  # if no valid elements at all, first=last=0
        if dxp[0] < 0:
            first, last = fp.shape[axis] - last - 1, fp.shape[axis] - first - 1
    else:
        first = 0
        last = fp.shape[axis] - 1

    # Look up upper bound of interval around each target coordinate
    # This will be 0 [invalid!] if first source coordinate < minimum target coordinate
    # This will be xp.size [invalid!] if last source coordinate >= maximum target
    # coordinate
    ix_right = xp.searchsorted(x, side="right")

    # Determine intervals left and right bounds).
    # These will be zero-width (ix_left == ix_right) at the boundaries
    ix_left = np.clip(ix_right - 1, first, last)
    np.clip(ix_right, first, last, out=ix_right)
    valid_interval = ix_left != ix_right
    xp_right = xp[ix_right]
    wx_left = np.true_divide(xp_right - x, xp_right - xp[ix_left], where=valid_interval)
    np.putmask(wx_left, ~valid_interval, 1.0)

    # If we reversed source coordinates, compute the correct indices
    if dxp[0] < 0:
        ix_left, ix_right = xp.size - ix_left - 1, xp.size - ix_right - 1

    f_left = np.take_along_axis(fp, ix_left, axis=axis)
    f_right = np.take_along_axis(fp, ix_right, axis=axis)
    return wx_left * f_left + (1.0 - wx_left) * f_right


def test():
    # Generate random x and y source grid (1D) and a corresponding random f (2D)
    # Draw random x and y from the source space, do linear interpolation,
    # and compare with authoratitive scipy result.
    import numpy.random
    import scipy.interpolate

    success = True

    # Maximum acceptable error
    eps = 5 * np.finfo(float).eps

    print("  Horizontal interpolation:")

    # Random source grid (1D x and y)
    # Still monotonically increasing because of the use of cumsum
    dx = numpy.random.random_sample((99,))
    dy = numpy.random.random_sample((100,))
    xp = np.cumsum(dx)
    yp = np.cumsum(dy)

    # Random source field
    fp = numpy.random.random_sample((xp.size, yp.size,))

    # Interpolate to source corner positions
    # and verify results are identical to source values
    def compare(name: str, a, b) -> bool:
        if a != b:
            print(
                "    ERROR interpolated value does not match source (%s vs. %s) at %s."
                % (a, b, name)
            )
        return a == b

    success = (
        compare(
            "min x, min y", Linear2DGridInterpolator(xp[0], yp[0], xp, yp)(fp), fp[0, 0]
        )
        and success
    )
    success = (
        compare(
            "min x, max y",
            Linear2DGridInterpolator(xp[0], yp[-1], xp, yp)(fp),
            fp[0, -1],
        )
        and success
    )
    success = (
        compare(
            "max x, min y",
            Linear2DGridInterpolator(xp[-1], yp[0], xp, yp)(fp),
            fp[-1, 0],
        )
        and success
    )
    success = (
        compare(
            "max x, max y",
            Linear2DGridInterpolator(xp[-1], yp[-1], xp, yp)(fp),
            fp[-1, -1],
        )
        and success
    )

    f_check = scipy.interpolate.interpn(
        (xp, yp), fp, (xp[:, np.newaxis], yp[np.newaxis, :])
    )
    print(
        "    Interpolate to original grid: maximum relative error: %s"
        % (np.abs(fp / f_check - 1).max(),)
    )
    success = (fp - f_check == 0).all() and success

    def compare_spatial(name: str, x, y):
        success = True

        f_check = scipy.interpolate.interpn((xp, yp), fp, (x, y))

        f = Linear2DGridInterpolator(x, y, xp, yp)(fp)
        print(
            "    %s: maximum relative error: %s"
            % (name, np.abs(f / f_check - 1).max(),)
        )
        success = np.isclose(f, f_check, rtol=eps, atol=eps).all() and success

        f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp)(fp[::-1, :])
        print(
            "    %s x reversed: maximum relative difference with original: %s"
            % (name, np.abs(f2 / f - 1).max(),)
        )
        success = (f - f2 == 0).all() and success

        f2 = Linear2DGridInterpolator(x, y, xp, yp[::-1])(fp[:, ::-1])
        print(
            "    %s y reversed: maximum relative difference with original: %s"
            % (name, np.abs(f2 / f - 1).max(),)
        )
        success = (f - f2 == 0).all() and success

        f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp[::-1])(fp[::-1, ::-1])
        print(
            "    %s xy reversed: maximum relative difference with original: %s"
            % (name, np.abs(f2 / f - 1).max(),)
        )
        success = (f - f2 == 0).all() and success

        return success

    # Interpolate to 2D target grid
    # and compare results with that of scipy.interpolate.interpn
    shape = (50, 51)
    x = numpy.random.uniform(xp[0], xp[-1], shape)
    y = numpy.random.uniform(yp[0], yp[-1], shape)
    success = compare_spatial("2D", x, y) and success

    # Interpolate to 1D target grid
    # and compare results with that of scipy.interpolate.interpn
    shape = (100,)
    x = numpy.random.uniform(xp[0], xp[-1], shape)
    y = numpy.random.uniform(yp[0], yp[-1], shape)
    success = compare_spatial("1D", x, y) and success

    print("  Vertical interpolation:")
    dx = numpy.random.random_sample((99,))
    xp = np.cumsum(dx)
    fp = numpy.random.random_sample((xp.size,))
    x = numpy.random.uniform(xp[0] - 5, xp[-1] + 5, (100,))
    f = interp_1d(x, xp, fp)
    f_check = np.interp(x, xp, fp)
    print("    1D: maximum relative error: %s" % (np.abs(f / f_check - 1).max(),))
    success = np.isclose(f, f_check, rtol=eps, atol=eps).all() and success

    f2 = interp_1d(x, xp[::-1], fp[::-1])
    print(
        "    1D reversed: maximum relative difference with original: %s"
        % (np.abs(f2 / f - 1).max(),)
    )
    success = (f == f2).all() and success

    nx, ny = 5, 6
    dx = numpy.random.random_sample((99,))
    xp = np.cumsum(dx)
    fp = numpy.random.random_sample((xp.size, ny, nx))
    x = -10 * numpy.random.random_sample((ny, nx)) + 1.5 * numpy.random.random_sample(
        (100, ny, nx)
    ).cumsum(axis=0)
    f = interp_1d(x, xp, fp)
    assert x.shape == f.shape
    current_success = True
    error = 0
    for i in range(x.shape[-1]):
        for j in range(x.shape[-2]):
            f_check = np.interp(x[:, j, i], xp, fp[:, j, i])
            current_success = (
                np.isclose(f[:, j, i], f_check, rtol=eps, atol=eps).all()
                and current_success
            )
            error = max(np.abs(f[:, j, i] / f_check - 1).max(), error)
    print("    3D: maximum relative error: %s" % (error,))
    success = current_success and success

    f2 = interp_1d(x, xp[::-1], fp[::-1, ...])
    print(
        "    3D reversed: maximum relative difference with original: %s"
        % (np.abs(f2 / f - 1).max(),)
    )
    success = (f == f2).all() and success

    start = numpy.random.randint(0, fp.shape[0] - 1, fp.shape[1:])
    stop = numpy.random.randint(start, fp.shape[0])
    ind = np.broadcast_to(np.arange(fp.shape[0])[:, np.newaxis, np.newaxis], fp.shape)
    mask = np.logical_or(ind < start, ind >= stop)
    masked_fp = np.ma.array(fp, mask=mask)
    f = interp_1d(x, xp, masked_fp)
    current_success = True
    error = 0
    for i in range(x.shape[-1]):
        for j in range(x.shape[-2]):
            if start[j, i] == stop[j, i]:
                # no valid points - all output should be masked
                assert np.isnan(f[:, j, i]).all(), "%s" % (f[:, j, i],)
            else:
                assert np.isfinite(fp[start[j, i] : stop[j, i], j, i]).all()
                f_check = np.interp(
                    x[:, j, i],
                    xp[start[j, i] : stop[j, i]],
                    fp[start[j, i] : stop[j, i], j, i],
                )
                current_success = (
                    np.isclose(f[:, j, i], f_check, rtol=eps, atol=eps).all()
                    and current_success
                )
                error = max(np.abs(f[:, j, i] / f_check - 1).max(), error)
                if not current_success:
                    print(
                        start[j, i],
                        stop[j, i],
                        mask[:, j, i],
                        f[:, j, i],
                        f_check,
                        f[:, j, i] - f_check,
                        error,
                    )
                    dsds
    print("    3D masked: maximum relative error: %s" % (error,))

    return success


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("N", type=int, help="number of tests")
    args = parser.parse_args()
    for itest in range(args.N):
        print("Test %i: " % itest)
        if not test():
            print("Large error found")
            sys.exit(1)
