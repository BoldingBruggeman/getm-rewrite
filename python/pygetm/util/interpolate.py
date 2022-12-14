import numpy as np
from numpy.typing import ArrayLike


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
        first = ind.min(axis=axis, where=valid, initial=fp.shape[axis], keepdims=True)
        last = ind.max(axis=axis, where=valid, initial=0, keepdims=True)
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
