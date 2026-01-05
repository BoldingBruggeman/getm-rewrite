from typing import Optional
import numpy as np
import numpy.typing as npt


class Linear2DGridInterpolator:
    def __init__(
        self,
        x: npt.ArrayLike,
        y: npt.ArrayLike,
        xp: npt.ArrayLike,
        yp: npt.ArrayLike,
        preslice=(Ellipsis,),
        ndim_trailing: int = 0,
        mask: Optional[npt.ArrayLike] = None,
    ):
        assert ndim_trailing >= 0
        xp = np.array(xp, dtype=float)
        yp = np.array(yp, dtype=float)
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        assert xp.ndim == 1, f"source x coordinate must be 1D but has shape {xp.shape}"
        assert yp.ndim == 1, f"source y coordinate must be 1D but has shape {yp.shape}"
        self.nxp, self.nyp = xp.size, yp.size
        assert (
            self.nxp > 1
        ), f"source x coordinate must have length > 1, but has length {self.nxp}"
        assert (
            self.nyp > 1
        ), f"source y coordinate must have length > 1, but has length {self.nyp}"
        x, y = np.broadcast_arrays(x, y)
        dxp = xp[1:] - xp[:-1]
        dyp = yp[1:] - yp[:-1]
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
            f"One or more target x coordinates ({x.min()} - {x.max()})"
            f" fall outside of source range ({xp[0]} - {xp[-1]})"
        )
        assert (y >= yp[0]).all() and (y <= yp[-1]).all(), (
            f"One or more target y coordinates ({y.min()} - {y.max()})"
            f" fall outside of source range ({yp[0]} - {yp[-1]})"
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
        assert np.allclose(self.w11 + self.w12 + self.w21 + self.w22, 1.0, 1e-14)

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
            # Force weights to zero for masked points and renormalize weights
            # so their sum is 1
            mask = np.asarray(mask)
            assert (
                mask.ndim >= 2 + ndim_trailing
            ), f"Mask should have at least {2 + ndim_trailing} dimensions"
            ndim_no_trail = mask.ndim - ndim_trailing
            mask_xy_shape = mask.shape[ndim_no_trail - 2 : ndim_no_trail]
            xy_shape = (self.nxp, self.nyp)
            assert (
                mask_xy_shape == xy_shape
            ), f"Bad mask shape for x, y: {mask_xy_shape} while expected {xy_shape}"
            target_shape = (
                mask.shape[: ndim_no_trail - 2] + x.shape + mask.shape[ndim_no_trail:]
            )
            use_nn = (
                mask[self.slice11]
                | mask[self.slice12]
                | mask[self.slice21]
                | mask[self.slice22]
            )
            if use_nn.any():
                import scipy.spatial

                # Build kd-tree for nearest-neighbor lookup, using only unmasked points
                if mask.all():
                    raise ValueError("All source points are masked, cannot interpolate")
                source_coords = np.indices(mask.shape)[:, ~mask]
                tree = scipy.spatial.KDTree(source_coords.T)

                ix = ix_right - wx_left
                iy = iy_right - wy_left
                sl = (Ellipsis,) + (np.newaxis,) * ndim_trailing
                target_indices = np.indices(target_shape)
                target_coords = []
                for i in range(len(target_shape) - 2 - ndim_trailing):
                    target_coords.append(target_indices[i])
                target_coords.append(np.broadcast_to(ix[sl], target_shape))
                target_coords.append(np.broadcast_to(iy[sl], target_shape))
                for i in range(ndim_trailing):
                    target_coords.append(target_indices[-ndim_trailing + i])
                target_coords = np.array(target_coords, dtype=float)

                _, inearest = tree.query(target_coords[:, use_nn].T, workers=-1)

                for k in ("slice11", "slice12", "slice21", "slice22"):
                    s = getattr(self, k)

                    # Broadcast original slices
                    slc_indices = []
                    for i in range(len(target_shape) - 2 - ndim_trailing):
                        slc_indices.append(target_indices[i])
                    ix = np.broadcast_to(s[-2 - ndim_trailing][sl], target_shape)
                    iy = np.broadcast_to(s[-1 - ndim_trailing][sl], target_shape)
                    slc_indices.extend([ix, iy])
                    for i in range(ndim_trailing):
                        slc_indices.append(target_indices[-ndim_trailing + i])

                    # Substitute nearest
                    for i in range(len(slc_indices)):
                        slc_indices[i] = np.array(slc_indices[i])
                        slc_indices[i][use_nn] = source_coords[i, inearest]

                    setattr(self, k, (Ellipsis,) + tuple(slc_indices))

            assert not mask[self.slice11].any()
            assert not mask[self.slice12].any()
            assert not mask[self.slice21].any()
            assert not mask[self.slice22].any()

        self.idim1 = -2 - ndim_trailing
        self.idim2 = -1 - ndim_trailing

    def __call__(self, fp: np.ndarray) -> np.ndarray:
        assert fp.shape[self.idim1] == self.nxp
        assert fp.shape[self.idim2] == self.nyp
        result = self.w11 * fp[self.slice11]
        result += self.w12 * fp[self.slice12]
        result += self.w21 * fp[self.slice21]
        result += self.w22 * fp[self.slice22]
        return result


class LinearVectorized1D:
    def __init__(
        self,
        x: npt.ArrayLike,
        xp: npt.ArrayLike,
        axis: int = 0,
        fill_value: float = np.nan,
    ):
        x = np.asarray(x)
        xp = np.asarray(xp)
        assert x.ndim == 1
        assert axis >= -xp.ndim and axis < xp.ndim
        xp_slice = [slice(None)] * xp.ndim
        final_shape = list(xp.shape)
        final_shape[axis] = x.size
        ix_left = np.empty(final_shape, dtype=np.intp)
        for ix, cur_x in enumerate(x):
            xp_slice[axis] = ix
            ix_left_cur = (xp < cur_x).sum(axis=axis) - 1
            ix_left_cur += (xp == cur_x).any(axis=axis)
            ix_left[tuple(xp_slice)] = ix_left_cur
        if (np.diff(xp, axis=axis) >= 0.0).all():
            # Source coordinate is monotonically INcreasing
            ix_right = np.minimum(ix_left + 1, xp.shape[axis] - 1)
            ix_left = np.maximum(ix_left, 0)
        else:
            # Source coordinate is monotonically DEcreasing
            ix_right = xp.shape[axis] - 1 - ix_left
            ix_left = np.maximum(ix_right - 1, 0)
            ix_right = np.minimum(ix_right, xp.shape[axis] - 1)
        valid = ix_left != ix_right
        xp_right = np.take_along_axis(xp, ix_right, axis=axis)
        xp_left = np.take_along_axis(xp, ix_left, axis=axis)
        dxp = xp_right - xp_left
        x_shape = [1] * xp.ndim
        x_shape[axis] = x.size
        x_bc = x.reshape(x_shape)
        w_left = np.ones(xp_right.shape)
        np.divide(xp_right - x_bc, dxp, out=w_left, where=valid)
        self.ix_left = ix_left
        self.ix_right = ix_right
        self.w_left = w_left
        self.axis = axis
        self.valid = valid
        self.fill_value = fill_value

    def __call__(self, yp) -> np.ndarray:
        yp = np.asarray(yp)
        yp_left = np.take_along_axis(yp, self.ix_left, axis=self.axis)
        yp_right = np.take_along_axis(yp, self.ix_right, axis=self.axis)
        y = self.w_left * yp_left + (1.0 - self.w_left) * yp_right
        y = np.where(self.valid, y, self.fill_value)
        return y


def interp_1d(
    x: npt.ArrayLike, xp: npt.ArrayLike, fp: npt.ArrayLike, axis: int = 0
) -> np.ndarray:
    """Vectorized 1D linear interpolation along a given axis

    Source values may contain NaNs or masked values at the beginning or end
    of the interpolated dimension; these will be skipped during interpolation.

    Where target coordinates fall outside the range of valid source coordinates,
    the corresponding output values will be equal to the nearest valid source value.

    Args:
        x: Target coordinate values (nD, matching shape of fp except at `axis`)
        xp: Source coordinate values (1D)
        fp: Source values to interpolate (nD, matching the size of xp at `axis`)
        axis: Axis along which to interpolate
    Returns:
        Interpolated values at target coordinates
    """
    x = np.asarray(x, dtype=float)
    xp = np.asarray(xp, dtype=float)
    fp = np.ma.filled(fp, np.nan)
    if fp.ndim != x.ndim:
        raise ValueError(
            f"Number of dimensions {fp.ndim} of source values"
            f" does not match {x.ndim} of target coordinate."
        )
    if xp.ndim != 1:
        raise ValueError(f"Source coordinate must be 1D but has shape {xp.shape}.")
    if fp.shape[:axis] != x.shape[:axis] or fp.shape[axis + 1 :] != x.shape[axis + 1 :]:
        raise ValueError(
            f"Shapes of source values {fp.shape} and target coordinate {x.shape}"
            f" should match everywhere except at the interpolated dimension ({axis})"
        )
    if fp.shape[axis] != xp.shape[0]:
        raise ValueError(
            f"Size of source values {fp.shape[axis]} and source coordinate {xp.shape[0]}"
            f" must match at the interpolated dimension ({axis})"
        )

    dxp = xp[1:] - xp[:-1]
    if not ((dxp > 0).all() or (dxp < 0).all()):
        raise ValueError("xp must be monotonically increasing or decreasing")

    invalid = np.isnan(fp)
    if invalid.any():
        # Source values include NaNs. Identify the inner valid (NaN-free) region.
        valid = ~invalid
        s = tuple(np.newaxis if i != axis else slice(None) for i in range(fp.ndim))
        ind = np.broadcast_to(np.arange(fp.shape[axis])[s], fp.shape)
        first = ind.min(axis=axis, where=valid, initial=fp.shape[axis], keepdims=True)
        last = ind.max(axis=axis, where=valid, initial=0, keepdims=True)
        first = np.minimum(first, last)  # if no valid elements at all, first=last=0
    else:
        first = 0
        last = xp.size - 1

    xp_reversed = dxp[0] < 0
    if xp_reversed:
        # Monotonically DEcreasing source coordinate. Flip it to make it
        # monotonically INcreasing for lookup of bounding intervals
        xp = xp[::-1]
        dxp = -dxp[::-1]
        first, last = (xp.size - 1) - last, (xp.size - 1) - first

    # Look up upper bound of interval around each target coordinate
    # This will be 0 [invalid!] if first source coordinate < minimum target coordinate
    # This will be xp.size [invalid!] if last source coordinate >= maximum target
    # coordinate
    i_right = xp.searchsorted(x, side="right")

    # Determine intervals (left and right bounds).
    # These will be zero-width (i_left == i_right) at the boundaries
    i_left = i_right - 1
    i_left.clip(first, last, out=i_left)
    i_right.clip(first, last, out=i_right)
    scale = np.zeros(x.shape, dtype=float)
    idxp = np.append(1.0 / dxp, 0.0)
    np.multiply(x - xp[i_left], idxp[i_left], where=i_left != i_right, out=scale)

    if xp_reversed:
        # Correct indices into fp for flipped source coordinate
        i_left, i_right = (xp.size - 1) - i_left, (xp.size - 1) - i_right
    f_left = np.take_along_axis(fp, i_left, axis=axis)
    f_right = np.take_along_axis(fp, i_right, axis=axis)
    return f_left + (f_right - f_left) * scale
