import unittest

import numpy as np
import scipy.interpolate

from pygetm.util.interpolate import (
    Linear2DGridInterpolator,
    LinearVectorized1D,
    interp_1d,
)
from pygetm.constants import EdgeTreatment


def generate_random_horizontal_grid(nx=99, ny=100):
    # Random source grid (1D x and y)
    # Still monotonically increasing because of the use of cumsum
    rng = np.random.default_rng()
    dx = rng.random(nx)
    dy = rng.random(ny)
    xp = np.cumsum(dx)
    yp = np.cumsum(dy)

    # Random source field
    fp = rng.random((xp.size, yp.size))
    return xp, yp, fp


class TestInterpolate(unittest.TestCase):
    # Maximum acceptable error
    EPS = 5 * np.finfo(float).eps

    def test_horizontal_corners(self):
        xp, yp, fp = generate_random_horizontal_grid()

        def compare_index(i: int, j: int):
            f = Linear2DGridInterpolator(xp[i], yp[j], xp, yp)(fp)
            self.assertEqual(f, fp[i, j], f"mismatch @ x={i}, y={j}")

        compare_index(0, 0)
        compare_index(0, -1)
        compare_index(-1, 0)
        compare_index(-1, -1)

    def test_horizontal_to_original_grid(self):
        xp, yp, fp = generate_random_horizontal_grid()
        x, y = np.broadcast_arrays(xp[:, np.newaxis], yp[np.newaxis, :])
        f = Linear2DGridInterpolator(x, y, xp, yp)(fp)
        self.assertTrue((f == fp).all(), "interpolate to original grid")

    def test_horizontal_compare_with_scipy(self):
        xp, yp, fp = generate_random_horizontal_grid()

        def compare_spatial(name: str, x, y):
            f_check = scipy.interpolate.interpn((xp, yp), fp, (x, y))
            f = Linear2DGridInterpolator(x, y, xp, yp)(fp)
            self.assertTrue(
                np.isclose(f, f_check, rtol=self.EPS, atol=self.EPS).all(),
                f"{name} - original order",
            )

            f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp)(fp[::-1, :])
            self.assertTrue((f - f2 == 0).all(), f"{name} - x reversed")

            f2 = Linear2DGridInterpolator(x, y, xp, yp[::-1])(fp[:, ::-1])
            self.assertTrue((f - f2 == 0).all(), f"{name} - y reversed")

            f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp[::-1])(fp[::-1, ::-1])
            self.assertTrue((f - f2 == 0).all(), f"{name} - xy reversed")

        # Interpolate to 2D target grid
        # and compare results with that of scipy.interpolate.interpn
        shape = (50, 51)
        rng = np.random.default_rng()
        x = rng.uniform(xp[0], xp[-1], shape)
        y = rng.uniform(yp[0], yp[-1], shape)
        compare_spatial("2D", x, y)

        # Interpolate to 1D target grid
        # and compare results with that of scipy.interpolate.interpn
        shape = (100,)
        x = rng.uniform(xp[0], xp[-1], shape)
        y = rng.uniform(yp[0], yp[-1], shape)
        compare_spatial("1D", x, y)

    def test_vertical_1d(self):
        rng = np.random.default_rng()
        dx = rng.random(99)
        xp = np.cumsum(dx)
        fp = rng.random(xp.size)
        x = rng.uniform(xp[0] - 5, xp[-1] + 5, (100,))
        f = interp_1d(x, xp, fp)
        f_check = np.interp(x, xp, fp)
        self.assertTrue(
            np.isclose(f, f_check, rtol=self.EPS, atol=self.EPS).all(),
            "1D - original order",
        )

        f2 = interp_1d(x, xp[::-1], fp[::-1])
        self.assertTrue((f == f2).all(), "1D - reversed")

        f3 = interp_1d(xp, xp, fp)
        self.assertTrue((f3 == fp).all(), "1D - to original grid")

        f4 = interp_1d(xp[::-1], xp, fp)
        self.assertTrue((f4 == fp[::-1]).all(), "1D - to reversed original grid")

    def test_vertical_3d(self):
        nx, ny = 5, 6
        rng = np.random.default_rng()
        dx = rng.random(99)
        xp = np.cumsum(dx)
        fp = rng.random((xp.size, ny, nx))

        x = -10 * rng.random((ny, nx)) + 1.5 * rng.random((100, ny, nx)).cumsum(axis=0)
        f = interp_1d(x, xp, fp)
        self.assertTrue(x.shape == f.shape)
        f_check = np.empty_like(f)
        for i in range(x.shape[-1]):
            for j in range(x.shape[-2]):
                f_check[:, j, i] = np.interp(x[:, j, i], xp, fp[:, j, i])
        self.assertTrue(
            np.isclose(f, f_check, rtol=self.EPS, atol=self.EPS).all(), "3D"
        )

        f2 = interp_1d(x, xp[::-1], fp[::-1, ...])
        self.assertTrue((f == f2).all(), "3d - xy reversed")

    def test_vertical_3d_masked(self):
        nx, ny = 5, 6
        rng = np.random.default_rng()
        dx = rng.random(99)
        xp = np.cumsum(dx)
        fp = rng.random((xp.size, ny, nx))

        x = -10 * rng.random((ny, nx)) + 1.5 * rng.random((100, ny, nx)).cumsum(axis=0)

        start = rng.integers(0, fp.shape[0] - 1, fp.shape[1:])
        stop = rng.integers(start, fp.shape[0])
        ind = np.broadcast_to(
            np.arange(fp.shape[0])[:, np.newaxis, np.newaxis], fp.shape
        )
        mask = np.logical_or(ind < start, ind >= stop)
        masked_fp = np.ma.array(fp, mask=mask)
        f = interp_1d(x, xp, masked_fp)
        for i in range(x.shape[-1]):
            for j in range(x.shape[-2]):
                if start[j, i] == stop[j, i]:
                    # no valid points - all output should be masked
                    self.assertTrue(np.isnan(f[:, j, i]).all())
                else:
                    self.assertTrue(
                        np.isfinite(fp[start[j, i] : stop[j, i], j, i]).all()
                    )
                    f_check = np.interp(
                        x[:, j, i],
                        xp[start[j, i] : stop[j, i]],
                        fp[start[j, i] : stop[j, i], j, i],
                    )
                    self.assertTrue(
                        np.isclose(
                            f[:, j, i], f_check, rtol=self.EPS, atol=self.EPS
                        ).all()
                    )

    def test_vertical_3d_masked_lastaxis(self):
        axis = 2
        nx, ny = 5, 6
        rng = np.random.default_rng()
        dx = rng.random(99)
        xp = np.cumsum(dx)
        shape = [ny, nx]
        shape.insert(axis, xp.size)
        fp = rng.random(shape)

        shape[axis] = 100
        shape2 = list(shape)
        shape2[axis] = 1
        x = -10 * rng.random(shape2) + 1.5 * rng.random(shape).cumsum(axis=axis)

        start = rng.integers(
            0, fp.shape[axis] - 1, fp.shape[:axis] + fp.shape[axis + 1 :]
        )
        stop = rng.integers(start, fp.shape[axis])
        ind = np.arange(fp.shape[axis])
        slc = [np.newaxis] * fp.ndim
        slc[axis] = slice(None)
        ind = np.broadcast_to(ind[tuple(slc)], fp.shape)
        slc = [slice(None)] * fp.ndim
        slc[axis] = np.newaxis
        fullstart = start[tuple(slc)]
        fullstop = stop[tuple(slc)]
        mask = np.logical_or(ind < fullstart, ind >= fullstop)
        masked_fp = np.ma.array(fp, mask=mask)
        f = interp_1d(x, xp, masked_fp, axis=axis)

        def make_slice(i, j, k=slice(None)):
            slc = [j, i]
            slc.insert(axis, k)
            return tuple(slc)

        for i in range(nx):
            for j in range(ny):
                if start[j, i] == stop[j, i]:
                    # no valid points - all output should be masked
                    self.assertTrue(np.isnan(f[make_slice(i, j)]).all())
                else:
                    self.assertTrue(
                        np.isfinite(
                            fp[make_slice(i, j, slice(start[j, i], stop[j, i]))]
                        ).all()
                    )
                    f_check = np.interp(
                        x[make_slice(i, j)],
                        xp[start[j, i] : stop[j, i]],
                        fp[make_slice(i, j, slice(start[j, i], stop[j, i]))],
                    )
                    self.assertTrue(
                        np.isclose(
                            f[make_slice(i, j)], f_check, rtol=self.EPS, atol=self.EPS
                        ).all()
                    )

    def test_linear_vectorized_1d(self):
        x = np.array([0.0, -50.0, -100.0])
        rng = np.random.default_rng()
        H_mean = 100.0
        h = H_mean / 0.5 / 30.0 * rng.random((30, 10, 11))
        elev = rng.random(h.shape[1:]) - 0.5
        H = h.sum(axis=0) - elev
        z = np.zeros((h.shape[0] + 1,) + h.shape[1:])
        z[1:] = h.cumsum(axis=0)
        assert ((z[1:] - z[:-1]) > 0.0).all()
        z -= H
        FILL_VALUE = -999.0
        Z_FILL_VALUE = -9999.0
        values = rng.random(z.shape)
        for mask in (True, False):
            masked = rng.random(H.shape) < 0.5 if mask else False
            masked = np.broadcast_to(masked, H.shape)
            z_masked = np.where(masked, Z_FILL_VALUE, z)
            for edge in (EdgeTreatment.MISSING, EdgeTreatment.CLAMP):
                with self.subTest(mask=mask, edge=edge):
                    ip = LinearVectorized1D(
                        x,
                        z_masked,
                        axis=0,
                        fill_value=FILL_VALUE,
                        mask=z_masked == Z_FILL_VALUE,
                        edges=edge,
                    )
                    f = ip(values)
                    assert f.shape[0] == x.size
                    for x_ip, v_ip in zip(x, f):
                        too_deep = (x_ip < z_masked[0]) & ~masked
                        too_shallow = (x_ip > z_masked[-1]) & ~masked
                        valid = ~(too_deep | too_shallow | masked)
                        deep_oob = (
                            FILL_VALUE if edge == EdgeTreatment.MISSING else values[0]
                        )
                        shallow_oob = (
                            FILL_VALUE if edge == EdgeTreatment.MISSING else values[-1]
                        )
                        self.assertTrue((v_ip == deep_oob).all(where=too_deep))
                        self.assertTrue((v_ip == shallow_oob).all(where=too_shallow))
                        self.assertTrue((v_ip == FILL_VALUE).all(where=masked))
                        self.assertTrue((v_ip != FILL_VALUE).all(where=valid))
                    for i in range(z.shape[-1]):
                        for j in range(z.shape[-2]):
                            if masked[j, i]:
                                self.assertTrue((f[:, j, i] == FILL_VALUE).all())
                            else:
                                oob = (
                                    None if edge == EdgeTreatment.CLAMP else FILL_VALUE
                                )
                                f_check = np.interp(
                                    x,
                                    z_masked[:, j, i],
                                    values[:, j, i],
                                    left=oob,
                                    right=oob,
                                )
                                self.assertTrue(
                                    np.isclose(
                                        f[:, j, i],
                                        f_check,
                                        rtol=self.EPS,
                                        atol=self.EPS,
                                    ).all()
                                )


if __name__ == "__main__":
    unittest.main()
