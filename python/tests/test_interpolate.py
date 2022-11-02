import unittest

import numpy as np
import scipy.interpolate

from pygetm.util.interpolate import Linear2DGridInterpolator, interp_1d


def generate_random_horizontal_grid(nx=99, ny=100):
    # Random source grid (1D x and y)
    # Still monotonically increasing because of the use of cumsum
    dx = np.random.random_sample((nx,))
    dy = np.random.random_sample((ny,))
    xp = np.cumsum(dx)
    yp = np.cumsum(dy)

    # Random source field
    fp = np.random.random_sample((xp.size, yp.size,))
    return xp, yp, fp


class TestInterpolate(unittest.TestCase):
    # Maximum acceptable error
    EPS = 5 * np.finfo(float).eps

    def test_horizontal_corners(self):
        xp, yp, fp = generate_random_horizontal_grid()

        def compare_index(i: int, j: int):
            f = Linear2DGridInterpolator(xp[i], yp[j], xp, yp)(fp)
            self.assertEqual(
                f, fp[i, j], "mismatch @ x=%i, y=%i" % (i, j),
            )

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
                "%s - original order" % name,
            )

            f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp)(fp[::-1, :])
            self.assertTrue((f - f2 == 0).all(), "%s - x reversed" % name)

            f2 = Linear2DGridInterpolator(x, y, xp, yp[::-1])(fp[:, ::-1])
            self.assertTrue((f - f2 == 0).all(), "%s - y reversed" % name)

            f2 = Linear2DGridInterpolator(x, y, xp[::-1], yp[::-1])(fp[::-1, ::-1])
            self.assertTrue((f - f2 == 0).all(), "%s - xy reversed" % name)

        # Interpolate to 2D target grid
        # and compare results with that of scipy.interpolate.interpn
        shape = (50, 51)
        x = np.random.uniform(xp[0], xp[-1], shape)
        y = np.random.uniform(yp[0], yp[-1], shape)
        compare_spatial("2D", x, y)

        # Interpolate to 1D target grid
        # and compare results with that of scipy.interpolate.interpn
        shape = (100,)
        x = np.random.uniform(xp[0], xp[-1], shape)
        y = np.random.uniform(yp[0], yp[-1], shape)
        compare_spatial("1D", x, y)

    def test_vertical_1d(self):
        dx = np.random.random_sample((99,))
        xp = np.cumsum(dx)
        fp = np.random.random_sample((xp.size,))
        x = np.random.uniform(xp[0] - 5, xp[-1] + 5, (100,))
        f = interp_1d(x, xp, fp)
        f_check = np.interp(x, xp, fp)
        self.assertTrue(
            np.isclose(f, f_check, rtol=self.EPS, atol=self.EPS).all(),
            "1D - original order",
        )

        f2 = interp_1d(x, xp[::-1], fp[::-1])
        self.assertTrue((f == f2).all(), "1D - reversed")

    def test_vertical_3d(self):
        nx, ny = 5, 6
        dx = np.random.random_sample((99,))
        xp = np.cumsum(dx)
        fp = np.random.random_sample((xp.size, ny, nx))

        x = -10 * np.random.random_sample((ny, nx)) + 1.5 * np.random.random_sample(
            (100, ny, nx)
        ).cumsum(axis=0)
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
        dx = np.random.random_sample((99,))
        xp = np.cumsum(dx)
        fp = np.random.random_sample((xp.size, ny, nx))

        x = -10 * np.random.random_sample((ny, nx)) + 1.5 * np.random.random_sample(
            (100, ny, nx)
        ).cumsum(axis=0)

        start = np.random.randint(0, fp.shape[0] - 1, fp.shape[1:])
        stop = np.random.randint(start, fp.shape[0])
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


if __name__ == "__main__":
    unittest.main()
