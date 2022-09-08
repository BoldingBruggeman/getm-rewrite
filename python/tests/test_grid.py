import unittest

import numpy as np

import pygetm


class Test(unittest.TestCase):
    def test_interpolation(self):
        domain = pygetm.domain.create_cartesian(
            np.linspace(0.0, 10000.0, 101),
            np.linspace(0.0, 10000.0, 100),
            nz=30,
            H=100,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        domain.initialize(pygetm.BAROCLINIC)

        for z in (False, pygetm.CENTERS, pygetm.INTERFACES):
            with self.subTest(z=z):
                # Random initialization of T
                t = domain.T.array(z=pygetm.CENTERS)
                t.all_values[...] = np.random.random(t.all_values.shape)

                # From T to U
                u = t.interp(domain.U)
                u_control = 0.5 * (t.all_values[..., :, :-1] + t.all_values[..., :, 1:])
                self.assertTrue((u.all_values[..., :, :-1] == u_control).all())

                # From T to V
                v = t.interp(domain.V)
                v_control = 0.5 * (t.all_values[..., :-1, :] + t.all_values[..., 1:, :])
                self.assertTrue((v.all_values[..., :-1, :] == v_control).all())

                # From T to X
                x = t.interp(domain.X)
                x_control = 0.25 * (
                    t.all_values[..., :-1, :-1]
                    + t.all_values[..., :-1, 1:]
                    + t.all_values[..., 1:, :-1]
                    + t.all_values[..., 1:, 1:]
                )
                self.assertTrue((x.all_values[..., 1:-1, 1:-1] == x_control).all())

                # Random initialization of X
                x.all_values[...] = np.random.random(x.all_values.shape)

                # From X to T
                t = x.interp(domain.T)
                t_control = 0.25 * (
                    x.all_values[..., :-1, :-1]
                    + x.all_values[..., :-1, 1:]
                    + x.all_values[..., 1:, :-1]
                    + x.all_values[..., 1:, 1:]
                )
                self.assertTrue((t.all_values == t_control).all())

                # Random initialization of U
                u.all_values[...] = np.random.random(u.all_values.shape)

                # From U to UU
                uu = u.interp(domain.UU)
                uu_control = 0.5 * (
                    u.all_values[..., :, :-1] + u.all_values[..., :, 1:]
                )
                self.assertTrue((uu.all_values[..., :, :-1] == uu_control).all())

                # From U to VU
                vu = u.interp(domain.VU)
                vu_control = 0.5 * (
                    u.all_values[..., :-1, :] + u.all_values[..., 1:, :]
                )
                self.assertTrue((vu.all_values[..., :-1, :] == vu_control).all())

                # From U to T
                t = u.interp(domain.T)
                t_control = 0.5 * (u.all_values[..., :, :-1] + u.all_values[..., :, 1:])
                self.assertTrue((t.all_values[..., :, 1:] == t_control).all())

                # Random initialization of V
                v.all_values[...] = np.random.random(v.all_values.shape)

                # From V to UV
                uv = v.interp(domain.UV)
                uv_control = 0.5 * (
                    v.all_values[..., :, :-1] + v.all_values[..., :, 1:]
                )
                self.assertTrue((uv.all_values[..., :, :-1] == uv_control).all())

                # From V to VV
                vv = v.interp(domain.VV)
                vv_control = 0.5 * (
                    v.all_values[..., :-1, :] + v.all_values[..., 1:, :]
                )
                self.assertTrue((vv.all_values[..., :-1, :] == vv_control).all())

                # From V to T
                t = v.interp(domain.T)
                t_control = 0.5 * (v.all_values[..., :-1, :] + v.all_values[..., 1:, :])
                self.assertTrue((t.all_values[..., 1:, :] == t_control).all())


if __name__ == "__main__":
    unittest.main()
