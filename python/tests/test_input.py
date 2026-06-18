import unittest
import datetime
import logging

import numpy as np
import cftime
import xarray as xr

import pygetm

NDEC = 14  # number of decimal places to check


class TestInput(unittest.TestCase):
    def test_temporal_interpolation(self):
        idays = np.arange(365)
        data = np.random.random(365)
        dates = cftime.num2date(np.arange(365), "days since 2015-01-01 00:00:00")
        array = xr.DataArray(data, {"time": dates}, name="test_series")

        xip = pygetm.input.temporal_interpolation(array)
        xdata: pygetm.input.TemporalInterpolation = xip.variable._data
        self.assertIsInstance(xdata, pygetm.input.TemporalInterpolation)
        self.assertTrue(xdata.is_time_varying())

        # Before time series starts
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        time = cftime.datetime(2014, 12, 31)
        with self.assertRaisesRegex(
            Exception,
            (
                "Cannot interpolate test_series to value at 2014-12-31 00:00:00,"
                " because time series starts only at 2015-01-01 00:00:00."
            ),
        ):
            xdata.update(time, time.toordinal(fractional=True))

        # At the exact start of the time series
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        self.assertTrue(xdata.update(cftime.datetime(2015, 1, 1)))

        # At the exact end of the time series
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        self.assertTrue(xdata.update(cftime.datetime(2015, 12, 31)))

        # Beyond the end of the time series
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        time = cftime.datetime(2016, 1, 1)
        with self.assertRaisesRegex(
            Exception,
            (
                "Cannot interpolate test_series to value at 2016-01-01 00:00:00"
                " because end of time series was reached \\(2015-12-31 00:00:00\\)."
            ),
        ):
            xdata.update(time, time.toordinal(fractional=True))

        # Specific points within the time series
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        time = cftime.datetime(2015, 1, 1)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertFalse(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xdata), data[0])

        time = cftime.datetime(2015, 5, 8)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertFalse(xdata.update(time, time.toordinal(fractional=True)))
        nday = (time - cftime.datetime(2015, 1, 1)).days
        self.assertEqual(np.asarray(xdata), np.interp(nday, idays, data))

        time = time + datetime.timedelta(hours=6)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertAlmostEqual(
            np.asarray(xdata), np.interp(nday + 0.25, idays, data), places=NDEC
        )

        # Try rewinding
        time = cftime.datetime(2015, 5, 8)
        xdata.update(time, time.toordinal(fractional=True))
        self.assertEqual(np.asarray(xdata), np.interp(nday, idays, data))

        time = cftime.datetime(2015, 12, 31)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xdata), data[-1])

        # Beyond the end of the time series
        time = cftime.datetime(2016, 1, 1)
        with self.assertRaisesRegex(
            Exception,
            (
                "Cannot interpolate test_series to value at 2016-01-01 00:00:00"
                " because end of time series was reached \\(2015-12-31 00:00:00\\)."
            ),
        ):
            xdata.update(time, time.toordinal(fractional=True))

        xip = pygetm.input.temporal_interpolation(array, climatology=True)
        xdata: pygetm.input.TemporalInterpolation = xip.variable._data
        self.assertIsInstance(xdata, pygetm.input.TemporalInterpolation)
        self.assertTrue(xdata.is_time_varying())

        for year in range(2000, 2030):
            self.assertTrue(xdata.update(cftime.datetime(year, 1, 1)))
            self.assertEqual(np.asarray(xdata), data[0])

            time = cftime.datetime(year, 5, 8)
            self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
            self.assertEqual(np.asarray(xdata), np.interp(nday, idays, data))

        data = np.random.random(12)
        dates = [cftime.datetime(2015, month, 16) for month in range(1, 13)]
        array = xr.DataArray(data, {"time": dates}, name="test_series")

        # Try start times just before and just at the start of the time series
        xdata = pygetm.input.temporal_interpolation(
            array, climatology=True
        ).variable._data
        self.assertTrue(xdata.update(cftime.datetime(year, 1, 1)))
        xdata = pygetm.input.temporal_interpolation(
            array, climatology=True
        ).variable._data
        self.assertTrue(xdata.update(cftime.datetime(year, 1, 16)))

        # Try start times just after and just at the end of the time series
        xdata = pygetm.input.temporal_interpolation(
            array, climatology=True
        ).variable._data
        self.assertTrue(xdata.update(cftime.datetime(year, 12, 31)))
        xdata = pygetm.input.temporal_interpolation(
            array, climatology=True
        ).variable._data
        self.assertTrue(xdata.update(cftime.datetime(year, 12, 16)))

        xip = pygetm.input.temporal_interpolation(array, climatology=True)
        xdata: pygetm.input.TemporalInterpolation = xip.variable._data
        self.assertIsInstance(xdata, pygetm.input.TemporalInterpolation)
        self.assertTrue(xdata.is_time_varying())

        for year in range(2000, 2030):
            self.assertTrue(xdata.update(cftime.datetime(year, 1, 1)))
            self.assertAlmostEqual(
                np.asarray(xdata),
                15.0 / 31.0 * data[-1] + 16.0 / 31.0 * data[0],
                places=NDEC,
            )

            self.assertTrue(xdata.update(cftime.datetime(year, 1, 16)))
            self.assertEqual(np.asarray(xdata), data[0])

            self.assertTrue(xdata.update(cftime.datetime(year, 5, 1)))
            self.assertAlmostEqual(
                np.asarray(xdata), 0.5 * (data[3] + data[4]), places=NDEC
            )

            self.assertTrue(xdata.update(cftime.datetime(year, 10, 1)))
            self.assertAlmostEqual(
                np.asarray(xdata), 0.5 * (data[8] + data[9]), places=NDEC
            )

        # Different calendars

        data = np.random.random(360)
        days = np.arange(data.size)
        ref = "days since 2015-01-01 00:00:00"

        dates = cftime.num2date(days, ref, calendar="proleptic_gregorian")
        array = xr.DataArray(data, {"time": dates}, name="test_series")
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        self.assertTrue(xdata.update(cftime.datetime(2015, 6, 1, calendar="standard")))
        self.assertTrue(xdata.update(cftime.datetime(2015, 8, 1, calendar="standard")))
        self.assertTrue(xdata.update(cftime.datetime(2015, 2, 1, calendar="standard")))

        dates = cftime.num2date(days, ref, calendar="proleptic_gregorian")
        array = xr.DataArray(data, {"time": dates}, name="test_series")
        ip = pygetm.input.temporal_interpolation(array, climatology=True)
        xdata = ip.variable._data
        self.assertTrue(xdata.update(cftime.datetime(2015, 6, 1, calendar="standard")))
        self.assertTrue(xdata.update(cftime.datetime(2018, 7, 1, calendar="standard")))
        self.assertTrue(xdata.update(cftime.datetime(2015, 2, 1, calendar="standard")))

        dates = cftime.num2date(days, ref, calendar="360_day")
        array = xr.DataArray(data, {"time": dates}, name="test_series")
        xdata = pygetm.input.temporal_interpolation(array).variable._data
        self.assertTrue(xdata.update(cftime.datetime(2015, 6, 1, calendar="360_day")))

        for clim in (True, False):
            with self.assertRaisesRegex(
                Exception,
                (
                    "Calendar 360_day used by test_series is incompatible"
                    " with simulation calendar standard."
                ),
            ):
                dates = cftime.num2date(np.arange(data.size), ref, calendar="360_day")
                array = xr.DataArray(data, {"time": dates}, name="test_series")
                ip = pygetm.input.temporal_interpolation(array, climatology=clim)
                xdata = ip.variable._data
                self.assertTrue(
                    xdata.update(cftime.datetime(2018, 6, 1, calendar="standard"))
                )

    def test_concatenate(self):
        axis = 1
        a1 = np.random.random((100, 10))
        a2 = np.random.random((100, 5))
        a = pygetm.input.Concatenate((a1, a2), axis=axis)
        self.assertTrue(a.shape == (100, a1.shape[axis] + a2.shape[1]))
        data = np.asarray(a)
        self.assertTrue((data == np.concatenate((a1, a2), axis=axis)).all())
        self.assertTrue((data[:, : a1.shape[axis]] == a1).all())
        self.assertTrue((data[:, a1.shape[axis] :] == a2).all())
        data2 = a[:, :]
        self.assertTrue((data == data2).all())
        self.assertTrue((a[:, 5] == a1[:, 5]).all())
        self.assertTrue((a[:, 12] == a2[:, 2]).all())
        data_slc = a[0, ...]
        self.assertTrue((data_slc == data[0, ...]).all())

        axis = 0
        a1 = np.random.random((10, 100))
        a2 = np.random.random((5, 100))
        a = pygetm.input.Concatenate((a1, a2), axis=axis)
        self.assertTrue(a.shape == (a1.shape[axis] + a2.shape[axis], 100))
        data = np.asarray(a)
        self.assertTrue((data == np.concatenate((a1, a2), axis=axis)).all())
        self.assertTrue((data[: a1.shape[axis], :] == a1).all())
        self.assertTrue((data[a1.shape[axis] :, :] == a2).all())
        data2 = a[:, :]
        self.assertTrue((data == data2).all())
        self.assertTrue((a[5, :] == a1[5, :]).all())
        self.assertTrue((a[12, :] == a2[2, :]).all())
        data_slc = a[..., 0]
        self.assertTrue((data_slc == data[..., 0]).all())


class TestManager(unittest.TestCase):
    def test_open_boundary(self):
        logger = logging.getLogger()
        logger.setLevel("ERROR")
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100, 10),
            [0.0, 1.0],
            interfaces=True,
            f=0,
            H=10,
            logger=logger,
        )
        domain.open_boundaries.add_right_boundary(
            "mouth", 8, 0, 1, type_2d=pygetm.FLATHER_ELEV, type_3d=pygetm.ZERO_GRADIENT
        )
        grid = domain.create_grids(
            10,
            0,
            0,
            input_manager=pygetm.input.InputManager(
                domain.root_logger.getChild("input")
            ),
        )
        arr2d = grid.array("foo2d", on_boundary=True)
        arr2d.set(1.0)
        values = np.arange(10)[:, np.newaxis]
        time = [
            cftime.datetime(2015, 1, 1) + datetime.timedelta(days=i) for i in range(10)
        ]
        forcing2d = xr.DataArray(values, dims=["time", "bdy"], coords={"time": time})
        arr2d.set(forcing2d)
        grid.input_manager.update(time[0], time[0].toordinal(fractional=True))

        arr3d = grid.array("foo3d", z=pygetm.CENTERS, on_boundary=True)
        arr3d.set(1.0)
        values = np.broadcast_to(np.arange(10)[:, np.newaxis, np.newaxis], (10, 1, 20))
        z = xr.Variable(("z",), np.linspace(0, 10, 20), dict(axis="Z"))
        forcing3d = xr.DataArray(
            values, dims=["time", "bdy", "z"], coords={"time": time, "z": z}
        )
        arr3d.set(forcing3d)
        grid.input_manager.update(time[0], time[0].toordinal(fractional=True))

        values = np.broadcast_to(np.arange(10)[:, np.newaxis, np.newaxis], (10, 20, 1))
        forcing3d = xr.DataArray(
            values, dims=["time", "z", "bdy"], coords={"time": time, "z": z}
        )
        arr3d.set(forcing3d)
        grid.input_manager.update(time[0], time[0].toordinal(fractional=True))


if __name__ == "__main__":
    unittest.main()
