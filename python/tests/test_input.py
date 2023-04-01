import unittest
import datetime

import numpy as np
import cftime
import xarray

import pygetm


class TestInput(unittest.TestCase):
    def test_temporal_interpolation(self):
        data = np.arange(365, dtype=float)
        dates = cftime.num2date(np.arange(365), "days since 2015-01-01 00:00:00")
        array = xarray.DataArray(data, {"time": dates}, name="test_series")
        xip = pygetm.input.temporal_interpolation(array)

        xdata: pygetm.input.TemporalInterpolation = xip.variable._data
        self.assertIsInstance(xdata, pygetm.input.TemporalInterpolation)
        self.assertTrue(xdata.is_time_varying())

        time = cftime.datetime(2014, 12, 31)
        with self.assertRaisesRegex(
            Exception,
            (
                "Cannot interpolate test_series to value at 2014-12-31 00:00:00,"
                " because time series starts only at 2015-01-01 00:00:00."
            ),
        ):
            xdata.update(time, time.toordinal(fractional=True))

        time = cftime.datetime(2015, 1, 1)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertFalse(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xip), 0.0)

        time = cftime.datetime(2015, 5, 8)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertFalse(xdata.update(time, time.toordinal(fractional=True)))
        nday = (time - cftime.datetime(2015, 1, 1)).days
        self.assertEqual(np.asarray(xip), nday)

        time = time + datetime.timedelta(hours=6)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xip), nday + 0.25)

        time = cftime.datetime(2015, 5, 8)
        with self.assertRaisesRegex(
            Exception,
            (
                "Time can only increase, but previous time was 2015-05-08 06:00:00,"
                " new time 2015-05-08 00:00:00"
            ),
        ):
            xdata.update(time, time.toordinal(fractional=True))

        time = cftime.datetime(2015, 12, 31)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xip), 364.0)

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
            self.assertEqual(np.asarray(xip), 0.0)

            time = cftime.datetime(year, 5, 8)
            self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
            self.assertEqual(np.asarray(xip), nday)

        data = np.arange(12, dtype=float)
        dates = [cftime.datetime(2015, month, 16) for month in range(1, 13)]
        array = xarray.DataArray(data, {"time": dates}, name="test_series")
        xip = pygetm.input.temporal_interpolation(array, climatology=True)
        xdata: pygetm.input.TemporalInterpolation = xip.variable._data
        self.assertIsInstance(xdata, pygetm.input.TemporalInterpolation)
        self.assertTrue(xdata.is_time_varying())

        for year in range(2000, 2030):
            self.assertTrue(xdata.update(cftime.datetime(year, 1, 1)))
            self.assertAlmostEqual(np.asarray(xip), 15.0 / 31.0 * 11, places=14)

            self.assertTrue(xdata.update(cftime.datetime(year, 1, 16)))
            self.assertEqual(np.asarray(xip), 0.0)

            self.assertTrue(xdata.update(cftime.datetime(year, 5, 1)))
            self.assertEqual(np.asarray(xip), 3.5)

            self.assertTrue(xdata.update(cftime.datetime(year, 10, 1)))
            self.assertEqual(np.asarray(xip), 8.5)


if __name__ == "__main__":
    unittest.main()
