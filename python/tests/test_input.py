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
        self.assertEqual(np.asarray(xdata), 0.0)

        time = cftime.datetime(2015, 5, 8)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertFalse(xdata.update(time, time.toordinal(fractional=True)))
        nday = (time - cftime.datetime(2015, 1, 1)).days
        self.assertEqual(np.asarray(xdata), nday)

        time = time + datetime.timedelta(hours=6)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xdata), nday + 0.25)

        # Try rewinding
        time = cftime.datetime(2015, 5, 8)
        xdata.update(time, time.toordinal(fractional=True))
        self.assertEqual(np.asarray(xdata), nday)

        time = cftime.datetime(2015, 12, 31)
        self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
        self.assertEqual(np.asarray(xdata), 364.0)

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
            self.assertEqual(np.asarray(xdata), 0.0)

            time = cftime.datetime(year, 5, 8)
            self.assertTrue(xdata.update(time, time.toordinal(fractional=True)))
            self.assertEqual(np.asarray(xdata), nday)

        data = np.arange(12, dtype=float)
        dates = [cftime.datetime(2015, month, 16) for month in range(1, 13)]
        array = xarray.DataArray(data, {"time": dates}, name="test_series")

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
            self.assertAlmostEqual(np.asarray(xdata), 15.0 / 31.0 * 11, places=14)

            self.assertTrue(xdata.update(cftime.datetime(year, 1, 16)))
            self.assertEqual(np.asarray(xdata), 0.0)

            self.assertTrue(xdata.update(cftime.datetime(year, 5, 1)))
            self.assertEqual(np.asarray(xdata), 3.5)

            self.assertTrue(xdata.update(cftime.datetime(year, 10, 1)))
            self.assertEqual(np.asarray(xdata), 8.5)


if __name__ == "__main__":
    unittest.main()
