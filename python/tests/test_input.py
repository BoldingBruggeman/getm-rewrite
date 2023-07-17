import unittest
import datetime

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


if __name__ == "__main__":
    unittest.main()
