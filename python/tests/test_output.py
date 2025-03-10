import unittest
from typing import Optional
import datetime

import numpy as np
import cftime
import netCDF4

import pygetm


class Simulation(pygetm.simulation.BaseSimulation):
    def __init__(self, domain: pygetm.domain.Domain):
        super().__init__(domain, log_level="ERROR")

        self.T = domain.create_grids(1, halox=2, haloy=2, fields=self._fields)

        self.T.array(
            name="month_micro_2d", units="Units", long_name="LongName", fill_value=-2e20
        )
        self.T.array(
            name="month_micro_3d",
            units="Units",
            long_name="LongName",
            z=pygetm.CENTERS,
            fill_value=-2e20,
        )
        self.T.array(
            name="month_macro_2d",
            units="Units",
            long_name="LongName",
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )
        self.T.array(
            name="month_macro_3d",
            units="Units",
            long_name="LongName",
            z=pygetm.CENTERS,
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )

    def _update_forcing_and_diagnostics(self, macro_active: bool):
        self["month_micro_2d"].values[...] = self.time.month
        self["month_micro_3d"].values[...] = self.time.month
        if macro_active:
            self["month_macro_2d"].values[...] = self.time.month
            self["month_macro_3d"].values[...] = self.time.month
        else:
            self["month_macro_2d"].values[...] = np.nan
            self["month_macro_3d"].values[...] = np.nan


class TestOutput(unittest.TestCase):
    def create(self) -> Simulation:
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            interfaces=True,
            f=0.0,
            H=10.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        return Simulation(domain)

    def test_vertical_interpolation(self):
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            interfaces=True,
            f=0.0,
            H=10.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        sim = pygetm.Simulation(
            domain,
            vertical_coordinates=pygetm.vertical_coordinates.Sigma(50),
            airsea=pygetm.airsea.Fluxes(),
        )
        sim.radiation.set_jerlov_type(pygetm.Jerlov.Type_II)
        centers = sim.T.array(
            name="centers",
            units="Units",
            long_name="LongName",
            z=pygetm.CENTERS,
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )
        interfaces = sim.T.array(
            name="interfaces",
            units="Units",
            long_name="LongName",
            z=pygetm.INTERFACES,
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )
        centers.values[...] = (
            np.arange(centers.shape[0])[::-1, np.newaxis, np.newaxis] + 0.5
        )
        interfaces.values[...] = np.arange(interfaces.shape[0])[
            ::-1, np.newaxis, np.newaxis
        ]

        start = cftime.datetime(2000, 1, 1)
        stop = cftime.datetime(2000, 1, 11)
        path = "test.nc"

        nc = sim.output_manager.add_netcdf_file(
            path, interval=datetime.timedelta(days=5)
        )
        nc.request("centers", "interfaces")
        nc.request("centers", z=[0, -5, -10], output_name="centers_ip")
        nc.request("interfaces", z=[0, -5, -10], output_name="interfaces_ip")
        sim.start(start, 3600.0, 24)
        while sim.time < stop:
            sim.advance()
        sim.finish()

        with netCDF4.Dataset(path) as ds:
            cip = ds["centers_ip"]
            self.assertEqual(cip.shape[1], 3)
            self.assertTrue(np.ma.getmaskarray(cip[:, 0, :, :]).all())
            cipmid = cip[:, 1, :, :]
            cipmin = cipmid.min()
            self.assertTrue((cipmid == cipmin).all())
            self.assertAlmostEqual(cipmin, 25.0, places=13)
            self.assertTrue(np.ma.getmaskarray(cip[:, 2, :, :]).all())

            iip = ds["interfaces_ip"]
            self.assertEqual(iip.shape[1], 3)
            # NB not testing k=0 because numerical roundoff when calculating
            # interface positions can cause top interface to be slightly above
            # 0, which causes interpolation to to return masked values.
            iipmid = iip[:, 1, :, :]
            iipmin = iipmid.min()
            self.assertTrue((iipmid == iipmin).all())
            self.assertAlmostEqual(iipmin, 25.0, places=13)
            self.assertTrue((iip[:, 2, :, :] == 50.0).all())

    def test_temporal_averaging(self):
        tol = 1e-14

        def add_file(
            path: str,
            save_initial: bool,
            start: Optional[cftime.datetime] = None,
            stop: Optional[cftime.datetime] = None,
        ):
            nc = sim.output_manager.add_netcdf_file(
                path,
                interval=1,
                interval_units=pygetm.output.TimeUnit.MONTHS,
                start=start,
                stop=stop,
                save_initial=save_initial,
            )
            nc.request(
                "month_micro_2d",
                "month_micro_3d",
                "month_macro_2d",
                "month_macro_3d",
                time_average=True,
            )

        def check_file(
            path: str,
            ntime: int,
            start: cftime.datetime,
            stop: cftime.datetime,
            save_initial: bool,
        ):
            def check_variable(ncvar):
                self.assertEqual(ncvar.long_name, "LongName")
                self.assertEqual(ncvar.units, "Units")
                values = ncvar[...]
                inds = (slice(None),) + (slice(1),) * (values.ndim - 1)
                if save_initial:
                    self.assertTrue(values[0, ...].mask.all())
                    values = values[1:, ...]
                self.assertEqual(values.shape[0], ntime)
                target = values[inds]
                self.assertTrue((values == target).all())
                diff = target.reshape((-1,)) - (np.arange(ntime) % 12 + 1)
                self.assertLess(np.abs(diff).max(), tol)

            with netCDF4.Dataset(path) as ds:
                times = netCDF4.num2date(
                    ds["time"], ds["time"].units, ds["time"].calendar
                )
                self.assertEqual(times[0], start)
                self.assertEqual(times[-1], stop)

                check_variable(ds["month_micro_2d"])
                check_variable(ds["month_macro_2d"])
                check_variable(ds["month_micro_3d"])
                check_variable(ds["month_macro_3d"])

        for save_initial in (True, False):
            for calendar in ("standard", "360_day"):
                with self.subTest(save_initial=save_initial, calendar=calendar):
                    start = cftime.datetime(2000, 1, 1, calendar=calendar)
                    stop = cftime.datetime(2003, 1, 1, calendar=calendar)

                    sim = self.create()

                    add_file("test.nc", save_initial)
                    add_file("test_stop.nc", save_initial, stop=stop.replace(year=2001))
                    add_file(
                        "test_start.nc", save_initial, start=start.replace(year=2002)
                    )
                    add_file(
                        "test_start_stop.nc",
                        save_initial,
                        start=start.replace(year=2001),
                        stop=stop.replace(year=2002),
                    )

                    sim.start(start, 3600.0, 24)
                    while sim.time < stop:
                        sim.advance()
                    sim.finish()

                    start_month = 1 if save_initial else 2

                    check_file(
                        "test.nc",
                        36,
                        start.replace(month=start_month),
                        stop,
                        save_initial,
                    )
                    check_file(
                        "test_stop.nc",
                        12,
                        start.replace(month=start_month),
                        stop.replace(year=2001),
                        save_initial,
                    )
                    check_file(
                        "test_start.nc",
                        12,
                        start.replace(year=2002, month=start_month),
                        stop,
                        save_initial,
                    )
                    check_file(
                        "test_start_stop.nc",
                        12,
                        start.replace(year=2001, month=start_month),
                        stop.replace(year=2002),
                        save_initial,
                    )


if __name__ == "__main__":
    unittest.main()
