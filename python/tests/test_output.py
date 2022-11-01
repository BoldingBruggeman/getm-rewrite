import unittest
import datetime

import numpy as np
import cftime
import netCDF4

import pygetm


class Test(unittest.TestCase):
    def create(self) -> None:
        domain = pygetm.domain.create_cartesian(
            np.linspace(0, 100e3, 50),
            np.linspace(0, 100e3, 51),
            10,
            interfaces=True,
            f=0.0,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )
        domain.logger.parent.setLevel("ERROR")
        domain.initialize(pygetm.BAROTROPIC_2D)

        domain.T.array(
            name="month_micro_2d", units="Units", long_name="LongName", fill_value=-2e20
        )
        domain.T.array(
            name="month_micro_3d",
            units="Units",
            long_name="LongName",
            z=pygetm.CENTERS,
            fill_value=-2e20,
        )
        domain.T.array(
            name="month_macro_2d",
            units="Units",
            long_name="LongName",
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )
        domain.T.array(
            name="month_macro_3d",
            units="Units",
            long_name="LongName",
            z=pygetm.CENTERS,
            fill_value=-2e20,
            attrs=dict(_time_varying=pygetm.TimeVarying.MACRO),
        )

        om = pygetm.output.OutputManager(
            domain.fields,
            logger=domain.logger.parent.getChild("output_manager"),
            rank=0,
        )
        return domain, om

    def loop(
        self,
        domain: pygetm.domain.Domain,
        om: pygetm.output.OutputManager,
        nyear=1,
        dt=3600,
        split_factor=24,
    ):
        def update_vars(macro):
            domain.fields["month_micro_2d"].values[...] = time.month
            domain.fields["month_micro_3d"].values[...] = time.month
            if macro:
                domain.fields["month_macro_2d"].values[...] = time.month
                domain.fields["month_macro_3d"].values[...] = time.month
            else:
                domain.fields["month_macro_2d"].values[...] = np.nan
                domain.fields["month_macro_3d"].values[...] = np.nan

        time = cftime.datetime(2000, 1, 1)
        itime = 0
        macro = True
        timedelta = datetime.timedelta(seconds=dt)
        update_vars(True)
        om.start(time=time)

        while True:
            om.prepare_save(itime * dt, itime, time, macro=macro)
            itime += 1
            time = time + timedelta
            macro = itime % split_factor == 0
            update_vars(macro)
            om.save(itime * dt, itime, time)
            if time.year == 2000 + nyear:
                break
        om.close(itime * dt, time)

    def test_temporal_averaging(self):
        domain, om = self.create()

        save_initial = True
        nc = om.add_netcdf_file(
            "test.nc",
            interval=1,
            interval_units=pygetm.output.TimeUnit.MONTHS,
            save_initial=save_initial,
        )
        nc.request(
            "month_micro_2d",
            "month_micro_3d",
            "month_macro_2d",
            "month_macro_3d",
            time_average=True,
        )

        self.loop(domain, om)

        tol = 1e-14
        ds = netCDF4.Dataset("test.nc")

        def check(ncvar):
            self.assertEqual(ncvar.long_name, "LongName")
            self.assertEqual(ncvar.units, "Units")
            values = ncvar[...]
            inds = (slice(None),) + (slice(1),) * (values.ndim - 1)
            if save_initial:
                self.assertTrue(values[0, ...].mask.all())
                values = values[1:, ...]
            self.assertEqual(values.shape[0], 12)
            target = values[inds]
            self.assertTrue((values == target).all())
            diff = target.reshape((-1,)) - np.arange(1, 13)
            self.assertLess(np.abs(diff).max(), tol)

        check(ds["month_micro_2d"])
        check(ds["month_macro_2d"])
        check(ds["month_micro_3d"])
        check(ds["month_macro_3d"])


if __name__ == "__main__":
    unittest.main()
