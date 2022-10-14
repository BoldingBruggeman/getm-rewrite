# This verifies whether all discrepancies between te 2D and 3D update disappear
# when split_factor=1 and nz=1. This means that:
# * the the difference between 2D and 3D transports at the end of the baroclinic
#   update (the source of the "correction" that normally is applied to 3D
#   transports) must be 0. Thus udev and vdev (the difference in depth-averaged
#   velocities before correction) must be 0.
# * advection/diffusion/bottom friction slow terms (S?A, S?D, S?F) must be 0
# * internal pressure slows terms (S?B) must equal their depth-explicit
#   equivalent (idpd?)
#
# Challenges:
# * 2D advection by default uses a type of splitting (either u/2-v-u/2 or u-v)
#   that is incompatible with the default type of splitting (u/2-v/2-w-v/2-u/2)
#   used in 3D. *Workaround*: set advection_split_2d to AdvectionSplit.HALF_ALWAYS
#   to force 2D to use u/2-v/2-v/2-u/2, which currently serves no purpose beyond
#   enabling 2D-3D convergence
# * The slow bottom friction term in the 2D momentum update is clipped to always
#   dampen the 2D velocities (its inclusion depends on the sign of u).
#   This behavior is hardcoded in master, and the default in iow. It prevents
#   2D-3D convergence/. *Workaround*: set Slr should be set equal SxF/SyF in
#   momentum_2d.F90 (as with iow's _SLR_NOCLIP_)
# * The implicit bottom friction term in the momentum update uses old thickness
#   ho in 2D but new thickness in 3D. 2D cannot use the new thickness because
#   this is not known yet. *Workaround* (in code): force 3D to use the old thickness
#   by multiplying rru and rrv with hn/ho when setting ea2 in momentum_3d.F90.
#   This allows convergence but could lower accuracy for the baroclinic update.
# * The external pressure gradient in the momentum equation is multiplied with
#   ho (1/2 before the start of the tracer timestep) in 2D but 0.5(ho+hn)
#   (centered at the start of the tracer time step) in 3D. 2D could potentially
#   reconstruct an equivalent from HU+0.5(zT(i) + zT(i+1)), but that is not fully
#   equivalent as it does not incorporate the future thickness hn. The long-term
#   solution may be to explicitly calculate layer thicknesses and water depth on U
#   and V grids at the start of the tracer timestep. These would be used in the
#   pressure gradient term, and from those (and old thicknesses) the values at t=-1/2
#   could be calculated. *Workaround* (in code): change the 3D external pressure
#   gradient term in momentum_3d.F90 to use ho instead of 0.5(ho+hn). This allows
#   convergence but likely reduces the accuracy of the baroclinic update.
# * For advection of velcoity, both 2D and 3D should use the same method to
#   interpolate velocities to the T and X points [advection grids]. Currently that
#   is not the case: 2D interpolates transports and then divides by thickness, 3D
#   interpolates velcoities directly
# * Freshwater fluxes break convergence because their update is split off in 3D,
#   but done as part of the sealevel update in 2D. *Workaround*: disable rivers,
#   precipitation, evaporation
# * delay_slow_ip must be off in order for 2D and 3D internal pressure terms to
#   converge

import sys
import os.path
import unittest
import datetime

import cftime
import numpy as np

import pygetm
import pygetm.util.compare_nc

sys.path.append(os.path.join(os.path.dirname(__file__), "../examples"))

import north_sea


class TestConvergence(unittest.TestCase):
    def setUp(self) -> None:
        setups_dir = "../../../getm-setups"
        if "GETM_SETUPS_DIR" in os.environ:
            setups_dir = os.environ["GETM_SETUPS_DIR"]
        self.setup_dir = os.path.join(setups_dir, "NorthSea")
        self.domain = north_sea.create_domain(
            self.setup_dir,
            nlev=1,
            use_rivers=False,
            use_boundaries=False,
            logger=pygetm.parallel.get_logger(level="ERROR"),
        )

    def test(self):
        start = cftime.datetime(2006, 1, 2)
        stop = cftime.datetime(2006, 1, 3)

        sim = north_sea.create_simulation(
            self.domain,
            pygetm.BAROCLINIC,
            self.setup_dir,
            delay_slow_ip=False,
            momentum=pygetm.momentum.Momentum(
                advection_split_2d=pygetm.operators.AdvectionSplit.HALF_ALWAYS
            ),
        )
        output = sim.output_manager.add_netcdf_file(
            "res.nc", interval=datetime.timedelta(hours=1)
        )
        output.request("sst")
        sim.start(start, timestep=60.0, split_factor=1, report=60, profile="ns")

        def print_range(variable, name=None):
            print(
                "  %s: %.3e - %.3e"
                % (name or variable.name, variable.ma.min(), variable.ma.max())
            )

        def print_difference(variable1, variable2, name=None):
            delta = variable1.ma - variable2.ma
            reldelta = delta / variable1.ma
            print(
                "  %s: %.3e - %.3e (relative: %.3e - %.3e)"
                % (
                    name or variable1.name + "-" + variable2.name,
                    delta.min(),
                    delta.max(),
                    reldelta.min(),
                    reldelta.max(),
                )
            )

        while sim.time < stop:
            sim.advance()
        sim.finish()

        TOLERANCE = 1e-14
        self.assertLess(np.abs(sim.momentum.SxA.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SyA.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SxD.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SyD.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SxF.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SyF.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SxB.ma - sim.idpdx).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.SyB.ma - sim.idpdy).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.udev.ma).max(), TOLERANCE)
        self.assertLess(np.abs(sim.momentum.vdev.ma).max(), TOLERANCE)


if __name__ == "__main__":
    unittest.main()

