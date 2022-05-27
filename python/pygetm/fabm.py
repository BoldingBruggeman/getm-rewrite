import itertools
import logging
from typing import List

import numpy
import cftime

from . import domain
from . import pyfabm
from . import core
from . import tracer


class FABM:
    def __init__(self, path: str = "fabm.yaml", repair: bool = True):
        self.path = path
        self.repair = repair

    def initialize(
        self,
        grid: domain.Grid,
        tracer_collection: tracer.TracerCollection,
        tracer_totals: List[core.Array],
        logger: logging.Logger,
    ):
        self.grid = grid
        pyfabm.logger = logger

        def variable_to_array(variable, send_data: bool = False, **kwargs):
            ar = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                fill_value=variable.missing_value,
                dtype=self.model.fabm.dtype,
                grid=grid,
                **kwargs
            )
            if send_data:
                ar.wrap_ndarray(variable.data)
            ar.register()
            return ar

        shape = grid.hn.all_values.shape  # shape including halos
        halo = grid.domain.halo
        self.model = pyfabm.Model(
            self.path,
            shape=shape,
            libname="fabm_c",
            start=(0, halo, halo),
            stop=(shape[0], shape[1] - halo, shape[2] - halo),
        )

        # State variables
        self.sources_interior = numpy.zeros_like(self.model.interior_state)
        self.sources_surface = numpy.zeros_like(self.model.surface_state)
        self.sources_bottom = numpy.zeros_like(self.model.bottom_state)
        self.vertical_velocity = numpy.zeros_like(self.model.interior_state)
        for i, variable in enumerate(self.model.interior_state_variables):
            ar_w = core.Array(grid=grid)
            ar_w.wrap_ndarray(self.vertical_velocity[i, ...])
            tracer_collection.add(
                data=variable.data,
                vertical_velocity=ar_w,
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                fill_value=variable.missing_value,
                rivers_follow_target_cell=variable.no_river_dilution,
                precipitation_follows_target_cell=variable.no_precipitation_dilution,
            )
        for variable in self.model.surface_state_variables:
            variable_to_array(variable, send_data=True, attrs={"_part_of_state": True})
        for variable in self.model.bottom_state_variables:
            variable_to_array(variable, send_data=True, attrs={"_part_of_state": True})

        # Diagnostics
        self._interior_diagnostic_arrays = [
            variable_to_array(variable, shape=grid.hn.shape)
            for variable in self.model.interior_diagnostic_variables
        ]
        self._horizontal_diagnostic_arrays = [
            variable_to_array(variable, shape=grid.H.shape)
            for variable in self.model.horizontal_diagnostic_variables
        ]

        # Required inputs: mask and cell thickness
        self.model.link_mask(grid.mask.all_values)
        self.model.link_cell_thickness(grid.hn.all_values)

        # Conserved quantities
        self.conserved_quantity_totals = numpy.empty(
            (len(self.model.conserved_quantities),) + shape[1:],
            dtype=self.sources_interior.dtype,
        )
        for i, variable in enumerate(self.model.conserved_quantities):
            ar = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_name,
                fill_value=variable.missing_value,
                dtype=self.model.fabm.dtype,
                grid=grid,
            )
            ar.wrap_ndarray(self.conserved_quantity_totals[i, ...])
            tracer_totals.append(ar)

    def start(self, time: cftime.datetime):
        """Prepare FABM. This includes flagging which diagnostics need saving based on
        the output manager configuration, offering fields registered with the field
        manager to FABM if they have a standard name assigned, and subsequently
        verifying whether FABM has all its dependencies fulfilled.
        """
        # Tell FABM which diagnostics are saved. FABM will allocate and manage memory
        # only for those that are. This MUST be done before calling self.model.start
        for variable, ar in zip(
            itertools.chain(
                self.model.interior_diagnostic_variables,
                self.model.horizontal_diagnostic_variables,
            ),
            itertools.chain(
                self._interior_diagnostic_arrays, self._horizontal_diagnostic_arrays
            ),
        ):
            variable.save = ar.saved

        # Transfer GETM fields with a standard name to FABM
        for field in self.grid.domain.field_manager.fields.values():
            for standard_name in field.attrs.get("_fabm_standard_names", []):
                try:
                    variable = self.model.dependencies.find(standard_name)
                except KeyError:
                    continue
                field.saved = True
                variable.link(field.all_values)

        try:
            self._yearday = self.model.dependencies.find(
                "number_of_days_since_start_of_the_year"
            )
            self._yearday.value = (
                time - cftime.datetime(time.year, 1, 1)
            ).total_seconds() / 86400.0
        except KeyError:
            self._yearday = None

        # Start FABM. This verifies whether all dependencies are fulfilled and freezes
        # the set of diagnostics that will be saved.
        assert (
            self.model.start()
        ), "FABM failed to start. Likely its configuration is incomplete."

        # Fill GETM placeholder arrays for all FABM diagnostics that will be
        # computed/saved.
        for variable, ar in zip(
            itertools.chain(
                self.model.interior_diagnostic_variables,
                self.model.horizontal_diagnostic_variables,
            ),
            itertools.chain(
                self._interior_diagnostic_arrays, self._horizontal_diagnostic_arrays
            ),
        ):
            if ar.saved:
                ar.wrap_ndarray(variable.data)

        # Apply mask to all state variables (interior, bottom, surface)
        for variable in self.model.state_variables:
            variable.value[..., self.grid.mask.all_values == 0] = variable.missing_value

    def get_dependency(self, name: str) -> core.Array:
        """Retrieve the array that will hold values for the specified FABM dependency.
        This array can subsequently be assigned a value or be linked to a
        time/space-varying input with :attr:`~pygetm.core.Array.set`.

        Args:
            name: name of the dependency
        """
        variable = self.model.dependencies.find(name)
        if len(variable.shape) == 0:
            return variable
        arr = self.grid.array(
            name=variable.output_name,
            units=variable.units,
            long_name=variable.long_path,
            z=len(variable.shape) == 3,
        )
        variable.link(arr.all_values)
        return arr

    def update_sources(self, time: cftime.datetime):
        """Update sources, vertical velocities, and diagnostics.
        This does not update the state variables themselves; that is done by
        :meth:`advance`
        """
        if self._yearday:
            self._yearday.value = (
                time - cftime.datetime(time.year, 1, 1)
            ).total_seconds() / 86400.0
        self.model.get_sources(
            out=(self.sources_interior, self.sources_surface, self.sources_bottom)
        )
        self.model.get_vertical_movement(self.vertical_velocity)

    def update_totals(self):
        """Ensure sums of conserved quantities are up to date."""
        self.model.get_conserved_quantities(out=self.conserved_quantity_totals)

    def advance(self, timestep: float):
        """Time-integrate source terms of all state variables (3D pelagic tracers as
        well as bottom- and surface-attached variables).

        Args:
            timestep: time step (s)
        """
        self.sources_interior *= timestep
        self.sources_surface *= timestep
        self.sources_bottom *= timestep
        self.model.interior_state += self.sources_interior
        self.model.surface_state += self.sources_surface
        self.model.bottom_state += self.sources_bottom
        self.model.check_state(repair=self.repair)
