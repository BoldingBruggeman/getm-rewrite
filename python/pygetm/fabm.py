import logging
from typing import List, Optional, MutableMapping, Iterable

import numpy as np
import cftime

from . import domain
from . import pyfabm
from . import core
from . import tracer
from . import _pygetm
from .constants import TimeVarying


class FABM:
    def __init__(
        self,
        path: str = "fabm.yaml",
        repair: bool = True,
        bioshade_feedback: bool = False,
        libname: str = "fabm",
        time_varying: TimeVarying = TimeVarying.MACRO,
    ):
        self.path = path
        self.repair = repair
        self.bioshade_feedback: bool = bioshade_feedback
        self.libname = libname
        self.time_varying = time_varying

        self._variable2array: MutableMapping[pyfabm.Variable, core.Array] = {}

    def initialize(
        self,
        grid: domain.Grid,
        tracer_collection: tracer.TracerCollection,
        tracer_totals: List[tracer.TracerTotal],
        logger: logging.Logger,
    ):
        self.grid = grid
        pyfabm.logger = logger

        # Assign FABM standard names to grid arrays
        grid.hn.fabm_standard_name = "cell_thickness"
        grid.H.fabm_standard_name = "bottom_depth_below_geoid"
        grid.D.fabm_standard_name = "bottom_depth"
        if grid.lon is not None:
            grid.lon.fabm_standard_name = "longitude"
        if grid.lat is not None:
            grid.lat.fabm_standard_name = "latitude"

        def variable_to_array(
            variable: pyfabm.Variable, send_data: bool = False, **kwargs
        ):
            kwargs.setdefault("attrs", {})["_time_varying"] = self.time_varying
            ar = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                fill_value=variable.missing_value,
                dtype=self.model.fabm.numpy_dtype,
                grid=grid,
                **kwargs,
            )
            if send_data:
                ar.wrap_ndarray(variable.data, register=False)
            ar.register()
            self._variable2array[variable] = ar

        shape = grid.hn.all_values.shape  # shape including halos
        halo = grid.domain.halo
        model = self.model = pyfabm.Model(
            self.path,
            shape=shape,
            libname=self.libname,
            start=(0, halo, halo),
            stop=(shape[0], shape[1] - halo, shape[2] - halo),
        )

        # State variables
        self.sources_interior = np.zeros_like(model.interior_state)
        self.sources_surface = np.zeros_like(model.surface_state)
        self.sources_bottom = np.zeros_like(model.bottom_state)
        self.vertical_velocity = np.zeros_like(model.interior_state)
        for i, variable in enumerate(model.interior_state_variables):
            ar_w = core.Array(grid=grid)
            ar_w.wrap_ndarray(self.vertical_velocity[i, ...])
            ar = tracer_collection.add(
                data=variable.data,
                vertical_velocity=ar_w,
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                fill_value=variable.missing_value,
                rivers_follow_target_cell=variable.no_river_dilution,
                precipitation_follows_target_cell=variable.no_precipitation_dilution,
            )
            self._variable2array[variable] = ar
        for variable in model.surface_state_variables + model.bottom_state_variables:
            variable_to_array(variable, send_data=True, attrs=dict(_part_of_state=True))

        # Add diagnostics, initially without associated data
        # Data will be sent later, only if the variable is selected for output,
        # and thus, activated in FABM
        for variable in model.diagnostic_variables:
            current_shape = grid.H.shape if variable.horizontal else grid.hn.shape
            variable_to_array(variable, shape=current_shape)

        # Required inputs: mask and cell thickness
        if hasattr(grid.domain, "mask3d"):
            model.link_mask(grid.domain.mask3d.all_values, grid.mask.all_values)
        else:
            model.link_mask(grid.mask.all_values)
        if hasattr(grid.domain, "bottom_indices"):
            model.link_bottom_index(grid.domain.bottom_indices.all_values)
        model.link_cell_thickness(grid.hn.all_values)

        # Conserved quantities (depth-integrated)
        self.conserved_quantity_totals = np.empty(
            (len(model.conserved_quantities),) + shape[1:],
            dtype=self.sources_interior.dtype,
        )
        for i, variable in enumerate(model.conserved_quantities):
            ar = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_name,
                fill_value=variable.missing_value,
                dtype=model.fabm.numpy_dtype,
                grid=grid,
                attrs=dict(_time_varying=self.time_varying),
            )
            ar.wrap_ndarray(self.conserved_quantity_totals[i, ...], register=False)
            tracer_totals.append(tracer.TracerTotal(ar))

        # Optionally request PAR attenuation coefficient from FABM for
        # feedbacks to physics
        self.kc = None
        self.kc_variable = None
        if self.bioshade_feedback:
            self.kc_variable = self.model.find_standard_variable(
                "attenuation_coefficient_of_photosynthetic_radiative_flux"
            )
            if self.kc_variable is not None:
                model.require_data(self.kc_variable)
                self.kc = core.Array(
                    name="kc_fabm",
                    units="m-1",
                    long_name="attenuation of visible radiation by FABM",
                    shape=grid.hn.shape,
                    dtype=self.model.fabm.numpy_dtype,
                    grid=self.grid,
                    attrs=dict(_time_varying=self.time_varying),
                )
                self.kc.register()

    @property
    def default_outputs(self) -> Iterable[core.Array]:
        return [a for v, a in self._variable2array.items() if v.output]

    @property
    def state_variables(self) -> Iterable[core.Array]:
        return [self._variable2array[v] for v in self.model.state_variables]

    def start(self, time: Optional[cftime.datetime] = None):
        """Prepare FABM. This includes flagging which diagnostics need saving based on
        the output manager configuration, offering fields registered with the field
        manager to FABM if they have a standard name assigned, and subsequently
        verifying whether FABM has all its dependencies fulfilled.
        """
        # Tell FABM which diagnostics are saved. FABM will allocate and manage memory
        # only for those that are. This MUST be done before calling self.model.start
        for variable in self.model.diagnostic_variables:
            variable.save = self._variable2array[variable].saved

        # Transfer GETM fields with a standard name to FABM
        for field in self.grid.domain.fields.values():
            for standard_name in field.attrs.get("_fabm_standard_names", []):
                try:
                    variable = self.model.dependencies.find(standard_name)
                except KeyError:
                    continue
                if not variable.is_set:
                    field.saved = True
                    variable.link(field.all_values)

        try:
            self._yearday = self.model.dependencies.find(
                "number_of_days_since_start_of_the_year"
            )
            timedelta = time - cftime.datetime(time.year, 1, 1, calendar=time.calendar)
            self._yearday.value = timedelta.total_seconds() / 86400.0
        except KeyError:
            self._yearday = None

        # Start FABM. This verifies whether all dependencies are fulfilled and freezes
        # the set of diagnostics that will be saved.
        if not self.model.start():
            raise Exception(
                "FABM failed to start. Likely its configuration is incomplete."
            )

        # Fill GETM placeholder arrays for all FABM diagnostics that will be
        # computed/saved.
        for variable in self.model.diagnostic_variables:
            array = self._variable2array[variable]
            if array.saved:
                # Provide the array with data (NB it has been registered before)
                array.wrap_ndarray(variable.data, register=False)
            else:
                # Remove the array from the list of available fields
                del self.grid.domain.fields[array.name]

        # Apply mask to all state variables (interior, bottom, surface)
        for variable in self.model.interior_state_variables:
            mask = getattr(self.grid, "_land3d", self.grid._land)
            variable.value[..., mask] = variable.missing_value
        for variable in self.model.bottom_state_variables:
            variable.value[..., self.grid._land] = variable.missing_value
        for variable in self.model.surface_state_variables:
            variable.value[..., self.grid._land] = variable.missing_value

        if self.kc_variable is not None:
            self.kc.wrap_ndarray(self.kc_variable.value, register=False)

    def has_dependency(self, name: str) -> bool:
        try:
            self.model.dependencies.find(name)
        except KeyError:
            return False
        return True

    def get_dependency(
        self, name: str, array: Optional[core.Array] = None
    ) -> core.Array:
        """Retrieve the array that will hold values for the specified FABM dependency.
        This array can subsequently be assigned a value or be linked to a
        time/space-varying input with :attr:`~pygetm.core.Array.set`.

        Args:
            name: name of the dependency
        """
        variable = self.model.dependencies.find(name)
        if array is None:
            array = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                shape=variable.shape,
                dtype=self.model.fabm.numpy_dtype,
                grid=self.grid,
                attrs=dict(_time_varying=self.time_varying),
            )
            data = np.empty(variable.shape, dtype=self.model.fabm.numpy_dtype)
            array.wrap_ndarray(data, register=False)
        else:
            data = array.all_values.view()
            data.shape = variable.shape
        variable.link(data)
        return array

    def update_sources(self, time: Optional[cftime.datetime] = None):
        """Update sources, vertical velocities, and diagnostics.
        This does not update the state variables themselves; that is done by
        :meth:`advance`
        """
        if self._yearday:
            timedelta = time - cftime.datetime(time.year, 1, 1, calendar=time.calendar)
            self._yearday.value = timedelta.total_seconds() / 86400.0
        valid = self.model.check_state(self.repair)
        if not (valid or self.repair):
            raise Exception("FABM state contains invalid values.")
        self.model.get_sources(
            out=(self.sources_interior, self.sources_surface, self.sources_bottom)
        )
        self.model.get_vertical_movement(self.vertical_velocity)

    def add_vertical_movement_to_sources(self):
        h = self.grid.hn.all_values
        halo = self.grid.domain.halo
        mask = self.grid.domain.mask3d.all_values
        for itracer in range(self.model.interior_state.shape[0]):
            w = self.vertical_velocity[itracer, ...]
            if w.any():
                c = self.model.interior_state[itracer, ...]
                s = self.sources_interior[itracer, ...]
                _pygetm.vertical_advection_to_sources(halo, mask, c, w, h, s)

    def update_totals(self):
        """Ensure sums of conserved quantities are up to date."""
        self.model.get_conserved_quantities(out=self.conserved_quantity_totals)

    def advance(self, timestep: float):
        """Time-integrate source terms of all state variables (3D pelagic tracers as
        well as bottom- and surface-attached variables).

        Args:
            timestep: time step (s)
        """
        _pygetm.multiply_add(
            self.model.interior_state.ravel(), self.sources_interior.ravel(), timestep
        )
        _pygetm.multiply_add(
            self.model.surface_state.ravel(), self.sources_surface.ravel(), timestep
        )
        _pygetm.multiply_add(
            self.model.bottom_state.ravel(), self.sources_bottom.ravel(), timestep
        )
