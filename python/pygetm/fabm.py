import logging
from typing import Any, Optional, Union
from collections.abc import Iterable, Mapping, MutableMapping
import os
import contextlib

import numpy as np
import cftime
import pyfabm

from . import core
from . import tracer
from . import _pygetm
from .constants import TimeVarying


class FABM:
    """Interface to the Framework for Aquatic Biogeochemical Models (FABM)."""

    def __init__(
        self,
        path: Union[os.PathLike[str], str, Mapping[str, Any]] = "fabm.yaml",
        repair: bool = True,
        bioshade_feedback: bool = False,
        libname: str = os.path.join(os.path.dirname(__file__), "fabm"),
        time_varying: TimeVarying = TimeVarying.MACRO,
        squeeze: bool = False,
    ):
        """Initialize FABM interface.

        Args:
            path: Path to FABM configuration file, in YAML format, or a mapping
                representing the parsed content of such a file.
            repair: Whether to attempt to repair invalid FABM state values
                by clipping them to valid ranges.
            bioshade_feedback: Whether to obtain the attenuation coefficient
                for photosynthetic active radiation from FABM. This will be
                available as attribute :attr:`kc`. It can be used to contribute
                to light absorption, e.g., by providing it as argument `kc2_add`
                to :class:`pygetm.radiation.TwoBand`. There it in turn affects
                the heat distribution in the water column.
            libname: Path to FABM shared library, excluding any platform-specific
                prefix/suffix
            time_varying: Model step at which FABM variables will be updated.
                FABM arrays will be flagged with this attribute to tell the
                output manager when they are updated.
            squeeze: Whether FABM's internal representation of arrays has all
                singleton dimensions (with length 1) squeezed out.
        """
        if not isinstance(path, Mapping):
            path = str(path)
        self.path = path
        self.repair = repair
        self.bioshade_feedback = bioshade_feedback
        self.libname = libname
        self.time_varying = time_varying
        self.squeeze = squeeze

        self._variable2array: MutableMapping[pyfabm.Variable, core.Array] = {}
        self._yearday: Optional[pyfabm.Dependency] = None
        self._nyear: Optional[pyfabm.Dependency] = None
        self._yearstart: Optional[cftime.datetime] = None

    def initialize(
        self,
        grid: core.Grid,
        tracer_collection: tracer.TracerCollection,
        tracer_totals: list[tracer.TracerTotal],
        logger: logging.Logger,
    ):
        self.grid = grid
        pyfabm.logger = logger

        # Assign FABM standard names to grid arrays
        grid.hn.fabm_standard_name = "cell_thickness"
        grid.H.fabm_standard_name = "bottom_depth_below_geoid"
        grid.Dclip.fabm_standard_name = "bottom_depth"
        if grid.lon is not None:
            grid.lon.fabm_standard_name = "longitude"
        if grid.lat is not None:
            grid.lat.fabm_standard_name = "latitude"

        def variable_to_array(
            variable: pyfabm.Variable, data: Optional[np.ndarray] = None, **kwargs
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
            if data is not None:
                ar.wrap_ndarray(data, register=False)
            ar.register()
            self._variable2array[variable] = ar

        shape = grid.hn.all_values.shape  # shape including halos
        fabm_shape = shape
        halox = grid.halox
        haloy = grid.haloy
        domain_start = (0, haloy, halox)
        domain_stop = (shape[0], shape[1] - haloy, shape[2] - halox)
        if self.squeeze:
            # Squeeze out singleton dimensions (lengh 1 AND no halos)
            # If there are any such dimensions, the domain rank for FABM
            # becomes 0D, 1D or 2D rather than 3D
            keep = [n != 1 for n in grid.hn.all_values.shape]
            fabm_shape = tuple(n for (n, k) in zip(shape, keep) if k)
            domain_start = tuple(n for (n, k) in zip(domain_start, keep) if k)
            domain_stop = tuple(n for (n, k) in zip(domain_stop, keep) if k)
        model = self.model = pyfabm.Model(
            self.path,
            shape=fabm_shape,
            libname=self.libname,
            start=domain_start,
            stop=domain_stop,
        )

        # State variables
        self.sources_interior = np.zeros_like(model.interior_state)
        self.sources_surface = np.zeros_like(model.surface_state)
        self.sources_bottom = np.zeros_like(model.bottom_state)
        self.vertical_velocity = np.zeros_like(model.interior_state)
        for i, variable in enumerate(model.interior_state_variables):
            ar_w = core.Array(grid=grid)
            ar_w.wrap_ndarray(self.vertical_velocity[i, ...].reshape(shape))
            ar = tracer_collection.add(
                data=variable.value.reshape(shape),
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
            data = variable.value.reshape(shape[1:])
            variable_to_array(variable, data=data, attrs=dict(_part_of_state=True))

        # Add diagnostics, initially without associated data
        # Data will be sent later, only if the variable is selected for output,
        # and thus, activated in FABM
        for variable in model.interior_diagnostic_variables:
            variable_to_array(variable, shape=grid.hn.shape)
        for variable in model.horizontal_diagnostic_variables:
            variable_to_array(variable, shape=grid.H.shape)

        # Required inputs: mask and cell thickness
        if model.fabm.mask_type != 0:
            if hasattr(grid, "mask3d"):
                model.link_mask(
                    grid.mask3d.all_values.reshape(model.interior_domain_shape),
                    grid.mask.all_values.reshape(model.horizontal_domain_shape),
                )
            else:
                model.link_mask(
                    grid.mask.all_values.reshape(model.horizontal_domain_shape)
                )
        if hasattr(grid, "bottom_indices"):
            bottom_indices = grid.bottom_indices.all_values
            model.link_bottom_index(
                bottom_indices.reshape(model.horizontal_domain_shape)
            )
        model.link_cell_thickness(
            grid.hn.all_values.reshape(model.interior_domain_shape)
        )

        # Conserved quantities (depth-integrated)
        self.conserved_quantity_totals = np.empty(
            (len(model.conserved_quantities),) + model.horizontal_domain_shape,
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
            data = self.conserved_quantity_totals[i, ...]
            ar.wrap_ndarray(data.reshape(shape[1:]), register=False)
            tracer_totals.append(tracer.TracerTotal(ar))

        # Optionally request PAR attenuation coefficient from FABM for
        # feedbacks to physics
        self.kc = None
        self._kc_variable = None
        if self.bioshade_feedback:
            self._kc_variable = self.model.find_standard_variable(
                "attenuation_coefficient_of_photosynthetic_radiative_flux"
            )
            if self._kc_variable is not None:
                model.require_data(self._kc_variable)
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

    def start(self, timestep: float, time: Optional[cftime.datetime] = None):
        """Prepare FABM. This includes flagging which diagnostics need saving based on
        the output manager configuration, offering fields registered with the field
        manager to FABM if they have a standard name assigned, and subsequently
        verifying whether FABM has all its dependencies fulfilled.
        """
        # Tell FABM which diagnostics are saved. FABM will allocate and manage memory
        # only for those that are. This MUST be done before calling self.model.start
        for variable in self.model.diagnostic_variables:
            if variable in self._variable2array:
                variable.save = self._variable2array[variable].saved

        # Transfer GETM fields with a standard name to FABM
        for field in self.grid.fields.values():
            for standard_name in field.attrs.get("_fabm_standard_names", []):
                try:
                    variable = self.model.dependencies.find(standard_name)
                except KeyError:
                    continue
                if variable.value is None:
                    field.saved = True
                    if field.z:
                        shape = self.model.interior_domain_shape
                    else:
                        shape = self.model.horizontal_domain_shape
                    variable.link(field.all_values.reshape(shape))

        with contextlib.suppress(KeyError):
            self._yearday = self.model.dependencies[
                "number_of_days_since_start_of_the_year"
            ]

        with contextlib.suppress(KeyError):
            self._nyear = self.model.dependencies["number_of_days_in_year"]

        with contextlib.suppress(KeyError):
            self.model.dependencies["maximum_time_step"].value = timestep

        if self._yearday or self._nyear:
            self._yearstart = cftime.datetime(time.year, 1, 1, calendar=time.calendar)
        if self._yearday:
            self._yearday.value = (time - self._yearstart).total_seconds() / 86400.0
        if self._nyear:
            yearstop = cftime.datetime(time.year + 1, 1, 1, calendar=time.calendar)
            self._nyear.value = (yearstop - self._yearstart).total_seconds() / 86400.0

        # Start FABM. This verifies whether all dependencies are fulfilled and freezes
        # the set of diagnostics that will be saved.
        if not self.model.start():
            raise Exception(
                "FABM failed to start. Likely its configuration is incomplete."
            )

        # Fill GETM placeholder arrays for all FABM diagnostics that will be
        # computed/saved.
        def get_diagnostic_values(diagnostic_variables, shape):
            for variable in diagnostic_variables:
                array = self._variable2array[variable]
                if array.saved:
                    # Provide the array with data (NB it has been registered before)
                    array.wrap_ndarray(variable.value.reshape(shape), register=False)
                else:
                    # Remove the array from the list of available fields
                    del self.grid.fields[array.name]

        get_diagnostic_values(
            self.model.interior_diagnostic_variables, self.grid.hn.all_values.shape
        )
        get_diagnostic_values(
            self.model.horizontal_diagnostic_variables, self.grid.H.all_values.shape
        )

        # Apply mask to all state variables (interior, bottom, surface)
        for variable in self.model.state_variables:
            array = self._variable2array[variable]
            array.all_values[array.all_mask] = variable.missing_value

        if self._kc_variable is not None:
            data = self._kc_variable.value
            data = data.reshape(self.grid.hn.all_values.shape)
            self.kc.wrap_ndarray(data, register=False)

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
        time/space-varying input with :meth:`pygetm.core.Array.set`.

        Args:
            name: name of the dependency
        """
        variable = self.model.dependencies.find(name)
        if array is None:
            shape = self.grid.hn.shape
            shape_ = self.grid.hn.all_values.shape
            if variable in self.model.horizontal_dependencies:
                shape = shape[1:]
                shape_ = shape_[1:]
            elif variable in self.model.scalar_dependencies:
                shape = shape_ = ()
            array = core.Array(
                name=variable.output_name,
                units=variable.units,
                long_name=variable.long_path,
                shape=shape,
                dtype=self.model.fabm.numpy_dtype,
                grid=self.grid,
                attrs=dict(_time_varying=self.time_varying),
            )
            data = np.empty(shape_, dtype=self.model.fabm.numpy_dtype)
            array.wrap_ndarray(data, register=False)
            array.register()
        data = array.all_values.view()
        data.shape = variable.shape
        variable.link(data)
        return array

    def update_sources(
        self, seconds_passed: float, time: Optional[cftime.datetime] = None
    ):
        """Update sources, vertical velocities, and diagnostics.
        This does not update the state variables themselves; that is done by
        :meth:`advance`
        """
        if self._yearstart and self._yearstart.year != time.year:
            # Year has changed
            self._yearstart = cftime.datetime(time.year, 1, 1, calendar=time.calendar)
            if self._nyear:
                yearstop = cftime.datetime(time.year + 1, 1, 1, calendar=time.calendar)
                timedelta = yearstop - self._yearstart
                self._nyear.value = timedelta.total_seconds() / 86400.0
        if self._yearday:
            self._yearday.value = (time - self._yearstart).total_seconds() / 86400.0
        valid = self.model.check_state(self.repair)
        if not (valid or self.repair):
            raise Exception("FABM state contains invalid values.")
        self.model.get_sources(
            seconds_passed,
            out=(self.sources_interior, self.sources_surface, self.sources_bottom),
        )
        self.model.get_vertical_movement(self.vertical_velocity)

    def add_vertical_movement_to_sources(self):
        if self.grid.nz == 1:
            return
        h = self.grid.hn.all_values
        halox = self.grid.halox
        haloy = self.grid.haloy
        mask = self.grid.mask3d.all_values
        for itracer in range(self.model.interior_state.shape[0]):
            w = self.vertical_velocity[itracer, ...].reshape(mask.shape)
            c = self.model.interior_state[itracer, ...].reshape(mask.shape)
            s = self.sources_interior[itracer, ...].reshape(mask.shape)
            _pygetm.vertical_advection_to_sources(halox, haloy, mask, c, w, h, s)

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
