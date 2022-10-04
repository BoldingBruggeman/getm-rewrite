import numbers
from typing import Optional, Union, Tuple, Literal, Mapping, Any, TYPE_CHECKING
import logging
from functools import partial

import numpy as np
import numpy.lib.mixins
from numpy.typing import DTypeLike, ArrayLike
import xarray

from . import _pygetm
from . import parallel
from .constants import CENTERS, INTERFACES, SPONGE

if TYPE_CHECKING:
    from . import domain


def _noop(*args, **kwargs):
    pass


class Array(_pygetm.Array, numpy.lib.mixins.NDArrayOperatorsMixin):
    __slots__ = (
        "_xarray",
        "_scatter",
        "_gather",
        "_name",
        "attrs",
        "_fill_value",
        "_ma",
        "saved",
        "_shape",
        "_ndim",
        "_size",
        "_dtype",
        "update_halos",
        "update_halos_start",
        "update_halos_finish",
        "compare_halos",
    )
    grid: "domain.Grid"

    def __init__(
        self,
        name: Optional[str] = None,
        units: Optional[str] = None,
        long_name: Optional[str] = None,
        fill_value: Optional[Union[float, int]] = None,
        shape: Optional[Tuple[int]] = None,
        dtype: Optional[DTypeLike] = None,
        grid: "domain.Grid" = None,
        fabm_standard_name: Optional[str] = None,
        attrs: Mapping[str, Any] = {},
    ):
        _pygetm.Array.__init__(self, grid)
        self._xarray: Optional[xarray.DataArray] = None
        self._scatter: Optional[parallel.Scatter] = None
        self._gather: Optional[parallel.Gather] = None
        assert (
            fill_value is None or np.ndim(fill_value) == 0
        ), "fill_value must be a scalar value"
        self._name = name
        self.attrs = attrs.copy()
        if units:
            self.attrs["units"] = units
        if long_name:
            self.attrs["long_name"] = long_name
        if fabm_standard_name:
            self.attrs.setdefault("_fabm_standard_names", set()).add(fabm_standard_name)
        self._fill_value = (
            fill_value
            if fill_value is None or dtype is None
            else np.array(fill_value, dtype=dtype)
        )
        self._ma = None
        self.saved = False  #: to be set if this variable is requested for output
        self._shape = shape
        self._ndim = None if shape is None else len(shape)
        self._size = None if shape is None else np.prod(shape)
        self._dtype = dtype
        self.values = None

    def set_fabm_standard_name(self, fabm_standard_name):
        self.attrs.setdefault("_fabm_standard_names", set()).add(fabm_standard_name)

    fabm_standard_name = property(fset=set_fabm_standard_name)

    def finish_initialization(self):
        """This is called by the underlying cython implementation after the array
        receives a value (:attr:`all_values` is valid)
        """
        assert self.grid is not None
        self._dtype = self.all_values.dtype
        self._ndim = self.all_values.ndim
        if self._fill_value is not None:
            # Cast fill value to dtype of the array
            self._fill_value = np.array(self._fill_value, dtype=self._dtype)
        if self.on_boundary or self._ndim == 0:
            # boundary array or scalar
            self.values = self.all_values[...]
        else:
            halox, haloy = self.grid.domain.halox, self.grid.domain.haloy
            self.values = self.all_values[..., haloy:-haloy, halox:-halox]

        self._shape = self.values.shape
        self._size = self.values.size

        if not self.grid.domain.tiling:
            self.update_halos = _noop
            self.update_halos_start = _noop
            self.update_halos_finish = _noop
            self.compare_halos = _noop
        else:
            self.update_halos = partial(self._distribute, "update_halos")
            self.update_halos_start = partial(self._distribute, "update_halos_start")
            self.update_halos_finish = partial(self._distribute, "update_halos_finish")
            self.compare_halos = partial(self._distribute, "compare_halos")

    def register(self):
        assert self.grid is not None
        if self._name is not None:
            if self._name in self.grid.domain.fields:
                raise Exception(
                    "A field with name '%s' has already been registered"
                    " with the field manager." % self._name
                )
            self.grid.domain.fields[self._name] = self

    def __repr__(self) -> str:
        return super().__repr__() + self.grid.postfix

    def _distribute(self, method: str, *args, **kwargs) -> parallel.DistributedArray:
        dist = parallel.DistributedArray(
            self.grid.domain.tiling,
            self.all_values,
            self.grid.halo,
            overlap=self.grid.overlap,
        )
        self.update_halos = dist.update_halos
        self.update_halos_start = dist.update_halos_start
        self.update_halos_finish = dist.update_halos_finish
        self.compare_halos = dist.compare_halos
        getattr(self, method)(*args, **kwargs)

    def scatter(self, global_data: Optional["Array"]):
        if self.grid.domain.tiling.n == 1:
            self.values[...] = global_data
            return
        if self._scatter is None:
            self._scatter = parallel.Scatter(
                self.grid.domain.tiling,
                self.all_values,
                halo=self.grid.halo,
                fill_value=self._fill_value,
            )
        self._scatter(None if global_data is None else global_data.all_values)

    def gather(self, out: Optional["Array"] = None, slice_spec=()):
        if self.grid.domain.tiling.n == 1:
            if out is not None:
                out[slice_spec + (Ellipsis,)] = self.values
            return self
        if self._gather is None:
            self._gather = parallel.Gather(
                self.grid.domain.tiling,
                self.values.shape,
                self.dtype,
                fill_value=self._fill_value,
            )
        result = self._gather(
            self.values,
            out.values if isinstance(out, Array) else out,
            slice_spec=slice_spec,
        )
        if result is not None and out is None:
            out = self.grid.domain.glob.grids[self.grid.type].array(dtype=self.dtype)
            out[...] = result
        return out

    def allgather(self) -> np.ndarray:
        if self.grid.domain.tiling.n == 1:
            return self.values
        if self.on_boundary:
            comm = self.grid.domain.tiling.comm
            open_boundaries = self.grid.domain.open_boundaries

            # Gather the number of open boundary points in each subdomain
            np_bdy = np.empty((comm.size,), dtype=int)
            np_local = np.array(open_boundaries.np, dtype=int)
            comm.Allgather(np_local, np_bdy)

            # Gather the global indices of open boundary points from each subdomain
            indices = np.empty((np_bdy.sum(),), dtype=int)
            comm.Allgatherv(open_boundaries.local_to_global_indices, (indices, np_bdy))
            assert frozenset(indices) == frozenset(range(open_boundaries.np_glob))

            # Gather the values at the open boundary points from each subdomain
            values = np.empty((np_bdy.sum(),), dtype=self.values.dtype)
            comm.Allgatherv(self.values, (values, np_bdy))

            # Map retrieved values to the appropriate indices in the global array
            all_values = np.empty((open_boundaries.np_glob,), dtype=values.dtype)
            all_values[indices] = values
            return all_values

    def global_sum(
        self, reproducible: bool = False, where: Optional["Array"] = None
    ) -> Optional[np.ndarray]:
        if reproducible:
            all = self.gather()
            if where is not None:
                where = where.gather()
            if all is not None:
                return all.values.sum(
                    where=np._NoValue if where is None else where.values
                )
        else:
            local_sum = self.values.sum(
                where=np._NoValue if where is None else where.values
            )
            return parallel.Sum(self.grid.domain.tiling, local_sum)()

    def global_mean(
        self, reproducible: bool = False, where: Optional["Array"] = None
    ) -> Optional[np.ndarray]:
        sum = self.global_sum(reproducible=reproducible, where=where)
        if where is not None:
            count = where.global_sum()
        else:
            count = parallel.Sum(self.grid.domain.tiling, self.values.size)()
        if sum is not None:
            return sum / count

    @staticmethod
    def create(
        grid: "domain.Grid",
        fill: Optional[ArrayLike] = None,
        z: Literal[None, True, False, CENTERS, INTERFACES] = None,
        dtype: DTypeLike = None,
        copy: bool = True,
        on_boundary: bool = False,
        register: bool = True,
        **kwargs
    ) -> "Array":
        """Create a new :class:`Array`

        Args:
            grid: grid associated with the new array
            fill: value to set the new array to
            z: vertical dimension of the new array.
                ``False`` for a 2D array, ``CENTERS`` (or ``True``) for an array
                defined at the layer centers, ``INTERFACES`` for an array defined at
                the layer interfaces. ``None`` to detect from ``fill``.
            dtype: data type
            copy: whether to create a copy of ``fill``, if provided
            on_boundary: whether to describe data along the open boundaries (1D),
                instead of the 2D x-y model domain
            register: whether to register the array as field available for output
            **kwargs: additional keyword arguments passed to :class:`Array`
        """
        ar = Array(grid=grid, **kwargs)
        if fill is None and ar.fill_value is not None:
            fill = ar.fill_value
        if fill is not None:
            fill = np.asarray(fill)
            if z is None and not on_boundary:
                if fill.ndim != 3:
                    z = False
                elif fill.shape[0] == grid.nz_ + 1:
                    z = INTERFACES
                else:
                    z = CENTERS
        if dtype is None:
            dtype = float if fill is None else fill.dtype
        shape = (
            [grid.domain.open_boundaries.np] if on_boundary else [grid.ny_, grid.nx_]
        )
        if z:
            shape.insert(
                1 if on_boundary else 0, grid.nz_ + 1 if z == INTERFACES else grid.nz_,
            )
        if copy or fill is None:
            data = np.empty(shape, dtype=dtype)
            if fill is not None:
                data[...] = fill
        else:
            data = np.broadcast_to(fill, shape)
        ar.wrap_ndarray(data, on_boundary=on_boundary, register=register)
        return ar

    def fill(self, value):
        """Set array to specified value, while respecting the mask: masked points are
        set to :attr:`fill_value`
        """
        try:
            self.all_values[...] = value
        except ValueError:
            self.values[...] = value
            self.update_halos()
        if self.fill_value is not None and not (self.ndim == 0 or self.on_boundary):
            self.all_values[..., self.grid._land] = self.fill_value

    @property
    def ma(self) -> np.ma.MaskedArray:
        """Masked array representation that combines the data and the mask associated
        with the array's native grid
        """
        if self._ma is None:
            if self.size == 0 or self.on_boundary:
                mask = False
            else:
                mask = self.grid.mask.values == 0
            self._ma = np.ma.array(self.values, mask=np.broadcast_to(mask, self._shape))
        return self._ma

    def plot(self, mask: bool = True, **kwargs):
        """Plot the array with :meth:`xarray.DataArray.plot`

        Args:
            **kwargs: additional keyword arguments passed to
                :meth:`xarray.DataArray.plot`
        """
        if "x" not in kwargs and "y" not in kwargs:
            kwargs["x"] = (
                "lon" if self.grid.domain.spherical else "x"
            ) + self.grid.postfix
            kwargs["y"] = (
                "lat" if self.grid.domain.spherical else "y"
            ) + self.grid.postfix
        if "shading" not in kwargs:
            kwargs["shading"] = "auto"
        return self.as_xarray(mask=mask).plot(**kwargs)

    def interp(
        self,
        target: Union["Array", "domain.Grid"],
        z: Literal[None, True, False, CENTERS, INTERFACES] = None,
    ) -> "Array":
        """Interpolate the array to another grid.

        Args:
            target: either the :class:`Array` that will hold the interpolated data,
                or the :class:`~pygetm.domain.Grid` to interpolate to. If a ``Grid`` is
                provided, a new array will be created to hold the interpolated values.
        """
        if not isinstance(target, Array):
            # Target must be a grid; we need to create the array
            target_z = z if z is not None else self.z
            target = Array.create(target, dtype=self._dtype, z=target_z)
        source_array, target_array = self.all_values, target.all_values
        if self.grid is target.grid:
            if self.z == INTERFACES and target.z == CENTERS:
                # vertical interpolation from layer interfaces to layer centers
                _pygetm.interp_z(source_array, target_array, offset=0)
            elif self.z == CENTERS and target.z == INTERFACES:
                # vertical interpolation from layer centers to layer interfaces
                # (top and bottom interfaces will be left untouched)
                _pygetm.interp_z(source_array, target_array, offset=1)
        else:
            if self._ndim == 2:
                source_array = source_array[None, ...]
                target_array = target_array[None, ...]
            interpolate = self.grid.interpolator(target.grid)
            interpolate(source_array, target_array)
        return target

    def __array__(self, dtype: Optional[DTypeLike] = None) -> np.ndarray:
        """Return interior of the array as a NumPy array.
        No copy will be made unless the requested data type differs from that
        of the underlying array.

        Args:
            dtype: data type
        """
        return np.asarray(self.values, dtype=dtype)

    def isel(self, *, z: int, **kwargs) -> "Array":
        """Select a single depth level. The data in the returned 2D :class:`Array`
        will be a view of the relevant data of the original 3D array. Thus, changes
        in one will affect the other.
        """
        if self._ndim != 3:
            raise NotImplementedError
        if self.units is not None:
            kwargs.setdefault("units", self.units)
        if self.long_name is not None:
            kwargs.setdefault("long_name", "%s @ k=%i" % (self.long_name, z))
        kwargs["attrs"] = kwargs.get("attrs", {}).copy()
        for att in ("_mask_output",):
            if att in self.attrs:
                kwargs["attrs"][att] = self.attrs[att]
        ar = Array(grid=self.grid, fill_value=self.fill_value, **kwargs)
        ar.wrap_ndarray(self.all_values[z, ...])
        return ar

    def __getitem__(self, key) -> np.ndarray:
        """Retrieve values from the interior of the array (excluding halos).
        For access to the halos, use :attr:`all_values`.
        """
        return self.values[key]

    def __setitem__(self, key, values):
        """Assign values to the interior of the array (excluding halos).
        For access to the halos, use :attr:`all_values`.
        """
        self.values[key] = values

    @property
    def shape(self) -> Tuple[int]:
        """Shape excluding halos"""
        return self._shape

    @property
    def ndim(self) -> int:
        """Number of dimensions"""
        return self._ndim

    @property
    def z(self):
        """Vertical dimension: ``False`` if the array has no vertical dimension,
        ``CENTERS`` for layer centers, ``INTERFACES`` for layer interfaces.
        """
        if self._ndim != (2 if self.on_boundary else 3):
            return False
        nz = self._shape[1 if self.on_boundary else 0]
        return INTERFACES if nz == self.grid.nz_ + 1 else CENTERS

    @property
    def size(self) -> int:
        """Total number of values, excluding halos"""
        return self._size

    @property
    def dtype(self) -> DTypeLike:
        """Data type"""
        return self._dtype

    @property
    def name(self) -> Optional[str]:
        """Name"""
        return self._name

    @property
    def units(self) -> Optional[str]:
        """Units"""
        return self.attrs.get("units")

    @property
    def long_name(self) -> Optional[str]:
        """Long name"""
        return self.attrs.get("long_name")

    @property
    def fill_value(self) -> Optional[Union[int, float]]:
        """Fill value"""
        return self._fill_value

    # Below based on https://np.org/devdocs/reference/generated/np.lib.mixins.NDArrayOperatorsMixin.html#np.lib.mixins.NDArrayOperatorsMixin
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__":
            return NotImplemented

        out = kwargs.get("out", ())
        for x in inputs + out:
            # Only support operations with instances of _HANDLED_TYPES.
            # Use ArrayLike instead of type(self) for isinstance to
            # allow subclasses that don't override __array_ufunc__ to
            # handle ArrayLike objects.
            if not isinstance(x, (np.ndarray, numbers.Number, Array)):
                return NotImplemented
            if isinstance(x, Array) and x.grid is not self.grid:
                return NotImplemented

        # Defer to the implementation of the ufunc on unwrapped values.
        inputs = tuple(x.all_values if isinstance(x, Array) else x for x in inputs)
        if out:
            kwargs["out"] = tuple(
                x.all_values if isinstance(x, Array) else x for x in out
            )
        result = getattr(ufunc, method)(*inputs, **kwargs)

        if type(result) is tuple:
            # multiple return values
            return tuple(self.create(self.grid, x) for x in result)
        elif method == "at":
            # no return value
            return None
        else:
            # one return value
            return self.create(self.grid, result)

    def set(self, value: Union[float, np.ndarray, xarray.DataArray], **kwargs):
        """Link this array to a field or value using
        :attr:`~pygetm.domain.Domain.input_manager`, which will perform temporal and
        spatial interpolation as required.

        Args:
            value: value to assign to this array. If it is time-dependent (if you pass
                an instance of :class:`xarray.DataArray` with a time dimension),
                the array's value will be updated during the simulation whenever
                :meth:`pygetm.input.InputManager.update` is called.
            **kwargs: keyword arguments passed to :meth:`pygetm.input.InputManager.add`
        """
        self.grid.domain.input_manager.add(self, value, **kwargs)

    def require_set(self, logger: Optional[logging.Logger] = None):
        """Assess whether all non-masked cells of this field have been set. If not, an
        error message is written to the log and False is returned.
        """
        valid = True
        if self._fill_value is not None:
            invalid = self.ma == self._fill_value
            if invalid.any():
                (logger or logging.getLogger()).error(
                    "%s is masked (%s) in %i active grid cells."
                    % (self.name, self._fill_value, invalid.sum())
                )
                valid = False
        return valid

    def as_xarray(self, mask: bool = False) -> xarray.DataArray:
        """Return this array wrapped in an :class:`xarray.DataArray` that includes
        coordinates and can be used for plotting
        """
        if self._xarray is not None and not mask:
            return self._xarray
        attrs = {}
        for key in ("units", "long_name"):
            value = getattr(self, key)
            if value is not None:
                attrs[key] = value
        dom = self.grid.domain
        coords = {}
        if self.name not in (
            "x" + self.grid.postfix,
            "y" + self.grid.postfix,
            "lon" + self.grid.postfix,
            "lat" + self.grid.postfix,
        ):
            if dom.x_is_1d:
                coords["x%s" % self.grid.postfix] = self.grid.x.xarray[0, :]
            if dom.y_is_1d:
                coords["y%s" % self.grid.postfix] = self.grid.y.xarray[:, 0]
            coords["x%s2" % self.grid.postfix] = self.grid.x.xarray
            coords["y%s2" % self.grid.postfix] = self.grid.y.xarray
            if dom.lon is not None:
                coords["lon%s" % self.grid.postfix] = self.grid.lon.xarray
            if dom.lat is not None:
                coords["lat%s" % self.grid.postfix] = self.grid.lat.xarray
        dims = ("y" + self.grid.postfix, "x" + self.grid.postfix)
        if self.ndim == 3:
            dims = ("zi" if self.z == INTERFACES else "z",) + dims
        values = self.values if not mask else self.ma
        _xarray = xarray.DataArray(
            values, coords=coords, dims=dims, attrs=attrs, name=self.name
        )
        if not mask:
            self._xarray = _xarray
        return _xarray

    xarray = property(as_xarray)

