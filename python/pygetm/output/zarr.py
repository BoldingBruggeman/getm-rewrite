from typing import Optional, Mapping
import logging
import asyncio

import cftime
import zarr.api.asynchronous as zarr
from zarr.storage import StoreLike
import numpy as np

from . import File
from . import operators
import pygetm.core


async def _create_arrays(
    root: zarr.AsyncGroup, fields: Mapping[str, pygetm.core.Array], chunk_size: int
) -> tuple[list[zarr.AsyncArray], bool, bool]:
    needs_time = False
    needs_time_av = False
    tasks: list[asyncio.Task] = []
    for output_name, field in fields.items():
        shape = field.shape
        chunks = field.shape
        dims = field.dims
        coords = field.coordinates
        if field.time_varying:
            shape = (0,) + shape
            chunks = (chunk_size,) + chunks
            dims = ("time",) + dims
            if "time: mean" in field.attrs.get("cell_methods", ""):
                coords = coords + ["time_av"]
                needs_time_av = True
            else:
                needs_time = True
        # Variable attributes
        attrs = field.attrs.copy()
        attrs["expression"] = field.expression
        if coords:
            attrs["coordinates"] = " ".join(coords)
        tasks.append(
            root.create_array(
                name=output_name,
                shape=shape,
                chunks=chunks,
                dtype=field.dtype,
                config=dict(write_empty_chunks=True),
                fill_value=field.fill_value,
                attributes=attrs,
                dimension_names=dims,
            )
        )
    arrays = await asyncio.gather(*tasks)
    return arrays, needs_time, needs_time_av


async def _add_time_coordinate(
    root: zarr.AsyncGroup, name: str, attrs: dict[str, str], chunk_size: int
) -> zarr.AsyncArray:
    return await root.create_array(
        name=name,
        shape=(0,),
        chunks=(chunk_size,),
        dtype=float,
        config=dict(write_empty_chunks=True),
        attributes=attrs,
        dimension_names=("time",),
    )


class ZarrGroup(File):
    def __init__(
        self,
        available_fields: Mapping[str, pygetm.core.Array],
        logger: logging.Logger,
        store: StoreLike,
        rank: int,
        chunk_size: int = 10,
        **kwargs,
    ):
        """Create a zarr group for output

        Args:
            available_fields: collection of model fields that may be added
            logger: target for log messages
            store: zarr store to create. If it exists it will be clobbered.
            rank: MPI rank of the current process
            chunk_size: number of time steps to buffer before writing to disk
                This is equivalent to the chunk size along the time dimension.
            **kwargs: additional keyword arguments passed to :class:`pygetm.output.File`
        """
        super().__init__(available_fields, logger, **kwargs)
        self.store = store
        self.root: Optional[zarr.AsyncGroup] = None
        self.is_root = rank == 0
        self._varying_fields: list[tuple[operators.Base, zarr.Array]] = []
        self.itime = 0
        self.time_offset = 0.0
        self._time_cache: Optional[np.ndarray] = None
        self._time_av_cache: Optional[np.ndarray] = None
        self._cache: list[tuple[zarr.AsyncArray, np.ndarray]] = []
        self._ncache = 0
        self._chunk_size = chunk_size
        self.loop = asyncio.new_event_loop()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.path}')"

    def start_now(
        self,
        seconds_passed: float,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ) -> bool:
        loop = self.loop
        field2array: dict[operators.Base, zarr.AsyncArray] = {}
        if self.is_root or self.sub:
            included_fields = {}
            time_array: Optional[zarr.AsyncArray] = None
            for output_name, field in self.fields.items():
                if 0 in field.shape:
                    self._logger.warning(
                        f"Skipping {output_name} because it contains no data"
                        f" (shape={field.shape})"
                    )
                else:
                    included_fields[output_name] = field

            if included_fields:
                self.root = loop.run_until_complete(
                    zarr.create_group(store=self.store, overwrite=True)
                )
                arrays, needs_time, needs_time_av = loop.run_until_complete(
                    _create_arrays(self.root, included_fields, self._chunk_size)
                )
                field2array.update(zip(included_fields.values(), arrays))
                attrs, self.time_offset = self.get_cf_time_attrs(
                    time, seconds_passed, default_time_reference or time
                )
                if needs_time:
                    time_array = loop.run_until_complete(
                        _add_time_coordinate(self.root, "time", attrs, self._chunk_size)
                    )
                # if needs_time_av:
                #     self.time_av_array = _add_time_coordinate(
                #         self.root, "time_av", attrs
                #     )

            if time_array is not None:
                self._time_cache = np.empty((self._chunk_size,), dtype=float)
                self._cache.append((time_array, self._time_cache))

        for field in self.fields.values():
            array = field2array.get(field)
            if field.time_varying:
                # Store field for update at each time step
                cache = np.empty((self._chunk_size,) + field.shape, dtype=field.dtype)
                self._varying_fields.append((field, cache))
                self._cache.append((array, cache))
            else:
                # Write static field now
                loop.run_until_complete(array.setitem(Ellipsis, field.get()))

        return len(self._varying_fields) > 0

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self._time_cache is not None:
            self._time_cache[self._ncache] = self.time_offset + seconds_passed
        for field, cache in self._varying_fields:
            field.get(cache[self._ncache, ...])
        self._ncache += 1
        self.itime += 1
        if self._ncache == self._chunk_size:
            self._empty_cache()

    def _empty_cache(self):
        istop = self.itime
        istart = self.itime - self._ncache

        async def resize_and_set(array: zarr.AsyncArray, data: np.ndarray):
            await array.resize((istop,) + data.shape[1:])
            await array.setitem((slice(istart, istop), Ellipsis), data[: self._ncache])

        tasks = []
        for array, data in self._cache:
            tasks.append(self.loop.create_task(resize_and_set(array, data)))
        self.loop.run_until_complete(asyncio.gather(*tasks))

        self._ncache = 0

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self._ncache > 0:
            self._empty_cache()
