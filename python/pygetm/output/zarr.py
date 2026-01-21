from typing import Optional, TypeVar
from collections.abc import Coroutine, Mapping
import logging
import asyncio
import concurrent.futures
import threading
import atexit

import cftime
import zarr.api.asynchronous as zarr
from zarr.storage import StoreLike
import numpy as np

from . import File
from . import operators
import pygetm.core

loop: list[Optional[asyncio.AbstractEventLoop]] = [None]
iothread: list[Optional[threading.Thread]] = [None]


def cleanup_resources() -> None:
    if loop[0] is not None:
        loop[0].call_soon_threadsafe(loop[0].stop)  # Stop loop from another thread
        iothread[0].join(timeout=0.2)  # Add a timeout to avoid hanging
        loop[0].close()


atexit.register(cleanup_resources)


async def _create_array(
    root: zarr.AsyncGroup, name: str, field: operators.Base, chunk_size: int
) -> zarr.AsyncArray:
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
    # Variable attributes
    attrs = field.attrs.copy()
    attrs["expression"] = field.expression
    if coords:
        attrs["coordinates"] = " ".join(coords)
    return await root.create_array(
        name=name,
        shape=shape,
        chunks=chunks,
        dtype=field.dtype,
        config=dict(write_empty_chunks=True),
        fill_value=field.fill_value,
        attributes=attrs,
        dimension_names=dims,
    )


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


async def _resize_and_set(array: zarr.AsyncArray, data: np.ndarray, istop: int):
    await array.resize((istop,) + data.shape[1:])
    istart = istop - data.shape[0]
    await array.setitem((slice(istart, istop),), data)


T = TypeVar("T")


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
        if loop[0] is None:
            loop[0] = asyncio.new_event_loop()
            iothread[0] = threading.Thread(target=loop[0].run_forever, daemon=True)
            iothread[0].start()
        self._futures: list[concurrent.futures.Future] = []

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.path}')"

    def start_now(
        self,
        seconds_passed: float,
        time: Optional[cftime.datetime],
        default_time_reference: Optional[cftime.datetime],
    ) -> bool:
        field2array: dict[operators.Base, zarr.AsyncArray] = {}
        if self.is_root or self.sub:
            included_fields = self.select_nonempty_fields()
            if included_fields:
                self.root = self._run(
                    zarr.create_group(store=self.store, overwrite=True)
                )

                needs_time = False
                needs_time_av = False
                for name, field in included_fields.items():
                    if field.time_varying:
                        if "time: mean" in field.attrs.get("cell_methods", ""):
                            needs_time_av = True
                        else:
                            needs_time = True
                    self._schedule(
                        _create_array(self.root, name, field, self._chunk_size)
                    )

                attrs, self.time_offset = self.get_cf_time_attrs(
                    time, seconds_passed, default_time_reference or time
                )

                if needs_time:
                    self._schedule(
                        _add_time_coordinate(self.root, "time", attrs, self._chunk_size)
                    )
                if needs_time_av:
                    self._schedule(
                        _add_time_coordinate(
                            self.root, "time_av", attrs, self._chunk_size
                        )
                    )
                    self.previous_seconds_passed = seconds_passed

                arrays = self._complete_scheduled_tasks()
                if needs_time_av:
                    self._time_av_cache = np.empty((self._chunk_size,), dtype=float)
                    self._cache.append((arrays.pop(), self._time_av_cache))
                if needs_time:
                    self._time_cache = np.empty((self._chunk_size,), dtype=float)
                    self._cache.append((arrays.pop(), self._time_cache))
                field2array.update(zip(included_fields.values(), arrays))

        for field in self.fields.values():
            array = field2array.get(field)
            if field.time_varying:
                # Store field for update at each time step
                cache = np.empty((self._chunk_size,) + field.shape, dtype=field.dtype)
                self._varying_fields.append((field, cache))
                self._cache.append((array, cache))
            else:
                # Write static field now
                self._schedule(array.setitem((), field.get()))

        return len(self._varying_fields) > 0

    def save_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self._time_cache is not None:
            self._time_cache[self._ncache] = self.time_offset + seconds_passed
        if self._time_av_cache is not None:
            self._time_av_cache[self._ncache] = self.time_offset + 0.5 * (
                self.previous_seconds_passed + seconds_passed
            )
            self.previous_seconds_passed = seconds_passed
        for field, cache in self._varying_fields:
            field.get(cache[self._ncache, ...])
        self.itime += 1
        self._ncache += 1
        if self._ncache == self._chunk_size:
            self._empty_cache()

    def close_now(self, seconds_passed: float, time: Optional[cftime.datetime]):
        if self._ncache > 0:
            self._empty_cache()
        self._complete_scheduled_tasks()

    def _empty_cache(self):
        self._complete_scheduled_tasks()
        for array, data in self._cache:
            self._schedule(
                _resize_and_set(array, data[: self._ncache].copy(), self.itime)
            )

        self._ncache = 0

    def _run(self, coro: Coroutine[None, None, T]) -> T:
        return asyncio.run_coroutine_threadsafe(coro, loop[0]).result()

    def _schedule(self, coro: Coroutine[None, None, T]):
        self._futures.append(asyncio.run_coroutine_threadsafe(coro, loop[0]))

    def _complete_scheduled_tasks(self) -> list:
        results = [future.result() for future in self._futures]
        self._futures.clear()
        return results
