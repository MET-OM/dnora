from __future__ import annotations
import numpy as np
from abc import ABC, abstractmethod
from typing import Union
from geo_skeletons import PointSkeleton

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dnora.grid import Grid, TriGrid

from dnora import msg, utils


class PointPicker(ABC):
    """PointPickers take in the grid, a PointSkeleton of all available points and a PointSkeleton of interest points and returns indeces
    of the chosen points."""

    @abstractmethod
    def __call__(
        self,
        grid: Union[Grid, TriGrid],
        all_points: PointSkeleton,
        selected_points: PointSkeleton,
        **kwargs,
    ) -> np.ndarray:
        pass


class Trivial(PointPicker):
    """Choose all the points in the list."""

    def __call__(
        self, grid: Union[Grid, TriGrid], all_points: PointSkeleton, **kwargs
    ) -> np.ndarray:
        return all_points.inds()


class NearestGridPoint(PointPicker):
    """Choose the nearest grid point to each boundary point in the grid.
    Set a maximum allowed distance using `max_dist` (in km) at instantiation time.
    """

    def __call__(
        self,
        grid: Union[Grid, TriGrid],
        all_points: PointSkeleton,
        selected_points: PointSkeleton,
        max_dist: float = None,
        fast: bool = True,
        **kwargs,
    ) -> np.ndarray:
        if selected_points is None:
            msg.plain("No interest points provided. Returning empty list of indeces.")
            return np.array([])

        lon, lat = selected_points.lonlat()

        # Go through all points where we want output and find the nearest available point
        chosen_inds = []
        dx = []
        for lo, la in zip(lon, lat):
            ind_dict = all_points.yank_point(lon=lo, lat=la, fast=fast, unique=True)
            chosen_inds.append(ind_dict.get("inds")[0])
            dx.append(ind_dict.get("dx")[0])
        inds = []
        _, index = np.unique(chosen_inds, return_index=True)

        chosen_inds = [chosen_inds[i] for i in index]
        dx = [dx[i] for i in index]
        lon = [lon[i] for i in index]
        lat = [lat[i] for i in index]

        for n, (x, y, ind) in enumerate(zip(lon, lat, chosen_inds)):

            ms = f"Point {n}: lat: {y:10.7f}, lon: {x:10.7f} <<< ({all_points.lat()[ind]: .7f}, {all_points.lon()[ind]: .7f}). Distance: {dx[n]/1000:.1f} km"

            if max_dist is None or dx[n] / 1000 <= max_dist:
                msg.plain(ms)
                inds.append(ind)
            else:
                msg.plain("DISCARDED, too far: " + ms)

        inds = np.array(inds)
        return inds


class Area(PointPicker):
    """Choose all the points within a certain area around the grid."""

    def __call__(
        self,
        grid: Union[Grid, TriGrid],
        all_points: PointSkeleton,
        expansion_factor: float = 1.5,
        **kwargs,
    ) -> np.ndarray:

        msg.process(f"Using expansion_factor = {expansion_factor:.2f}")

        # Define area to search in
        if grid.core.is_cartesian():
            number, zone = grid.utm()
            all_points.set_utm(number, zone)
            x, y = utils.grid.expand_area(
                grid.edges("x"), grid.edges("y"), expansion_factor
            )
            x_all, y_all = all_points.xy()
        else:
            x, y = utils.grid.expand_area(
                grid.edges("lon"), grid.edges("lat"), expansion_factor
            )
            x_all, y_all = all_points.lonlat()

        maskx = np.logical_and(x_all >= x[0], x_all <= x[1])
        masky = np.logical_and(y_all >= y[0], y_all <= y[1])
        mask = np.logical_and(maskx, masky)

        inds = np.where(mask)[0]

        msg.info(
            f"Found {len(inds)} points inside {x[0]:10.7f}-{x[1]:10.7f}, {y[0]:10.7f}-{y[1]:10.7f}."
        )

        return inds
