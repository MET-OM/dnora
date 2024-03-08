from __future__ import annotations
import pandas as pd
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .tiles import TileObject


def single_patch(
    tiles: TileObject,
) -> tuple[list[tuple], list[tuple], list[tuple], str]:
    patch_dates = _find_patches_in_time(tiles)
    lon, lat = tiles.lonlat(tiles.covering_files())
    lon = np.array(lon)
    lat = np.array(lat)
    if not patch_dates:
        return patch_dates

    patch_dates = pd.to_datetime(patch_dates)

    patch_dimension = "time"
    return (
        [(min(patch_dates), max(patch_dates) + pd.Timedelta(hours=23.99))],
        [(np.min(lon), np.max(lon))],
        [(np.min(lat), np.max(lat))],
        patch_dimension,
    )


def patch_in_time(
    tiles: TileObject,
) -> tuple[list[tuple], list[tuple], list[tuple], str]:
    patch_dates = _find_patches_in_time(tiles)
    lon, lat = tiles.lonlat(tiles.covering_files())
    lon = np.array(lon)
    lat = np.array(lat)
    if not patch_dates:
        return patch_dates

    patch_dates = pd.to_datetime(patch_dates)
    patch_time = [(start, start + pd.Timedelta(hours=23.99)) for start in patch_dates]

    patch_lon = [(np.min(lon), np.max(lon)) for _ in range(len(patch_time))]
    patch_lat = [(np.min(lat), np.max(lat)) for _ in range(len(patch_time))]
    patch_dimension = "time"
    return patch_time, patch_lon, patch_lat, patch_dimension


def _find_patches_in_time(tiles: TileObject) -> list[str]:
    all_lon, all_lat = tiles.lonlat(tiles.covering_files())

    number_of_spatial_tiles = len(
        np.unique([f"{la[0]:03.0f}N{lo[0]:03.0f}E" for lo, la in zip(all_lon, all_lat)])
    )
    all_year, all_month, all_day = tiles.times(tiles.covering_files())
    year, month, day = tiles.times(tiles.relevant_files())
    patch_dates = []

    for y, m, d in zip(all_year, all_month, all_day):
        tiles_for_day = sum(
            np.logical_and(
                np.logical_and(np.array(year) == y, np.array(month) == m),
                np.array(day) == d,
            )
        )
        # If even one tile is missing for that day, we will read the day for all area
        if tiles_for_day < number_of_spatial_tiles:
            patch_dates.append(f"{y:04.0f}-{m:02.0f}-{d:02.0f}")

    return patch_dates
