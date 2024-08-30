from __future__ import annotations
from typing import TYPE_CHECKING

from abc import ABC, abstractmethod

if TYPE_CHECKING:
    from dnora.type_manager.dnora_objects import DnoraObject
import xarray as xr


class GriddedDataProcessor(ABC):
    @abstractmethod
    def __call__(self, obj: DnoraObject) -> xr.Dataset:
        """Takes in a DnoraObject (GriddedSkeleton) and returns an xarray Dataset containing the new data"""
        pass

    @abstractmethod
    def __str__(self) -> str:
        """Describes how the spectral values as processed"""
        pass


class FillNaNs(GriddedDataProcessor):
    def __init__(self, pad: float) -> None:
        self.pad = pad

    def __call__(self, obj) -> xr.Dataset:
        """Replaces NaN-values with the given pad values"""
        return obj.ds().fillna(self.pad)

    def __str__(self) -> str:
        return f"Replacing all NaN-values with {self.pad}..."
