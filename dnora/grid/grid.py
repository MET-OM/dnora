from __future__ import annotations
from typing import TYPE_CHECKING, Union
from geo_skeletons import GriddedSkeleton
import numpy as np
import xarray as xr
from geo_skeletons.decorators import add_mask, add_datavar

from dnora import msg

from .mesh import Mesher, Interpolate
from .process import GridProcessor
from pathlib import Path

if TYPE_CHECKING:
    from dnora.read.abstract_readers import DataReader
from dnora.type_manager.data_sources import DataSource, data_source_from_string
from dnora.type_manager.dnora_types import DnoraDataType
import matplotlib.pyplot as plt
from .topo import import_topo
import cmocean.cm

import geo_parameters as gp

from dnora.defaults import read_environment_variable


@add_mask(name="waveseries", coord_group="grid", default_value=0)
@add_mask(name="boundary", coord_group="grid", default_value=0)
@add_mask(name="output", coord_group="grid", default_value=0)
@add_mask(
    name="sea",
    coord_group="grid",
    default_value=1,
    opposite_name="land",
    triggered_by="topo",
    valid_range=(0, None),
    range_inclusive=False,
)
@add_datavar(gp.ocean.WaterDepth("topo"), default_value=999.0)
class Grid(GriddedSkeleton):
    _default_reader = None

    @classmethod
    def from_swan_bot(cls, filename: str, lon: tuple, lat: tuple, folder: str = ""):
        """Recreate the bathymetrical grid from a SWAN bot-file"""
        filepath = Path(folder).joinpath(filename)
        msg.from_file(filepath)
        topo = np.loadtxt(filepath)

        grid = cls(lon=lon, lat=lat)
        grid.set_spacing(nx=topo.shape[1], ny=topo.shape[0])
        grid.set_topo(topo)
        return grid

    def plot(self) -> None:
        vmin, vmax = np.min(self.topo()), np.max(self.topo())
        if vmax - vmin < 20:
            levels = np.linspace(vmin, vmax, np.floor(vmax - vmin + 1).astype(int))
        else:
            levels = np.linspace(vmin, vmax, 11)
        xx, yy = np.meshgrid(self.lon(native=True), self.lat(native=True))

        if len(levels) > 1:
            cont = plt.contourf(xx, yy, self.topo(), cmap=cmocean.cm.deep)

        cbar = plt.colorbar(cont, label=f"Water depth [m]")
        plt.xlabel(self.core.x_str)
        plt.ylabel(self.core.y_str)

        plt.show()

    def import_topo(
        self,
        topo_reader: DataReader = None,
        source: Union[str, DataSource] = DataSource.LOCAL,
        folder: str = None,
        **kwargs,
    ) -> None:
        source = data_source_from_string(source)
        folder = read_environment_variable(DnoraDataType.GRID, source)
        topo_reader = topo_reader or self._default_reader
        raw_topo = import_topo(self, topo_reader, source, folder, **kwargs)
        self._raw = raw_topo

    def mesh_grid(self, mesher: Mesher = Interpolate(), **kwargs) -> None:
        """Meshes the raw data down to the grid definitions."""
        if self.raw() is None:
            msg.warning("Import topography using .import_topo() before meshing!")
            return

        msg.header(mesher, "Meshing grid bathymetry...")

        xQ, yQ = np.meshgrid(self.x(native=True), self.y(native=True))

        if self.core.is_cartesian():
            x, y = self.raw().xy()
        else:
            x, y = self.raw().lonlat()

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ, **kwargs)
        print(mesher)

        self.set_topo(topo)
        grid_name = self.name
        self.meta.set(
            self.raw().meta.get()
        )  # This might overwrite the name of the grid with the name of the raw topography
        self.name = grid_name
        self.meta.set(self.raw().ds().topo.attrs, name="topo")

    def process_grid(
        self, grid_processor: GridProcessor = None, raw: bool = False, **kwargs
    ) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing meshed grid...")

        obj = self.raw() if raw else self

        topo = grid_processor(obj, **kwargs)
        print(grid_processor)

        obj.set_topo(topo)

    def set_boundary_points(self, mask_setter) -> None:
        boundary_mask = mask_setter(self)
        self.set_boundary_mask(boundary_mask)

    def set_output_points(self, mask_setter) -> None:
        mask = mask_setter(self)
        self.set_output_mask(mask)

    def set_waveseries_points(self, mask_setter) -> None:
        mask = mask_setter(self)
        self.set_waveseries_mask(mask)

    def time(self) -> tuple:
        return (None, None)

    def raw(self):
        if hasattr(self, "_raw"):
            return self._raw
        return None

    def boundary_nx(self) -> int:
        """Return approximate number of grid points in the longitude direction"""
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff = np.median(abs_diff[abs_diff > 0]).astype(int)

        return np.ceil(self.nx() / abs_diff).astype(int)

    def boundary_ny(self) -> int:
        """Return approximate number of grid points in the longitude direction"""
        abs_diff = np.abs(np.diff(np.where(self.boundary_mask())))
        if abs_diff.size == 0:
            return 0
        abs_diff = np.median(abs_diff[abs_diff > 0]).astype(int)

        return np.ceil(self.ny() / abs_diff).astype(int)

    def to_netcdf(self, filename: str = "dnora_grid.nc", folder: str = "") -> None:
        """Exports grid to netcdf file"""
        filepath = Path(folder).joinpath(filename)
        msg.to_file(filepath)
        self.ds().to_netcdf(filepath)

    def cfl(self, dx=None, f0=0.041180):
        """Calculates approximate time step [s].
        Based on grid resolution and given lowest frequency [Hz] (default=0.041180)
        """
        if dx is None:
            dx = min(self.dx(), self.dy())  # Grid spacing [m]

        cg = 1.56 / f0 * 0.5  # Deep water group velocity [m/s]
        dt = dx / cg

        print(f"Grid spacing dx = {dx:.0f} m and f0 = {f0:.8f} Hz")
        print(
            f"Approximate minimum time step: dt = dx/cg = {dx:.0f}/{cg:.1f} = {dt:.1f} s"
        )

        return dt

    def coord_dict(self, strict: bool = True):
        return {
            "lon": self.lon(strict=strict),
            "lat": self.lat(strict=strict),
            "x": self.x(strict=strict),
            "y": self.y(strict=strict),
        }
