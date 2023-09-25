from geo_skeletons import GriddedSkeleton, PointSkeleton
from geo_skeletons.decorators import add_datavar
from .. import aux_funcs
from .. import msg
import xarray as xr
from .read import TopoReader
from .mesh import Mesher, Interpolate
from .process import GridProcessor
from copy import copy
import numpy as np


@add_datavar(name="topo", default_value=999.0)
class GriddedTopo(GriddedSkeleton):
    pass


@add_datavar(name="topo", default_value=999.0)
class PointTopo(PointSkeleton):
    pass


class GridMethods:
    @classmethod
    def from_ds(cls, ds: xr.Dataset, topo_var_name: str = "topo"):
        lon, lat, x, y = aux_funcs.get_coordinates_from_ds(ds)
        grid = cls(lon=lon, lat=lat, x=x, y=y)

        if grid.is_gridded():
            grid.set_spacing(nx=len(x or lon), ny=len(y or lat))
        topo = ds.get(topo_var_name)
        sea_mask = ds.get("sea_mask")
        boundary_mask = ds.get("boundary_mask")

        if topo is not None:
            grid.set_topo(topo.values)

        if sea_mask is not None:
            grid.set_sea_mask(sea_mask.values)

        if boundary_mask is not None:
            grid.set_boundary_mask(boundary_mask.values)

        return grid

    @classmethod
    def from_netcdf(cls, filename: str):
        ds = xr.open_dataset(filename)
        grid = cls.from_ds(ds, topo_var_name="topo")
        return grid

    def export_grid(self, filename: str = "dnora_grid") -> None:
        """Exports grid to netcdf file"""
        self.ds().to_netcdf(filename + ".nc")

    def import_topo(self, topo_reader: TopoReader) -> None:
        """Reads the raw bathymetrical data."""

        msg.header(topo_reader, "Importing topography...")
        print(topo_reader)
        topo, lon, lat, x, y, zone_number, zone_letter = topo_reader(
            self.edges("lon"), self.edges("lat"), self.edges("x"), self.edges("y")
        )

        if 0 in topo.shape:
            msg.warning(
                "Imported topography seems to be empty. Maybe using wrong tile?"
            )
            return

        if aux_funcs.is_gridded(topo, lon, lat) or aux_funcs.is_gridded(topo, x, y):
            self._raw = GriddedTopo(lon=lon, lat=lat, x=x, y=y)
            x = x or lon
            y = y or lat
            self._raw.set_spacing(nx=len(x), ny=len(y))
        else:
            self._raw = PointTopo(lon=lon, lat=lat, x=x, y=y)

        if (
            self.edges("lon", native=True)[0] < self.raw().edges(self.x_str)[0]
            or self.edges("lon", native=True)[1] > self.raw().edges(self.x_str)[1]
        ):
            msg.warning(
                f"The data gotten from the TopoReader doesn't cover the grid in the {self.x_str} direction. Grid: {self.edges('lon')}, imported topo: {self.raw().edges(self.x_str)}"
            )

        if (
            self.edges("lat", native=True)[0] < self.raw().edges(self.y_str)[0]
            or self.edges("lat", native=True)[1] > self.raw().edges(self.y_str)[1]
        ):
            msg.warning(
                f"The data gotten from the TopoReader doesn't cover the grid in the {self.y_str} direction. Grid: {self.edges('lat')}, imported topo: {self.raw().edges(self.y_str)}"
            )

        if zone_number is not None:
            self._raw.set_utm(zone_number, zone_letter)

        self.raw().set_topo(topo)

    def mesh_grid(self, mesher: Mesher = Interpolate(method="nearest")) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.raw() is None:
            msg.warning("Import topography using .import_topo() before meshing!")
            return

        msg.header(mesher, "Meshing grid bathymetry...")
        print(mesher)

        if self.is_gridded():
            xQ, yQ = np.meshgrid(self.x(native=True), self.y(native=True))
        else:
            xQ, yQ = self.xy(native=True)

        if self.is_cartesian():
            x, y = self.raw().xy()
        else:
            x, y = self.raw().lonlat()

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ)

        self.set_topo(topo)
        self.set_sea_mask(self.topo() > 0)

    def process_grid(self, grid_processor: GridProcessor = None) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing meshed grid...")
        print(grid_processor)
        if self.is_gridded():
            topo = grid_processor.grid(
                self.topo(),
                self.lon(),
                self.lat(),
                self.sea_mask(),
                self.boundary_mask(),
            )
            if topo is None:
                msg.warning(
                    "Filtering of gridded topography is not implemented in this GridProcessor."
                )
                return
        else:
            topo = grid_processor.topo(
                self.topo(), self.lon(), self.lat(), self.sea_mask()
            )
            if topo is None:
                msg.warning(
                    "Filtering of unstructured topography is not implemented in this GridProcessor."
                )
                return

        self.set_topo(topo)
        self.set_sea_mask(self.topo() > 0)

    def set_boundary_points(self, mask_setter) -> None:
        boundary_mask = mask_setter(self)
        self.set_boundary_mask(boundary_mask)

    def set_output_points(self, mask_setter) -> None:
        mask = mask_setter(self)
        self.set_output_mask(mask)

    def time(self) -> tuple:
        return (None, None)

    def tri(self):
        if hasattr(self, "_tri"):
            return copy(self._tri)
        else:
            return None

    def raw(self):
        if hasattr(self, "_raw"):
            return self._raw
        return None

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
