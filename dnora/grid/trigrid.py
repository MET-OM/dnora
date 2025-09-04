from geo_skeletons import PointSkeleton
import numpy as np
import xarray as xr
from geo_skeletons.decorators import add_mask, add_datavar, add_coord

from dnora import msg
from typing import Union
from .mesh import Mesher, Interpolate
from .process import GridProcessor
from pathlib import Path
from dnora.read.grid.grid_readers import MshFile as topo_MshFile
from dnora.read.triang import MshReader
from .tri_arangers import TriAranger
from .mesh import Trivial as TrivialMesher
from dnora.read.abstract_readers import DataReader
import cmocean.cm
from dnora.type_manager.data_sources import DataSource

from pathlib import Path
from .topo import import_topo

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import geo_parameters as gp


@add_datavar(name="triangles", coord_group="gridpoint")
@add_coord(name="corner", grid_coord=False)
@add_coord(name="ntriang", grid_coord=False)
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
@add_datavar(gp.ocean.WaterDepth("topo"), default_value=999.0, coord_group="grid")
class TriGrid(PointSkeleton):
    _default_reader = None

    def __init__(
        self,
        x=None,
        y=None,
        lon=None,
        lat=None,
        ntriang=range(1),
        corner=range(3),
        **kwargs,
    ):
        super().__init__(
            x=x, y=y, lon=lon, lat=lat, ntriang=ntriang, corner=corner, **kwargs
        )

    @classmethod
    def generate(
        cls,
        triang_reader: DataReader,
        folder: str = None,
        name: str = "LonelyGrid",
        **kwargs,
    ):
        (
            tri,
            coord_dict,
            edge_nodes,
            zone_number,
            zone_letter,
        ) = triang_reader(source=DataSource.LOCAL, folder=folder, **kwargs)

        tri_grid = cls(
            x=coord_dict.get("x"),
            y=coord_dict.get("y"),
            lon=coord_dict.get("lon"),
            lat=coord_dict.get("lat"),
            name=name,
            ntriang=range(tri.shape[0]),
            corner=range(3),
        )
        if zone_number is not None:
            tri_grid.set_utm((zone_number, zone_letter))
        edge_nodes = np.array(edge_nodes)
        edge_nodes = edge_nodes.astype(int)
        tri_grid._update_boundary(edge_nodes)
        tri_grid.set_triangles(tri)
        return tri_grid

    @classmethod
    def from_msh(cls, filename: str, read_topo: bool = True, **kwargs):
        tri_grid = cls.generate(triang_reader=MshReader(), filename=filename, **kwargs)

        if read_topo:
            tri_grid.import_topo(topo_MshFile(), filename=filename)
            tri_grid.mesh_grid(TrivialMesher())

        return tri_grid

    @classmethod
    def from_netcdf(cls, filename: str, folder: str = ""):
        filepath = Path(folder).joinpath(filename)
        msg.from_file(filepath)
        ds = xr.open_dataset(filepath)
        grid = cls.from_ds(ds)
        return grid

    def import_topo(
        self,
        topo_reader: DataReader = None,
        source: Union[str, DataSource] = None,
        folder: str = None,
        **kwargs,
    ) -> None:
        topo_reader = topo_reader or self._default_reader
        raw_topo = import_topo(self, topo_reader, source, folder, **kwargs)
        self._raw = raw_topo

    def mesh_grid(self, mesher: Mesher = Interpolate(), **kwargs) -> None:
        """Meshes the raw data down to the grid definitions."""
        if self.raw() is None:
            msg.warning("Import topography using .import_topo() before meshing!")
            return

        msg.header(mesher, "Meshing grid bathymetry...")

        xQ, yQ = self.xy(native=True)

        x, y = self.raw().xy(native=True)

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ, **kwargs)
        print(mesher)

        self.set_topo(topo)

        self.meta.set(self.raw().meta.get())
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

    def plot(self) -> None:
        vmin, vmax = np.min(self.topo()), np.max(self.topo())
        if vmax - vmin < 20:
            levels = np.linspace(vmin, vmax, np.floor(vmax - vmin + 1).astype(int))
        else:
            levels = np.linspace(vmin, vmax, 11)

        tri = mtri.Triangulation(
            self.lon(native=True), self.lat(native=True), triangles=self.triangles()
        )

        if len(levels) > 1:
            cont = plt.tricontourf(
                tri, self.topo(), cmap=cmocean.cm.deep, levels=levels
            )
        cbar = plt.colorbar(cont, label=f"Water depth [m]")
        plt.xlabel(self.core.x_str)
        plt.ylabel(self.core.y_str)

        plt.show()

    def arange_triangulation(self, tri_aranger: TriAranger) -> None:
        print(tri_aranger)
        bnd_nodes, tri, nodes, x, y = tri_aranger(
            self.inds(),
            np.where(self.boundary_mask())[0],
            self.tri(),
            self.x(native=True),
            self.y(native=True),
        )

        x, y = self.xy(strict=True)
        lon, lat = self.lonlat(strict=True)
        self._init_structure(
            x=x,
            y=y,
            lon=lon,
            lat=lat,
            name=self.name,
            ntriang=range(tri.shape[0]),
            corner=range(3),
        )

        self._update_boundary(bnd_nodes)
        self.set_triangles(tri)

    def _update_boundary(self, boundary_inds):
        mask = np.array([ind in boundary_inds for ind in self.inds()])
        self.set_boundary_mask(mask)

    def to_netcdf(self, filename: str = "dnora_grid.nc", folder: str = "") -> None:
        """Exports grid to netcdf file"""
        filepath = Path(folder).joinpath(filename)
        msg.to_file(filepath)
        self.ds().to_netcdf(filepath)

    def time(self) -> tuple:
        return (None, None)

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

    def coord_dict(self, strict: bool = True):
        return {
            "lon": self.lon(strict=strict),
            "lat": self.lat(strict=strict),
            "x": self.x(strict=strict),
            "y": self.y(strict=strict),
        }
