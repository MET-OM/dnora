from geo_skeletons import GriddedSkeleton, PointSkeleton
import numpy as np
import xarray as xr
from geo_skeletons.decorators import add_mask, add_datavar

from dnora import aux_funcs, msg
from copy import copy
from .mesh import Mesher, Interpolate
from .process import GridProcessor
from pathlib import Path
from .read_tr import TriangReader
from .read import MshFile as topo_MshFile
from .read_tr import MshFile as triang_MshFile
from .tri_arangers import TriAranger
from .mesh import Trivial as TrivialMesher
from dnora.readers.abstract_readers import DataReader
import os
from dnora.dnora_types import DnoraDataType, DataSource, data_source_from_string

from dnora.defaults import read_environment_variable

from pathlib import Path


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
    def from_netcdf(cls, filename: str, folder: str = ""):
        filepath = Path(folder).joinpath(filename)
        msg.from_file(filepath)
        ds = xr.open_dataset(filepath)
        grid = cls.from_ds(ds, topo_var_name="topo")
        return grid

    def to_netcdf(self, filename: str = "dnora_grid.nc", folder: str = "") -> None:
        """Exports grid to netcdf file"""
        filepath = Path(folder).joinpath(filename)
        msg.to_file(filepath)
        self.ds().to_netcdf(filepath)

    def import_topo(
        self,
        topo_reader: DataReader = None,
        source: str | DataSource = None,
        folder: str = None,
        **kwargs,
    ) -> None:
        """Reads the raw bathymetrical data."""
        topo_reader = self._default_reader or topo_reader
        if topo_reader is None:
            raise ValueError("Define a DataReader!")
        msg.header(topo_reader, "Importing topography...")

        source = source or topo_reader.default_data_source()
        source = data_source_from_string(source)

        folder = folder or read_environment_variable(DnoraDataType.GRID, source)
        if folder and source == DataSource.LOCAL:
            if not os.path.exists(os.path.expanduser(folder)):
                os.mkdir(folder)
            if not os.path.exists(aux_funcs.get_url(folder, topo_reader.name())):
                os.mkdir(aux_funcs.get_url(folder, topo_reader.name()))

        topo, coord_dict, zone_number, zone_letter, metadata = topo_reader(
            self, source=source, folder=folder, **kwargs
        )
        print(topo_reader)

        lon, lat, x, y = (
            coord_dict.get("lon"),
            coord_dict.get("lat"),
            coord_dict.get("x"),
            coord_dict.get("y"),
        )

        if 0 in topo.shape:
            msg.warning(
                "Imported topography seems to be empty. Maybe using wrong tile?"
            )
            return

        if aux_funcs.is_gridded(topo, lon, lat) or aux_funcs.is_gridded(topo, x, y):
            self._raw = Grid(lon=lon, lat=lat, x=x, y=y)
            x = x or lon
            y = y or lat
            self._raw.set_spacing(nx=len(x), ny=len(y))
        else:
            self._raw = UnstrGrid(lon=lon, lat=lat, x=x, y=y)

        if (
            self.edges("lon", native=True)[0] < self.raw().edges(self.x_str)[0]
            or self.edges("lon", native=True)[1] > self.raw().edges(self.x_str)[1]
        ):
            msg.warning(
                f"The data gotten from the DataReader doesn't cover the grid in the {self.x_str} direction. Grid: {self.edges('lon', native=True)}, imported topo: {self.raw().edges(self.x_str)}"
            )

        if (
            self.edges("lat", native=True)[0] < self.raw().edges(self.y_str)[0]
            or self.edges("lat", native=True)[1] > self.raw().edges(self.y_str)[1]
        ):
            msg.warning(
                f"The data gotten from the DataReader doesn't cover the grid in the {self.y_str} direction. Grid: {self.edges('lat', native=True)}, imported topo: {self.raw().edges(self.y_str)}"
            )

        if zone_number is not None:
            self._raw.set_utm((zone_number, zone_letter))

        self.raw().set_topo(topo)
        self.raw().set_metadata(metadata)

    def mesh_grid(self, mesher: Mesher = Interpolate(), **kwargs) -> None:
        """Meshes the raw data down to the grid definitions."""

        if self.raw() is None:
            msg.warning("Import topography using .import_topo() before meshing!")
            return

        msg.header(mesher, "Meshing grid bathymetry...")

        if self.is_gridded():
            xQ, yQ = np.meshgrid(self.x(native=True), self.y(native=True))
        else:
            xQ, yQ = self.xy(native=True)

        if self.is_cartesian():
            x, y = self.raw().xy()
        else:
            x, y = self.raw().lonlat()

        topo = mesher(self.raw().topo().ravel(), x, y, xQ, yQ, **kwargs)
        print(mesher)

        self.set_topo(topo)
        self.set_sea_mask(self.topo() > 0)
        self.set_metadata(self.raw().metadata())

    def process_grid(self, grid_processor: GridProcessor = None, **kwargs) -> None:
        """Processes the gridded bathymetrical data, e.g. with a filter."""
        if grid_processor is None:
            return

        msg.header(grid_processor, "Processing meshed grid...")

        topo = grid_processor(self, **kwargs)
        print(grid_processor)

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


@add_datavar(name="topo", default_value=999.0)
@add_mask(name="boundary", coords="grid", default_value=0)
@add_mask(name="output", coords="grid", default_value=0)
@add_mask(name="sea", coords="grid", default_value=1, opposite_name="land")
class Grid(GriddedSkeleton, GridMethods):
    _default_reader = None

    @classmethod
    def from_ww3_grid(cls, gridname: str, folder: str = ""):
        """Recreate a WW3 grid object based on the _info, _bathy and _mapsta files"""

        filename = Path(folder) / f"{gridname}_info.txt"

        print(filename)
        (
            lon_min,
            lon_max,
            lat_min,
            lat_max,
            dlon,
            dlat,
            NX,
            NY,
        ) = aux_funcs.read_ww3_info(filename)

        filename = Path(folder).joinpath(f"{gridname}_bathy.txt")
        topo = np.loadtxt(filename).reshape((NY, NX))
        filename = Path(folder).joinpath(f"{gridname}_mapsta.txt")
        mask = (
            np.loadtxt(filename).reshape((NY, NX)) == 2
        )  # Boundary points given as value 2

        grid = cls(lon=(lon_min, lon_max), lat=(lat_min, lat_max), name=gridname)
        grid.set_spacing(nx=NX, ny=NY)
        grid.set_topo(topo)
        grid.set_boundary_mask(mask)

        return grid

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


@add_datavar(name="topo", default_value=999.0)
@add_mask(name="boundary", coords="grid", default_value=0)
@add_mask(name="output", coords="grid", default_value=0)
@add_mask(name="sea", coords="grid", default_value=1, opposite_name="land")
class UnstrGrid(PointSkeleton, GridMethods):
    pass


class TriGrid(UnstrGrid):
    @classmethod
    def from_msh(cls, filename: str, name: str = "LonelyGrid"):
        tri_grid = cls(name=name)
        tri_grid.import_triang(triang_MshFile(filename))
        tri_grid.import_topo(topo_MshFile(filename))
        tri_grid.mesh_grid(TrivialMesher())

        return tri_grid

    def __init__(self, x=None, y=None, lon=None, lat=None, name="LonelyGrid"):
        self.name = name
        # Only initialize if x, y, lon, lat given
        if [a for a in (x, y, lon, lat) if a is not None]:
            self._init_structure(x, y, lon, lat)

    def import_triang(self, triang_reader: TriangReader):
        """Reads a triangular mesh."""
        (
            tri,
            nodes,
            lon,
            lat,
            x,
            y,
            types,
            edge_nodes,
            zone_number,
            zone_letter,
        ) = triang_reader()

        self._init_structure(x, y, lon, lat)

        self.set_utm(zone_number, zone_letter)
        edge_nodes = np.array(edge_nodes)
        edge_nodes = edge_nodes.astype(int)
        self._update_boundary(edge_nodes)
        self._tri = tri
        # self._nodes = nodes # These are now in self.inds()
        self._types = types  # ???

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
        self._init_structure(x=x, y=y, lon=lon, lat=lat)

        self._update_boundary(bnd_nodes)
        self._tri = tri

    def _update_boundary(self, boundary_inds):
        mask = np.array([ind in boundary_inds for ind in self.inds()])
        self.set_boundary_mask(mask)
