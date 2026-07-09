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
from dnora.read.triang import MshReader, TxtReader
from .triangulation import TriangulationProcessor
from .triangulation.funcs import reindex_triangulation_when_nodes_removed
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
@add_mask(name="bad", coord_group="grid", default_value=0)
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
            tri_grid.utm.set((zone_number, zone_letter))
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
    def from_txt(cls, filename: str, boundary_filename: str = None, **kwargs):
        tri_grid = cls.generate(triang_reader=TxtReader(), filename=filename, boundary_filename=boundary_filename,**kwargs)

        return tri_grid

    def to_lonlat(self):
        """Gives an identical grid but in lon, lat coordinates"""
        cls = self.__class__
        lon, lat = self.lonlat()
        new_grid = cls(lon=lon, lat=lat, name=self.name, ntriang=range(len(self.triangles())), corner=range(3))
        new_grid.set_triangles(self.triangles())
        new_grid.set_boundary_mask(self.boundary_mask())
        new_grid.utm.set(self.utm.zone())
        new_grid.topo(self.topo(strict=True))
        
        return new_grid

    def to_xy(self):
        """Gives an identical grid but in x, y coordinates"""
        cls = self.__class__
        x, y = self.xy()
        new_grid = cls(x=x, y=y, name=self.name, ntriang=range(len(self.triangles())), corner=range(3))
        new_grid.set_triangles(self.triangles())
        new_grid.set_boundary_mask(self.boundary_mask())
        new_grid.utm.set(self.utm.zone())
        new_grid.topo(self.topo(strict=True))
        
        return new_grid

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
        grid_name = self.name
        self.meta.set(self.raw().meta.get())
        self.meta.set(self.raw().ds().topo.attrs, name="topo")
        self.name = grid_name

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

    def plot(self, edge_numbers: bool=False, show_node: int=None) -> None:
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
        plt.triplot(tri, color='black', linewidth=0.5)

        blon, blat = self.boundary_points()
        plt.scatter(blon, blat, 4,'b')
        
        blon, blat = self.bad_points()
        plt.scatter(blon, blat, 4,'r')

        if edge_numbers:
            from .triangulation.funcs import get_boundary_edges, get_boundary_nodes
            bnd_edges = get_boundary_edges(self.triangles())
            bnd_nodes = get_boundary_nodes(bnd_edges)
            bnd_nodes = bnd_nodes[0::20]
            
            lon,lat=self.lonlat(native=True)
            
            for i, b in enumerate(bnd_nodes):
                plt.text(lon[b], lat[b], str(b), fontsize=10, color="red")
        if show_node is not None:
            lon,lat=self.lonlat(native=True)
            plt.annotate(
                '', xy=(lon[show_node+1], lat[show_node+1]), xytext=(lon[show_node], lat[show_node]),
                arrowprops=dict(arrowstyle='->', lw=2, color='m', mutation_scale=15)
            )
            plt.text(lon[show_node], lat[show_node], str(show_node), fontsize=10, color="m")
            #plt.scatter(self.lon()[0], self.lat()[0], color='m')
        plt.xlabel(self.core.x_str)
        plt.ylabel(self.core.y_str)

        plt.show()

    def arange_triangulation(self, tri_aranger: TriangulationProcessor, flag_points: bool=False) -> None:
        """Use flag_points to not change the triangulation, but mark points as bad points in case they needs to be removed"""
        msg.process(tri_aranger)
        if flag_points:
            nodes = tri_aranger(
                self.inds(),
                np.where(self.boundary_mask())[0],
                self.triangles(),
                self.x(native=True),
                self.y(native=True),
                flag_points,
            )
            if nodes is None:
                msg.plain(f"{tri_aranger} doesn't plan on removing nodes, so none are flagged.")
                return
            self._update_bad(nodes)
            return
        
        nodes, bnd_nodes, tri, x, y = tri_aranger(
            self.inds(),
            np.where(self.boundary_mask())[0],
            self.triangles(),
            self.x(native=True),
            self.y(native=True),
            flag_points,
        )


        if len(nodes) < self.ny():
            msg.process("Nodes removed, updating triangulation indexing...")
            tri = reindex_triangulation_when_nodes_removed(tri, nodes)

        if self.core.is_cartesian():
            coords = {'x':x, 
                      'y':y}
        else:
            coords = {'lon':x, 'lat':y}
        self._init_structure(
            name=self.name,
            ntriang=range(tri.shape[0]),
            corner=range(3),
            **coords
        )

        self._update_boundary(bnd_nodes)
        self.set_triangles(tri)

    def set_boundary_inds(self, boundary_inds: list[int]) -> None:
        self._update_boundary(boundary_inds)

    def set_bad_inds(self, bad_inds: list[int]) -> None:
        self._update_bad(bad_inds)

    def _update_boundary(self, boundary_inds):
        mask = np.array([ind in boundary_inds for ind in self.inds()])
        self.set_boundary_mask(mask)

    def _update_bad(self, bad_inds):
        mask = np.array([ind in bad_inds for ind in self.inds()])
        self.set_bad_mask(mask)

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
