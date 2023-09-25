from geo_skeletons import PointSkeleton
import numpy as np
from geo_skeletons.decorators import add_mask, add_datavar
from .grid_methods import GridMethods

from .read import MshFile as topo_MshFile
from .read_tr import TriangReader
from .read_tr import MshFile as triang_MshFile
from .tri_arangers import TriAranger
from .mesh import Trivial as TrivialMesher


@add_datavar(name="topo", default_value=999.0)
@add_mask(name="boundary", coords="grid", default_value=0)
@add_mask(name="sea", coords="grid", default_value=1, opposite_name="land")
class TriGrid(PointSkeleton, GridMethods):
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


# if self.x_str == 'x':
#             self._init_structure(x=x, y=y, lon=None, lat=None)
#         else:
#             self._init_structure(x=None, y=None, lon=x, lat=y)
