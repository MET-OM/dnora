import numpy as np
from abc import ABC, abstractmethod
from dnora import msg
from geo_skeletons import PointSkeleton


class MaskSetter(ABC):
    """Set points (boundary, spec etc.) in the grid.

    The dimensions and orientation of the boolean array [True = boundary point]
    that is returned to the object should be:

    rows = latitude and colums = longitude (i.e.) shape = (nr_lat, nr_lon).

    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """

    @abstractmethod
    def __call__(self, grid) -> np.ndarray:
        """This method is called from within the Grid-object."""
        return mask

    @abstractmethod
    def __str__(self):
        """Describes how the boundary points are set.

        This is called by the Grid-objeect to provide output to the user.
        """
        pass


class Clear(MaskSetter):
    """Clears all boundary points by setting a mask with False values."""

    def __init__(self):
        pass

    def __call__(self, grid):
        mask_size = grid.shape("topo")
        return np.full(mask_size, False)

    def __str__(self):
        return "Clearing all possible mask points and setting an empty mask."


class All(MaskSetter):
    """Set all points to boundary points."""

    def __call__(self, grid):
        return np.full(grid.size("spatial"), True)

    def __str__(self):
        return f"Setting all points to boundary points."


class LonLat(MaskSetter):
    """Sets a list of lon, lat points to interest points"""

    def __init__(self, lon=np.ndarray, lat=np.ndarray):
        self._points = PointSkeleton(lon=lon, lat=lat)

    def __call__(self, grid):
        ind_dict = grid.yank_point(
            lon=self._points.lon(), lat=self._points.lat(), fast=True
        )

        mask = np.full(grid.shape("topo"), False)

        if grid.is_gridded():
            mask[ind_dict.get("inds_y"), ind_dict.get("inds_x")] = True
        else:
            mask[ind_dict.get("inds")] = True

        return mask

    def __str__(self):
        return f"Setting given lon, lat points to mask points"


class XY(MaskSetter):
    """Sets a list of x, y points to interest points"""

    def __init__(self, x=np.ndarray, y=np.ndarray):
        self._points = PointSkeleton(x=x, y=y)

    def __call__(self, grid):
        self._points.set_utm(grid.utm(), silent=True)
        ind_dict = grid.yank_point(
            lon=self._points.y(), lat=self._points.y(), fast=True
        )

        mask = np.full(grid.shape("topo"), False)

        if grid.is_gridded():
            mask[ind_dict.get("inds_y"), ind_dict.get("inds_x")] = True
        else:
            mask[ind_dict.get("inds")] = True

        return mask

    def __str__(self):
        return f"Setting given x, y points to mask points"


class Edges(MaskSetter):
    """Set the grid edges as mask points.

    Any combination of North, South, East, West ['N', 'S', 'E', 'W'] edges
    can be set.

    If step is e.g. 5, then only every fifth point of the edges are set. This
    is useful if the boundary spectra are coarse and we want to let the wave
    model interpolate the spectra.
    """

    def __init__(self, edges: list[str] = ["N", "S", "E", "W"], step: int = 1) -> None:
        self.edges = [edge.upper() for edge in edges]
        if step < 1:
            raise ValueError("step cannot be smaller than 1")
        else:
            self.step = int(step)
        return

    def __call__(self, grid):
        mask_size = grid.shape("topo")

        if mask_size == (1, 1):
            return np.full(mask_size, True)

        mask = np.full(mask_size, False)
        # --------- North boundary ----------
        if "N" in self.edges:
            mask[-1, :: self.step] = True
        ## --------- South boundary ----------
        if "S" in self.edges:
            mask[0, :: self.step] = True
        ## --------- East boundary ----------
        if "E" in self.edges:
            mask[:: self.step, -1] = True
        ## --------- West boundary ----------
        if "W" in self.edges:
            mask[:: self.step, 0] = True

        return mask

    def __str__(self):
        return f"Setting all edges {self.edges} to mask points using step {self.step}."


class MidPoint(MaskSetter):
    """Set the middle point of grid edges as mask points.

    Any combination of North, South, East, West ['N', 'S', 'E', 'W'] edges
    can be set.
    """

    def __init__(self, edges: list[str] = ["N", "S", "E", "W"]) -> None:
        self.edges = [edge.upper() for edge in edges]

    def __call__(self, grid):
        mask_size = grid.shape("topo")
        if mask_size == (1, 1):
            return np.full(mask_size, True)

        mask = np.full(mask_size, False)
        ny = np.round(mask_size[0] / 2).astype(int)
        nx = np.round(mask_size[1] / 2).astype(int)

        sea_mask = grid.sea_mask()

        # --------- North boundary ----------
        if "N" in self.edges:
            edge = sea_mask[-1, :]
            nx = np.round(np.median(np.where(edge))).astype(int)
            mask[-1, nx] = True
        ## --------- South boundary ----------
        if "S" in self.edges:
            edge = sea_mask[0, :]
            nx = np.round(np.median(np.where(edge))).astype(int)
            mask[0, nx] = True
        ## --------- East boundary ----------
        if "E" in self.edges:
            edge = sea_mask[:, -1]
            ny = np.round(np.median(np.where(edge))).astype(int)
            mask[ny, -1] = True
        ## --------- West boundary ----------
        if "W" in self.edges:
            edge = sea_mask[:, 0]
            ny = np.round(np.median(np.where(edge))).astype(int)
            mask[ny, 0] = True

        return mask

    def __str__(self):
        return f"Setting mid point of edges {self.edges} to mask point."
