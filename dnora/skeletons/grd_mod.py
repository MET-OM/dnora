from skeletons import GriddedSkeleton
import numpy as np
import xarray as xr



class Grid(GriddedSkeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, name='AnonymousGrid'):
        self.name = name
        self.data = self._create_structure(x, y, lon, lat)
        self._reset_vars()

    def set_spacing(self, dlon: float=0., dlat: float=0.,
                    dx: float=0., dy: float=0., dm: float=0,
                    nx: int=0, ny: int=0, floating_edge: bool=False) -> None:
        """Defines longitude and latitude vectors based on desired spacing.

        Options (priority in this order)
        nx, ny [grid pnts]: Grid resolution is set to have nx points in
                            longitude and ny points in latitude direction.

        dlon, dlat [deg]:   Grid spacing is set as close to the given resolution
                            as possible (edges are fixed).

                            Set floating_edge=True to force exact dlon, dlat
                            and instead possibly move lon_max, lat_max slightly
                            to make it work.

        dm [m]:             Grid spacing is set as close as dm metres as
                            possible.

        dx, dy [m]:         Grid spacing is set as close as dx and xy metres as
                            possible.
        """
        def determine_nx_ny(dlon, dlat, dx, dy, dm, nx, ny):
            if nx and ny:
                return int(nx), int(ny)

            if dlon and dlat:
                nx = np.round((self.lon_edges()[1]-self.lon_edges()[0])/dlon) + 1
                ny = np.round((self.lat_edges()[1]-self.lat_edges()[0])/dlat) + 1
                return int(nx), int(ny)

            if dm:
                dx = dm
                dy = dm

            if dx and dy:
                nx = np.round((self.x_edges()[1]-self.x_edges()[0])/dx) + 1
                ny = np.round((self.y_edges()[1]-self.y_edges()[0])/dy) + 1
                return int(nx), int(ny)

            raise ValueError('Give a combination of nx/xy, dlon/dlat, dx/dy or dm')

        nx, ny = determine_nx_ny(dlon, dlat, dx, dy, dm, nx, ny)

        x = np.linspace(self.native_x_edges()[0], self.native_x_edges()[1], nx)
        y = np.linspace(self.native_y_edges()[0], self.native_y_edges()[1], ny)
        self.data = self._init_ds(x=x, y=y)
        self._reset_vars()

    def _reset_vars(self):
        self.merge_in_ds(self._list_of_empty_ds())

    def _list_of_empty_ds(self) -> list[xr.Dataset]:
        ds_topo = self.compile_to_ds(999*np.ones(self.size()),'topo')
        ds_sea = self.compile_to_ds(np.ones(self.size()),'sea_mask')
        ds_bnd = self.compile_to_ds(np.zeros(self.size()),'boundary_mask')
        return [ds_topo, ds_sea, ds_bnd]

    def sea_mask(self, logical=True) -> np.ndarray:
        """Returns bool array of the sea mask.
        Set logical=False to get 0 for land and 1 for sea. """

        if hasattr(self.data, 'sea_mask'):
            mask = self.data.sea_mask.values
        else:
            mask = np.full(self.size(), 1.)

        if logical:
            mask = mask.astype(bool)
        return mask

    def boundary_mask(self, logical=True) -> np.ndarray:
        """Returns bool array of the boundary mask.
        Set logical=False to get 1 for boundary points """

        if hasattr(self.data, 'boundary_mask'):
            mask = self.data.boundary_mask.values
        else:
            mask = np.full(self.size(), 1.)

        if logical:
            mask = mask.astype(bool)
        return mask


    def topo(self, land: float=-999) -> np.ndarray:
        if hasattr(self.data, 'topo'):
            topo = self.data.topo.values
            topo[np.logical_not(self.sea_mask())] = land
            return topo
        else:
            return np.array([])
