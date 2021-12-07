import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
import numpy as np
from .grd.grd_mod import Grid, GridProcessor
from .grd.process import TrivialFilter
from .bnd.bnd_mod import Boundary
from .wnd.wnd_mod import Forcing
from .aux import add_suffix
from . import msg
from typing import Tuple

class GridPlotter(ABC):
    """Plots data from Grid-object."""

    @abstractmethod
    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, filename: str, plain: bool) -> Tuple:
        return fig, filename

class TopoPlotter(GridPlotter):
    def __init__(self):
        return

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, filename: str, plain: bool) -> Tuple:
        """Creates a plot of the topography when a Grid-object is provided.

        Options
        boundary:       If a boundary object from dnora.bnd is given, then the
                        boundary spectra locations are printed on the plot.

        forcing:        If a forcing object from dnora.wnd is given, then the
                        forcing grid is printed on the plot.

        filename:       Filename for possibly saving the plot. This can be
                        modified and returned to the object that saves the plot.

        plain:          True / False(default). Hint to plotter how much detail
                        the user wants into the plot.
        """

        # Create basic plot
        fig = plt.figure()
        levels = np.linspace(0, np.max(grid.topo()), 100, endpoint=True)
        plt.contourf(grid.lon(),grid.lat(),grid.topo(),levels)

        # Plot boundary points if they exist
        lonlat=grid.boundary_points()
        if lonlat.shape[0] > 0:
            plt.plot(lonlat[:,0], lonlat[:,1],'k*', label='Set boundary points')

        # Plot locations of boundary spectra
        if not plain and boundary is not None:
            plt.plot(boundary.lon(), boundary.lat(),'kx', label=f"Available spectra from {boundary.name()}")

        # Plot locations of wind forcing data points
        if not plain and forcing is not None:
            lonlat=forcing._point_list(mask=np.full(forcing.size()[1:], True))
            plt.plot(lonlat[:,0], lonlat[:,1],'r.', markersize=1.5, label=f"Forcing from {forcing.name()}")

        plt.legend(loc="upper right")
        cbar = plt.colorbar()
        cbar.set_label('Depth (m)', rotation=90)
        plt.title(f"{grid.name()} topography")

        return fig, add_suffix(filename, 'topo')

class MaskPlotter(GridPlotter):
    def __init__(self):
        return

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, filename: str, plain: bool) -> Tuple:
        """Creates a plot of the land-sea mask when a Grid-object is provided.

        Options
        boundary:       If a boundary object from dnora.bnd is given, then the
                        boundary spectra locations are printed on the plot.

        forcing:        If a forcing object from dnora.wnd is given, then the
                        forcing grid is printed on the plot.

        filename:       Filename for possibly saving the plot. This can be
                        modified and returned to the object that saves the plot.

        plain:          True / False(default). Hint to plotter how much detail
                        the user wants into the plot.
        """

        fig = plt.figure()
        plt.contourf(grid.lon(),grid.lat(),grid.land_sea_mask())

        # Plot boundary points if they exist
        lonlat=grid.boundary_points()
        if lonlat.shape[0] > 0:
            plt.plot(lonlat[:,0], lonlat[:,1],'k*', label='Set boundary points')

        # Plot locations of boundary spectra
        if not plain and boundary is not None:
            plt.plot(boundary.lon(), boundary.lat(),'kx', label=f"Available spectra from {boundary.name()}")

        # Plot locations of wind forcing data points
        if not plain and forcing is not None:
            lonlat=forcing._point_list(mask=np.full(forcing.size()[1:], True))
            plt.plot(lonlat[:,0], lonlat[:,1],'r.', markersize=1.5, label=f"Forcing from {forcing.name()}")


        plt.legend(loc="upper right")
        cbar = plt.colorbar()
        cbar.set_label('0=Land, 1=Sea', rotation=90)
        plt.title(f"{grid.name()} land-sea mask")

        return fig, add_suffix(filename, 'mask')
