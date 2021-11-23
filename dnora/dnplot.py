import matplotlib.pyplot as plt
import numpy as np
from .grd.grd_mod import Grid, GridProcessor
from .grd.process import TrivialFilter
from .bnd.bnd_mod import Boundary
from .wnd.wnd_mod import Forcing

from . import msg


def grd_topo(grid: Grid, boundary_points: bool=False, save_fig: bool=False, filename: str='', boundary: Boundary=None, forcing: Forcing=None, grid_processor: GridProcessor=TrivialFilter()) -> None:
    """Creates a plot of the topography when a Grid-object is provided.

    Options
    boundary:       If a boundary object from dnora.bnd is given, then the
                    boundary spectra locations are printed on the plot.

    forcing:        If a forcing object from dnora.wnd is given, then the
                    forcing grid is printed on the plot.

    boundary_points: Show grid boundary points (default: False)

    save_fig:       If set to true, plot is saved.

    filename:       Used if save_fig=True, e.g. 'My_plot.pdf'.
                    Default is ''{grid.name()}_topo.pdf'

    grid_processor: A GridProcessor might be provided to modify the data
                    before plotting. These modifications are not applied
                    to the data, but are only used for plotting.
    """

    # Allows a modification of the topography before plotting
    topo = grid_processor(grid.topo(), grid.lon(), grid.lat(), grid.land_sea_mask(), grid.boundary_mask())

    if grid.topo().size > 0:
        if not filename:
            filename = f"{grid.name()}_topo.pdf"
        levels = np.linspace(0, np.max(topo), 100, endpoint=True)

        plt.figure()
        plt.contourf(grid.lon(),grid.lat(),topo,levels)

        if boundary_points:
            lonlat=grid.boundary_points()
            plt.plot(lonlat[:,0], lonlat[:,1],'k*', label='Set boundary points')

        if boundary is not None:
            plt.plot(boundary.lon(), boundary.lat(),'kx', label=f"Available spectra from {boundary.name()}")

        if forcing is not None:
            lonlat=forcing._point_list(mask=np.full(forcing.size()[1:], True))
            plt.plot(lonlat[:,0], lonlat[:,1],'r.', markersize=1.5, label=f"Forcing from {forcing.name()}")

        plt.legend(loc="upper right")
        cbar = plt.colorbar()
        cbar.set_label('Depth (m)', rotation=90)
        plt.title(f"{grid.name()} topography")
        plt.show()

        if save_fig:
            plt.savefig(filename, dpi=300)
            msg.to_file(filename)
    else:
        msg.templates('no_topo')

    return

def grd_mask(grid, boundary_points: bool=True, save_fig: bool=False, filename: str='', boundary: Boundary=None, forcing: Forcing=None) -> None:
    """Creates a plot of the land-sea mask based on the provided Grid object.

    Options
    boundary:       If a boundary object from dnora.bnd is given, then the
                    boundary spectra locations are printed on the plot.

    forcing:        If a forcing object from dnora.wnd is given, then the
                    forcing grid is printed on the plot.

    boundary_points: Show grid boundary points (default: True)

    save_fig:       If set to true, plot is saved.

    filename:       Used if save_fig=True, e.g. 'My_plot.pdf'.
                    Default is ''{grid.name()}_mask.pdf'
    """
    if grid.boundary_mask().size > 0:
        if not filename:
            filename = f"{grid.name()}_mask.pdf"

        plt.figure()
        plt.contourf(grid.lon(),grid.lat(),grid.land_sea_mask())

        if boundary_points:
            lonlat=grid.boundary_points()
            plt.plot(lonlat[:,0], lonlat[:,1],'k*', label='Set boundary points')

        if boundary is not None:
            plt.plot(boundary.lon(), boundary.lat(),'kx', label=f"Available spectra from {boundary.name()}")

        if forcing is not None:
            lonlat=forcing._point_list(mask=np.full(forcing.size()[1:], True))
            plt.plot(lonlat[:,0], lonlat[:,1],'r.', markersize=1.5, label=f"Forcing from {forcing.name()}")


        plt.colorbar()
        plt.title(f"{grid.name()} land-sea mask")
        plt.show()

        if save_fig:
            plt.savefig(filename, dpi=300)
            msg.to_file(filename)
    else:
        msg.templates('no_mask')

    return
