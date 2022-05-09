import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
import numpy as np
from ..grd.grd_mod import Grid, GridProcessor
from ..grd.process import TrivialFilter
from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing
from .. import file_module
from .. import msg
from typing import Tuple
import matplotlib.tri as mtri
from . import basic_funcs
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
class GridPlotter(ABC):
    """Plots data from Grid-object."""

    @abstractmethod
    def _extension(self) -> str:
        pass

    def grid(self, dict_of_objects: dict, plain: bool) -> Tuple:
        msg.warning('Plotting of meshed data not implemented!')
        return None, None

    def topo(self, dict_of_objects: dict, plain: bool) -> Tuple:
        msg.warning('Plotting of xyz data not implemented!')
        return None, None

class TopoPlotter(GridPlotter):
    def _extension(self):
        return 'pdf'

    def grid(self, dict_of_objects: dict, plain: bool=True):
        grid = dict_of_objects['Grid']
        boundary = dict_of_objects['Boundary']
        forcing = dict_of_objects['Forcing']
        #
        fig, ax = basic_funcs.plot_field(grid.lon(), grid.lat(), grid.topo(land=np.nan), var='topo', title_str=f"{grid.name()}")

        lonlat=grid.boundary_points()
        if lonlat.shape[0] > 0:
            ax.plot(lonlat[:,0], lonlat[:,1],'k*', label='Set boundary points')

        # Plot locations of boundary spectra
        if not plain and boundary is not None:
            ax.plot(boundary.lon(), boundary.lat(),'kx', label=f"Available spectra from {boundary.name()}")

        # Plot locations of wind forcing data points
        if not plain and forcing is not None:
            lonlat=forcing._point_list(mask=np.full(forcing.size()[1:], True))
            ax.plot(lonlat[:,0], lonlat[:,1],'r.', markersize=1.5, label=f"Forcing from {forcing.name()}")

        #ax.legend(loc="upper right")
        ax.legend(bbox_to_anchor=(0.3, 1.2))

        return fig


class ForcingPlotter(GridPlotter):
    def _extension(self):
        return 'pdf'

    def grid(self, dict_of_objects: dict, plain: bool=True):
        def update_plot(val):
            title_str = f"{forcing.name()}: {str(forcing.time()[val])}"
            vmax = np.max(forcing.magnitude())
            figure_dict = basic_funcs.plot_field(forcing.lon(), forcing.lat(), forcing.u()[val,:,:], forcing.v()[val,:,:],
                            var='ff', title_str=title_str, fig=self.fig, ax=self.ax, cbar=self._init_plot, vmin=0, vmax=vmax, barbs=True)
            self.fig = figure_dict.get('fig')
            self.ax = figure_dict.get('ax')
            self._init_plot = False

        forcing = dict_of_objects['Forcing']
        self.fig, self.ax = plt.subplots()

        ax_slider = plt.axes([0.17, 0.05, 0.65, 0.03])
        time_slider = Slider(ax_slider, 'time_index', 0, len(forcing.time())-1, valinit=0, valstep=1)
        time_slider.on_changed(update_plot)
        self._init_plot = True
        update_plot(0)
        #axnodenr = plt.axes([0.02,0.9,0.2,0.03])
        #giver    = Button(axnodenr,'Next')
        #giver.on_clicked(update_plot)
        plt.show(block=True)

        return self.fig

class OldTopoPlotter(GridPlotter):

    def _extension(self):
        return 'pdf'

    def grid(self, dict_of_objects: dict, plain: bool=True):
        """Creates a plot of the topography when a Grid-object is provided.

        Options

        plain:          True / False(default). Hint to plotter how much detail
                        the user wants into the plot.
        """

        grid = dict_of_objects['Grid']
        boundary = dict_of_objects['Boundary']
        forcing = dict_of_objects['Forcing']

        # Create basic plot
        fig = plt.figure()
        levels = np.linspace(0, np.max(grid.topo()), 20, endpoint=True)

        if grid.structured():
            plt.contourf(grid.lon(),grid.lat(),grid.topo(),levels)
        else:
            triang = mtri.Triangulation(grid.lon(), grid.lat(), triangles = grid.tri())
            topo=grid.topo()
            topo[np.isnan(topo)]=0 # Land to zero
            plt.tricontourf(triang, topo, levels)
            plt.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')

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



        return fig

    def topo(self, dict_of_objects: dict, plain: bool=True):
        """Creates a plot of the topography when a Grid-object is provided.

        Options

        plain:          True / False(default). Hint to plotter how much detail
                        the user wants into the plot.
        """

        grid = dict_of_objects['Grid']
        boundary = dict_of_objects['Boundary']
        forcing = dict_of_objects['Forcing']

        # Create basic plot
        fig = plt.figure()
        levels = np.linspace(0, max(grid.raw_topo()), 20, endpoint=True)
        mask = ~np.isnan(grid.raw_topo())
        #plt.tricontourf(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask],levels)
        msg.plain('Might take some time to create irregular plot...')
        #plt.tripcolor(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask])
        plt.tricontourf(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask], levels=levels)
        plt.plot(grid.raw_lon()[~mask],grid.raw_lat()[~mask], 'w.', markersize=0.5)

        # Plot grid edges


        x0=min(grid.lon())
        x1=max(grid.lon())
        y0=min(grid.lat())
        y1=max(grid.lat())

        x=[x0,x0,x1,x1,x0]
        y=[y0,y1,y1,y0,y0]

        plt.plot(x,y,'k', label='Grid')

        if len(grid.lon())>2:
            if grid.structured():
                lonQ, latQ = np.meshgrid(grid.lon(), grid.lat())
                lonQ=lonQ.ravel()
                latQ=latQ.ravel()
            elif grid.tri() is None: # Unstruct without triangulation
                plt.scatter(grid.lon(), grid.lat(),2,'k', label='Grid points')
            else:
                plt.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')

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
        #cbar.set_clim(0, 1500)
        cbar.set_label('Depth (m)', rotation=90)
        plt.title(f"{grid.name()} topography")



        return fig


class MaskPlotter(GridPlotter):
    def _extension(self):
        return 'pdf'

    def grid(self, grid: Grid, forcing: Forcing=None, boundary: Boundary=None, filename: str='', plain: bool=True) -> Tuple:
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

        return fig, file_module.add_suffix(filename, 'mask')
