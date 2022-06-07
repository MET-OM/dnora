import matplotlib.pyplot as plt
from abc import ABC, abstractmethod
import numpy as np
from ..grd.grd_mod import Grid, GridProcessor
from ..grd.process import TrivialFilter
from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing
from .. import file_module
from .. import msg
from .. import aux
from typing import Tuple
import matplotlib.tri as mtri
from . import basic_funcs
import cmocean.cm
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
from matplotlib.colors import LinearSegmentedColormap
class GridPlotter(ABC):
    """Plots data from Grid-object."""

    @abstractmethod
    def _extension(self) -> str:
        pass

    def grid(self, dict_of_objects: dict, plain: bool) -> dict:
        msg.warning('Plotting of meshed data not implemented!')
        return None

    def topo(self, dict_of_objects: dict, plain: bool) -> dict:
        msg.warning('Plotting of xyz data not implemented!')
        return None

class TopoPlotter(GridPlotter):
    def _extension(self):
        return 'pdf'

    def grid(self, dict_of_objects: dict, plain: bool=True):
        grid = dict_of_objects['Grid']
        boundary = dict_of_objects['Boundary']
        forcing = dict_of_objects['Forcing']
        #
        figure_dict = basic_funcs.plot_field(grid.lon(), grid.lat(), grid.topo(land=np.nan), var='topo', title_str=f"{grid.name()}", cbar=True)
        fig = figure_dict.get('fig')
        ax = figure_dict.get('ax')
        cbar = figure_dict.get('cbar')
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
        x0 = np.min([min(grid.lon()), min(forcing.lon()), min(boundary.lon())])
        x1 = np.max([max(grid.lon()), max(forcing.lon()), max(boundary.lon())])
        y0 = np.min([min(grid.lat()), min(forcing.lat()), min(boundary.lat())])
        y1 = np.max([max(grid.lat()), max(forcing.lat()), max(boundary.lat())])

        x0, x1, y0, y1 = aux.expand_area(x0, x1, y0, y1, 1.2)

        ax.set_xlim([x0,x1])
        ax.set_ylim([y0,y1])
        return {'fig': fig, 'ax': ax, 'cbar': cbar}


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
        fig, ax = plt.subplots()
        levels = np.linspace(0, max(grid.raw_topo()), 20, endpoint=True)
        mask = ~np.isnan(grid.raw_topo())
        #plt.tricontourf(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask],levels)
        msg.plain('Might take some time to create irregular plot...')
        #plt.tripcolor(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask])
        cont = ax.tricontourf(grid.raw_lon()[mask],grid.raw_lat()[mask],grid.raw_topo()[mask], levels=levels)
        ax.plot(grid.raw_lon()[~mask],grid.raw_lat()[~mask], 'w.', markersize=0.5)

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
                ax.scatter(grid.lon(), grid.lat(),2,'k', label='Grid points')
            else:
                ax.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')

        # Plot boundary points if they exist
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

        ax.legend(loc="upper right")
        cbar = fig.colorbar(cont)
        #cbar.set_clim(0, 1500)
        cbar.set_label('Depth (m)', rotation=90)
        ax.set(title=f"{grid.name()} topography")

        return {'fig': fig, 'ax': ax, 'cont': cont, 'cbar': cbar}


class ForcingPlotter(GridPlotter):
    def _extension(self):
        return 'pdf'

    def grid(self, dict_of_objects: dict, plain: bool=True):
        def update_plot(val):
            title_str = f"{forcing.name()}: {str(forcing.time()[val])}"
            vmax = np.max(forcing.magnitude())
            self.fig_dict = basic_funcs.plot_field(forcing.lon(), forcing.lat(), forcing.u()[val,:,:], forcing.v()[val,:,:],
                            var='ff', title_str=title_str, fig_dict=self.fig_dict, cbar=self._init_plot, vmin=0, vmax=vmax, barbs=False)
            # self.fig = fig_dict.get('fig')
            # self.ax = fig_dict.get('ax')
            self._init_plot = False

        forcing = dict_of_objects['Forcing']
        fig, ax = plt.subplots()
        self.fig_dict = {'fig': fig, 'ax': ax}
        ax_slider = plt.axes([0.17, 0.05, 0.65, 0.03])
        time_slider = Slider(ax_slider, 'time_index', 0, len(forcing.time())-1, valinit=0, valstep=1)
        time_slider.on_changed(update_plot)
        self._init_plot = True
        update_plot(0)
        #axnodenr = plt.axes([0.02,0.9,0.2,0.03])
        #giver    = Button(axnodenr,'Next')
        #giver.on_clicked(update_plot)
        plt.show(block=True)

        return self.fig_dict

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
        fig, ax = plt.subplots()


        if grid.structured():
            levels = np.linspace(0, np.max(grid.topo()), 20, endpoint=True)
            cont = ax.contourf(grid.lon(),grid.lat(),grid.topo(),levels)
        else:

            levels = np.linspace(0, np.max(grid.topo()), 20, endpoint=True)
            triang = mtri.Triangulation(grid.lon(), grid.lat(), triangles = grid.tri())
            topo=grid.topo()
            topo[np.isnan(topo)]=0 # Land to zero
            cmap = cmocean.cm.topo_r
            colors = cmap(np.linspace(0.5, 1, cmap.N // 2))
            cmap2 = LinearSegmentedColormap.from_list('Upper Half', colors)
            cont = ax.tricontourf(triang, topo, levels, cmap=cmap2)
            ax.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')

        # Plot boundary points if they exist
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

        ax.legend(loc="upper left")
        cbar = fig.colorbar(cont)
        cbar.set_label('Depth (m)', rotation=90)
        #ax.title(f"{grid.name()} topography")



        return {'fig': fig, 'ax': ax, 'cbar': cbar}

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
