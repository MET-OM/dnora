from __future__ import annotations

from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .trg_mod import Grid

class TrGridPlotter(ABC):
    """Plots data from Grid-object."""

    @abstractmethod
    def __call__(self, grid: Grid, filename: str) -> Tuple:


        return fig, filename


class TriPlotter(TrGridPlotter):
    def __call__(self, grid: Grid, filename: str='') -> Tuple:
        fig = plt.figure()
        plt.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')
        plt.plot(grid.lon()[grid.boundary_inds()],grid.lat()[grid.boundary_inds()],'rx')

        return fig, filename

class TriTopoPlotter(TrGridPlotter):
    def __call__(self, grid: Grid, filename: str='') -> Tuple:
        fig = plt.figure()
        triang = mtri.Triangulation(grid.lon(), grid.lat(), triangles = grid.tri())
        levels = np.linspace(0, np.max(grid.topo()), 100, endpoint=True)
        plt.tricontourf(triang, grid.topo(), levels)
        cbar = plt.colorbar()

        plt.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')
        plt.plot(grid.lon()[grid.boundary_inds()],grid.lat()[grid.boundary_inds()],'rx')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        cbar.set_label('Depth (m)', rotation=90)
        return fig, filename
