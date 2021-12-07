from __future__ import annotations

from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .trg_mod import TrGrid

class TrGridPlotter(ABC):
    """Plots data from Grid-object."""

    @abstractmethod
    def __call__(self, grid: Grid, filename: str) -> Tuple:


        return fig, filename


class TriPlotter(TrGridPlotter):
    def __call__(self, grid: TrGrid, filename: str='') -> Tuple:
        fig = plt.figure()
        plt.triplot(grid.lon(), grid.lat(), triangles=grid.tri(), linewidth=0.2, color='black')
        plt.plot(grid.lon()[grid.boundary()-1],grid.lat()[grid.boundary()-1],'rx')

        return fig, filename
