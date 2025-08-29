from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import add_frequency, add_direction


@add_direction(grid_coord=True)
@add_frequency(grid_coord=True)
class SpectralGrid(PointSkeleton):
    def __init__(self, **kwargs):
        super().__init__(x=0, y=0, **kwargs)

    def __repr__(self) -> str:
        string = "Spectral grid\n"
        string += f"freq: {self.freq()}\n"
        string += f"dirs: {self.dirs()}"
        return string
