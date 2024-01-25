from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_datavar, add_time

from dnora.metaparameter.parameters import SeaLevel


@add_datavar(name="eta", default_value=0.0)
@add_time(grid_coord=True)
class WaterLevel(GriddedSkeleton):
    meta_dict = {"eta": SeaLevel}
    pass
