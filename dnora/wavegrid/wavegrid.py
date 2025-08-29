from geo_skeletons import GriddedSkeleton
from geo_skeletons.decorators import add_mask, add_time, add_datavar, add_magnitude
import geo_parameters as gp


@add_magnitude(
    gp.ocean.Current, x="x_current", y="y_current", direction=gp.ocean.CurrentDir
)
@add_datavar(gp.ocean.YCurrent, default_value=0.1)
@add_datavar(gp.ocean.XCurrent, default_value=0.1)
@add_magnitude(gp.wind.Wind, x="x_wind", y="y_wind", direction=gp.wind.WindDir)
@add_datavar(gp.wind.YWind, default_value=0.1)
@add_datavar(gp.wind.XWind, default_value=0.1)
@add_datavar(gp.ocean.WaterDepth)
@add_datavar(gp.wave.Tm02Sea)
@add_datavar(gp.wave.Tm_10Sea)
@add_datavar(gp.wave.Tm01Sea)
@add_datavar(gp.wave.DirmSea)
@add_datavar(gp.wave.DirpSea)
@add_datavar(gp.wave.TpSea)
@add_datavar(gp.wave.HsSea)
@add_datavar(gp.wave.Tm02Swell)
@add_datavar(gp.wave.Tm_10Swell)
@add_datavar(gp.wave.Tm01Swell)
@add_datavar(gp.wave.DirmSwell)
@add_datavar(gp.wave.DirpSwell)
@add_datavar(gp.wave.TpSwell)
@add_datavar(gp.wave.HsSwell)
@add_datavar(gp.wave.Tm02)
@add_datavar(gp.wave.Tm_10)
@add_datavar(gp.wave.Tm01)
@add_datavar(gp.wave.Dirm)
@add_datavar(gp.wave.Dirp)
@add_datavar(gp.wave.Tp)
@add_datavar(gp.wave.Hs)
@add_mask(name="buoy", coord_group="grid", default_value=1)
@add_time(grid_coord=False)
class WaveGrid(GriddedSkeleton):
    pass
