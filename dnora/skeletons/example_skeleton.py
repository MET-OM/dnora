from skeletons import Skeleton, PointSkeleton, GriddedSkeleton
import numpy as np
import xarray as xr

psktn = PointSkeleton(lon=np.arange(10,20), lat=np.arange(60,70))

gsktn = GriddedSkeleton(lon=np.arange(10,20), lat=np.arange(60,70))
