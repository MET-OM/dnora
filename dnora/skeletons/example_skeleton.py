from skeletons import Skeleton, PointSkeleton, GriddedSkeleton
import numpy as np
import xarray as xr

coords_dict = {'x': np.arange(0,10), 'y': np.arange(11,20)}
#vars_dict = {'topo': (['points'], topo), 'lon': (['points'], lon), 'lat': (['points'], lat)}
data = xr.Dataset(coords=coords_dict)

sktn = Skeleton()
sktn.data = data

psktn = PointSkeleton(lon=np.arange(10,20), lat=np.arange(60,70))
#coords_dict = {'station': np.arange(0,10)}
#vars_dict = {'x': (['station'], np.arange(10,20)), 'y': (['station'], np.arange(20,30))}
#vars_dict = {'topo': (['points'], topo), 'lon': (['points'], lon), 'lat': (['points'], lat)}
#pdata = xr.Dataset(coords=coords_dict, data_vars=vars_dict)


gsktn = GriddedSkeleton()
gdata = gsktn._create_xr(x=np.arange(0,10), y=np.arange(10,20))
gsktn.data = gdata
