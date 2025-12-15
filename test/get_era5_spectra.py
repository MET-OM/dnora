import dnora as dn
import pandas as pd

boundaries = pd.read_csv('boundaries/cluster_CARRA2/boundary_cluster_0.txt',sep='\s+', names=['lon', 'lat'])
points = dn.grid.TriGrid(lon=boundaries.lon[:].to_list(), lat=boundaries.lat[:].to_list(), name='CARRA2WW')
points.set_boundary_points(dn.grid.mask.All())

# Option 1: Define start and end_time
model = dn.modelrun.ERA5(points, start_time='2017-01-01T00:00',end_time='2017-01-02T00:00')

# Option 2: Download for one day
#model = dn.modelrun.ERA5(points, year=2017,month=1, day=1) # 2017-01-01 00:00 - 2017-01-01 23:00
#model = dn.modelrun.ERA5(points, year=2017,month=1, day=1, hotstart_hour=True) # 2017-01-01 00:00 - 2017-01-02 00:00

# Option 3: Set start and number of hours
#model = dn.modelrun.ERA5(points)
#model.activate_forecast_mode(reference_time='2017-01-01T00:00', forecast_length=24) # 2017-01-01 00:00 - 2017-01-02 00:00


# post_process = False means that empty (possibly NaN) spectra are NOT removed
model.import_spectra(point_picker=dn.pick.NearestGridPoint(), post_process=False) 

exp = dn.export.WW3(model)
exp.export_spectra()