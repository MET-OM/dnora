# =============================================================================
# Read in a constant wind speed for an area
# =============================================================================
from dnora import grid, modelrun

area = grid.Grid(lon=(4.00, 11.0), lat=(62.0, 65.0), name="sula_trondheim")
model = modelrun.ERA5(area, start_time="2021-08-25T00:00", end_time="2021-08-25T05:00")
model.import_wind()
print(model.wind())
