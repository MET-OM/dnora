# =============================================================================
# Read in a constant wind speed for an area
# =============================================================================
from dnora import grid, modelrun, wind

area = grid.Grid(lon=(4.00, 11.0), lat=(62.0, 65.0), name="sula_trondheim")
model = modelrun.NORA3(area, start_time="2021-08-25T00:00", end_time="2021-08-25T05:00")
model.import_wind()
print(model.wind())


model2 = modelrun.ModelRun(
    area, start_time="2021-08-25T00:00", end_time="2021-08-25T05:00"
)
model2.import_wind(wind.read_metno.MyWave3km())
print(model2.wind())

model3 = modelrun.WAM4km(
    area, start_time="2021-08-25T00:00", end_time="2021-08-25T05:00"
)
model3.import_wind()
print(model3.wind())
