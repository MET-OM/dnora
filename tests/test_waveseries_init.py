from dnora.grd import Grid
from dnora.wsr import WaveSeries
def test_waveseries_one_point():
    grid = Grid(lon=5, lat=60)
    waveseries = WaveSeries(grid=grid)
    print(waveseries)

def test_waveseries_one_area():
    grid = Grid(lon=(5,6), lat=(60,61))
    waveseries = WaveSeries(grid=grid)
    print(waveseries)

def test_waveseries_gridded():
    grid = Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=10, ny=10)
    waveseries = WaveSeries(grid=grid)
    print(waveseries)

def test_waveseries_one_point_cartesian():
    grid = Grid(x=5, y=3)
    waveseries = WaveSeries(grid=grid)
    print(waveseries)

def test_waveseries_one_area_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    waveseries = WaveSeries(grid=grid)
    print(waveseries)

def test_waveseries_gridded_cartesian():
    grid = Grid(x=(5,6), y=(3,4))
    grid.set_spacing(nx=10, ny=10)
    waveseries = WaveSeries(grid=grid)
    print(waveseries)
