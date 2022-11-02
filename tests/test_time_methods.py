from dnora.grd import Grid
from dnora.mdl import ModelRun
from dnora.wnd.read import ConstantForcing
from dnora.wnd.write import Null
import pandas as pd
def test_one_time():
    grid = Grid(lon=(5,6), lat=(60,61))
    start_time = '2020-01-31 00:00:00'
    end_time = '2020-01-31 00:00:00'
    model = ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    model.import_forcing(ConstantForcing())
    time = pd.date_range(start=start_time, end=end_time, freq='H')

    assert time[0] == model.forcing().time()[0]
    assert time[0] == model.forcing().time()[1]

    assert model.forcing().years(datetime=False) == ['2020']
    assert model.forcing().months(datetime=False) == ['2020-01']
    assert model.forcing().days(datetime=False) == ['2020-01-31']

    pd.testing.assert_index_equal(model.forcing().years(), pd.to_datetime(['2020-01-01']))
    pd.testing.assert_index_equal(model.forcing().months(), pd.to_datetime(['2020-01-01']))
    pd.testing.assert_index_equal(model.forcing().days(), pd.to_datetime(['2020-01-31']))

    model.export_forcing(Null())
def test_time_over_month():
    grid = Grid(lon=(5,6), lat=(60,61))
    start_time = '2020-01-31 05:00:00'
    end_time = '2020-02-01 05:00:00'
    model = ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    model.import_forcing(ConstantForcing())
    time = pd.date_range(start=start_time, end=end_time, freq='H')

    pd.testing.assert_index_equal(time, model.forcing().time())

    assert model.forcing().years(datetime=False) == ['2020']
    assert model.forcing().months(datetime=False) == ['2020-01', '2020-02']
    assert model.forcing().days(datetime=False) == ['2020-01-31', '2020-02-01']

    pd.testing.assert_index_equal(model.forcing().years(), pd.to_datetime(['2020-01-01']))
    pd.testing.assert_index_equal(model.forcing().months(), pd.to_datetime(['2020-01-01', '2020-02-01']))
    pd.testing.assert_index_equal(model.forcing().days(), pd.to_datetime(['2020-01-31', '2020-02-01']))

def test_time_over_year():
    grid = Grid(lon=(5,6), lat=(60,61))
    start_time = '2020-12-31 05:00:00'
    end_time = '2021-02-01 05:00:00'
    model = ModelRun(grid=grid, start_time=start_time, end_time=end_time)
    model.import_forcing(ConstantForcing())
    time = pd.date_range(start=start_time, end=end_time, freq='H')

    pd.testing.assert_index_equal(time, model.forcing().time())

    assert model.forcing().years(datetime=False) == ['2020', '2021']
    assert model.forcing().months(datetime=False) == ['2020-12', '2021-01', '2021-02']
    assert model.forcing().days(datetime=False) == ['2020-12-31']+[f'2021-01-{dd+1:02.0f}' for dd in range(31)]+['2021-02-01']

    pd.testing.assert_index_equal(model.forcing().years(), pd.to_datetime(['2020-01-01', '2021-01-01']))
    pd.testing.assert_index_equal(model.forcing().months(), pd.to_datetime(['2020-12-01', '2021-01-01', '2021-02-01']))
    pd.testing.assert_index_equal(model.forcing().days(), pd.to_datetime(['2020-12-31']+[f'2021-01-{dd+1:02.0f}' for dd in range(31)]+['2021-02-01']))
