from dnora.modelrun import ModelRun
import pandas as pd


def test_normal_init():
    t0 = "2020-01-01 00:00"
    t1 = "2020-01-02 00:00"
    model = ModelRun(start_time=t0, end_time=t1)

    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())


def test_forecast_init():
    t0 = "2020-01-01 00:00"
    t1 = "2020-01-01 10:00"
    model = ModelRun()
    model.activate_forecast_mode(reference_time=t0, forecast_length=10)

    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())


def test_year_init():
    t0 = "2020-01-01 00:00"
    t1 = "2020-12-31 23:00"
    model = ModelRun(year=2020)
    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())


def test_month_init():
    t0 = "2012-02-01 00:00"
    t1 = "2012-02-29 23:00"
    model = ModelRun(year=2012, month=2)
    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())


def test_year_init_hotstart():
    t0 = "2020-01-01 00:00"
    t1 = "2021-01-01 00:00"
    model = ModelRun(year=2020, hotstart_hour=True)
    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())


def test_month_init_hotstart():
    t0 = "2012-02-01 00:00"
    t1 = "2012-03-01 00:00"
    model = ModelRun(year=2012, month=2, hotstart_hour=True)
    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())
