from dnora.modelrun import ModelRun
import pandas as pd

def test_forecast_init():
    t0 = "2020-01-01 00:00"
    t1 = "2020-01-01 10:00"
    model = ModelRun()
    model.activate_forecast_mode(reference_time=t0, forecast_length=10)
    assert model.forecast_mode()

    pd.testing.assert_index_equal(pd.date_range(t0, t1, freq="h"), model.time())

def test_set_stride():
    t0 = "2020-01-01 00:00"
    t1 = "2020-01-02 00:00"
    model = ModelRun()
    model.activate_forecast_mode(reference_time=t0, forecast_length=48, stride=24)

    assert model._stride == 24
    model.deactivate_forecast_mode()
    assert model._stride is None
    assert not model.forecast_mode()