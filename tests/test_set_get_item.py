import dnora as dn
from dnora.type_manager.dnora_types import DnoraDataType
import pytest


def test_get_grid():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    assert model.get(DnoraDataType.GRID) is not None
    assert model[DnoraDataType.GRID] is not None


def test_set_wind():
    model = dn.modelrun.ModelRun()
    assert model.wind() is None
    model[DnoraDataType.WIND] = "DummyWind"
    assert model.wind() is not None


def test_delete_grid_and_set_trigrid():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    del model[DnoraDataType.GRID]
    assert model.get(DnoraDataType.GRID) is None
    model[DnoraDataType.TRIGRID] = "DummyTriGrid"
    assert model.get(DnoraDataType.TRIGRID) is not None
    assert model[DnoraDataType.TRIGRID] is not None
    assert model.get(DnoraDataType.GRID) is None
    assert model.grid() == "DummyTriGrid"


def test_raises_keyerror():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    with pytest.raises(KeyError):
        model[DnoraDataType.WIND]
    assert model.get(DnoraDataType.WIND) is None


def test_get_grid_with_str():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    assert model.get("grid") is not None
    assert model["grid"] is not None


def test_set_wind_wit_str():
    model = dn.modelrun.ModelRun()
    assert model.wind() is None
    model["wind"] = "DummyWind"
    assert model.wind() is not None


def test_delete_grid_and_set_trigrid_with_str():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    del model["grid"]
    assert model.get("grid") is None
    model["trigrid"] = "DummyTriGrid"
    assert model.get("trigrid") is not None
    assert model["trigrid"] is not None
    assert model.get("grid") is None
    assert model.grid() == "DummyTriGrid"


def test_raises_keyerror_with_str():
    model = dn.modelrun.ModelRun()
    assert model.grid() is not None
    with pytest.raises(KeyError):
        model["wind"]
    assert model.get("wind") is None
