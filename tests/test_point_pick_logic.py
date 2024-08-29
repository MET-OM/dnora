import dnora as dn
from dnora.read.abstract_readers import SpectralDataReader
from dnora.type_manager.dnora_types import DnoraDataType
from dnora.type_manager.data_sources import DataSource
import pandas as pd
import numpy as np


class DummyReader(SpectralDataReader):
    _lon = np.arange(0, 10)
    _lat = np.arange(50, 60)

    def default_data_source(self) -> DataSource:
        return DataSource.CREATION

    def get_coordinates(self, grid, start_time, source, folder, **kwargs):
        return {"lon": self._lon, "lat": self._lat}

    def __call__(
        self,
        obj_type: DnoraDataType,
        grid: dn.grid.Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        inds: np.ndarray,
        **kwargs,
    ):
        freq = np.arange(0, 1, 0.1)
        dirs = np.arange(0, 360, 10)
        times = pd.date_range(start_time, end_time, freq="1h")
        spec = np.zeros((len(times), len(inds), len(freq), len(dirs)))
        coord_dict = {
            "lon": self._lon[inds],
            "lat": self._lat[inds],
            "time": times,
            "freq": freq,
            "dirs": dirs,
        }
        data_dict = {"spec": spec}

        return coord_dict, data_dict, {"meta": "test_meta"}


def test_no_grid_leads_to_trivial_picker():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    model = dn.modelrun.ModelRun(start_time=start_time, end_time=end_time)
    assert isinstance(model._point_picker, dn.pick.Trivial)
    model.import_spectra(DummyReader())
    np.testing.assert_array_equal(np.arange(0, 10), model.spectra().inds())


def test_single_point_grid_uses_nearest_picker():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    grid = dn.grid.Grid(lon=1.9, lat=52.2)
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_spectra(DummyReader())
    np.testing.assert_almost_equal(model.spectra().lon(), 2)
    np.testing.assert_almost_equal(model.spectra().lat(), 52)


def test_single_point_grid_uses_nearest_picker_but_user_overrides():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    grid = dn.grid.Grid(lon=1.9, lat=52.2)
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_spectra(DummyReader(), point_picker=dn.pick.Trivial())
    np.testing.assert_array_equal(np.arange(0, 10), model.spectra().inds())


def test_empty_area_uses_area_picker():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    grid = dn.grid.Grid(lon=(4, 5), lat=(54, 55))
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_spectra(DummyReader(), expansion_factor=1)
    np.testing.assert_almost_equal(model.spectra().lon(), np.array([4, 5]))
    np.testing.assert_almost_equal(model.spectra().lat(), np.array([54, 55]))


def test_empty_area_uses_area_picker_but_user_overrides():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    grid = dn.grid.Grid(lon=(4, 5), lat=(54, 55))
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_spectra(DummyReader(), point_picker=dn.pick.NearestGridPoint())
    assert model.spectra() is None


def test_empty_area_uses_area_picker_but_user_sets_points_and_overrides():
    """If no grid is given, we will read all the data points in the source"""
    start_time = "2020-01-31 00:00:00"
    end_time = "2020-02-01 00:00:00"
    grid = dn.grid.Grid(lon=(4, 5), lat=(54, 55))
    grid.set_boundary_points(dn.grid.mask.LonLat(lon=(4), lat=54))
    model = dn.modelrun.ModelRun(grid, start_time=start_time, end_time=end_time)
    model.import_spectra(DummyReader(), point_picker=dn.pick.NearestGridPoint())
    np.testing.assert_almost_equal(model.spectra().lon(), np.array([4]))
    np.testing.assert_almost_equal(model.spectra().lat(), np.array([54]))
