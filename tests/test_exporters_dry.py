from dnora import (
    grid,
    pick,
    spectra,
    wind,
    spectra1d,
    waveseries,
    export,
    modelrun,
)
import pytest


@pytest.fixture(scope="session")
def model():
    # grid = grd.Grid(lon=(10, 20), lat=(60, 65))
    area = grid.Grid(lon=5, lat=60)
    model = modelrun.ModelRun(grid=area, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.Trivial()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.Trivial(),
    )
    return model


def test_dnora(model):
    exporter = export.DataExporter(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_cache(model):
    exporter = export.Cacher(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_swan(model):
    exporter = export.SWAN(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_ww3(model):
    exporter = export.WW3(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_hos(model):
    exporter = export.HOS_Ocean(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_reef3d(model):
    exporter = export.REEF3D(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_swash(model):
    exporter = export.SWASH(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()
