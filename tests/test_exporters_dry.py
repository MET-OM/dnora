import dnora as dn
import pytest


@pytest.fixture(scope="session")
def model():
    # grid = grd.Grid(lon=(10, 20), lat=(60, 65))
    area = dn.grid.Grid(lon=5, lat=60)
    model = dn.modelrun.ModelRun(grid=area, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )
    return model


def test_dnora(model):
    exporter = dn.export.DataExporter(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_cache(model):
    exporter = dn.export.Cacher(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_swan(model):
    exporter = dn.export.SWAN(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_ww3(model):
    exporter = dn.export.WW3(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_hos(model):
    exporter = dn.export.HOS_Ocean(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_reef3d(model):
    exporter = dn.export.REEF3D(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()


def test_swash(model):
    exporter = dn.export.SWASH(model)
    exporter.export_spectra()
    exporter.export_wind()
    exporter.export_spectra1d()
    exporter.export_waveseries()
