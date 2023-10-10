from dnora import grd, wnd, bnd, spc, wsr, mdl, pick, exp, run
import pytest


@pytest.fixture(scope="session")
def model():
    grid = grd.Grid(lon=(10, 20), lat=(60, 65))
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid=grid, dry_run=True)

    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(
        spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker()
    )
    model.spectra_to_waveseries()

    return model


def test_dnora(model):
    exporter = exp.DataExporter(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()


def test_swan(model):
    exporter = exp.SWAN(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.write_input_file()

    exporter.cache_boundary()
    exporter.cache_forcing()
    exporter.cache_spectra()
    exporter.cache_waveseries()


def test_ww3(model):
    exporter = exp.WW3(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.write_input_file()


def test_ww3(model):
    exporter = exp.HOS_ocean(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.write_input_file()


def test_ww3(model):
    exporter = exp.REEF3D(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.write_input_file()


def test_ww3(model):
    exporter = exp.SWASH(model)
    exporter.export_boundary()
    exporter.export_forcing()
    exporter.export_spectra()
    exporter.export_waveseries()
    exporter.write_input_file()
