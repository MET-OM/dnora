from dnora.grid import Grid
from dnora import wind, spectra, spectra1d, waveseries, modelrun, pick, export, executer


def test_model_one_point():
    grid = Grid(lon=5, lat=60)
    model = modelrun.ModelRun(grid=grid, dry_run=True)

    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )

    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )

    DNORA = export.DataExporter(model)
    DNORA.export_spectra()
    DNORA.export_wind()
    DNORA.export_spectra1d()
    DNORA.export_waveseries()

    SWAN = export.SWAN(model)
    SWAN.export_spectra()
    SWAN.export_wind()
    SWAN.export_spectra1d()
    SWAN.export_waveseries()

    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_area():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    model = modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )
    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_gridded():
    grid = Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=10, ny=10)
    model = modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )

    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_point_cartesian():
    grid = Grid(x=5, y=3)
    model = modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )

    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_area_cartesian():
    grid = Grid(x=(5, 6), y=(3, 4))
    model = modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )
    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_gridded_cartesian():
    grid = Grid(x=(5, 6), y=(3, 4))
    grid.set_spacing(nx=10, ny=10)
    model = modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(spectra.read_metno.NORA3(), pick.Area())
    model.import_wind(wind.read_metno.NORA3())
    model.import_spectra1d(
        spectra1d.read.SpectraTo1D(model.spectra()), point_picker=pick.TrivialPicker()
    )
    model.import_waveseries(
        waveseries.read.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=pick.TrivialPicker(),
    )
    exe = executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)
