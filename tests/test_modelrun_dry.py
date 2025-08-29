import dnora as dn


def test_model_one_point():
    grid = dn.grid.Grid(lon=5, lat=60)
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)

    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )

    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )

    generic = dn.export.DataExporter(model)
    generic.export_spectra()
    generic.export_wind()
    generic.export_spectra1d()
    generic.export_waveseries()

    SWAN = dn.export.SWAN(model)
    SWAN.export_spectra()
    SWAN.export_wind()
    SWAN.export_spectra1d()
    SWAN.export_waveseries()

    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_area():
    grid = dn.grid.Grid(lon=(5, 6), lat=(60, 61))
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )
    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_gridded():
    grid = dn.grid.Grid(lon=(5, 6), lat=(60, 61))
    grid.set_spacing(nx=10, ny=10)
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )

    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_point_cartesian():
    grid = dn.grid.Grid(x=5, y=3)
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )

    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_one_area_cartesian():
    grid = dn.grid.Grid(x=(5, 6), y=(3, 4))
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )
    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)


def test_model_gridded_cartesian():
    grid = dn.grid.Grid(x=(5, 6), y=(3, 4))
    grid.set_spacing(nx=10, ny=10)
    model = dn.modelrun.ModelRun(grid=grid, dry_run=True)
    model.import_spectra(dn.read.spectra.metno.NORA3(), dn.pick.Area())
    model.import_wind(dn.read.wind.metno.NORA3())
    model.import_spectra1d(
        dn.read.spectra1d.SpectraTo1D(model.spectra()), point_picker=dn.pick.Trivial()
    )
    model.import_waveseries(
        dn.read.waveseries.Spectra1DToWaveSeries(model.spectra1d()),
        point_picker=dn.pick.Trivial(),
    )
    exe = dn.executer.SWAN(model)
    exe.write_input_file()
    exe.run_model()

    print(model)
