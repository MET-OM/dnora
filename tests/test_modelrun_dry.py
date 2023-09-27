from dnora import grd, wnd, bnd, spc, wsr, mdl, pick, exp, run


def test_model_one_point():
    grid = grd.Grid(lon=5, lat=60)
    model = mdl.ModelRun(grid=grid, dry_run=True)
    
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    
    

    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    DNORA = exp.DataExporter(model)
    DNORA.export_boundary()
    DNORA.export_forcing()
    DNORA.export_spectra()
    DNORA.export_waveseries()
    
    SWAN = exp.SWAN(model)
    SWAN.export_boundary()
    SWAN.export_forcing()
    SWAN.export_spectra()
    SWAN.export_waveseries()
    
        

    SWAN.write_input_file()

    model.run_model(run.SWAN())

    print(model)

def test_model_one_area():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    model = mdl.ModelRun(grid=grid, dry_run=True)
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    model.run_model(run.SWAN())
    print(model)

def test_model_gridded():
    grid = grd.Grid(lon=(5,6), lat=(60,61))
    grid.set_spacing(nx=10, ny=10)
    model = mdl.ModelRun(grid=grid, dry_run=True)
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    model.run_model(run.SWAN())
    print(model)

def test_model_one_point_cartesian():
    grid = grd.Grid(x=5, y=3)
    model = mdl.ModelRun(grid=grid, dry_run=True)
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    model.run_model(run.SWAN())
    print(model)

def test_model_one_area_cartesian():
    grid = grd.Grid(x=(5,6), y=(3,4))
    model = mdl.ModelRun(grid=grid, dry_run=True)
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    model.run_model(run.SWAN())
    print(model)

def test_model_gridded_cartesian():
    grid = grd.Grid(x=(5,6), y=(3,4))
    grid.set_spacing(nx=10, ny=10)
    model = mdl.ModelRun(grid=grid, dry_run=True)
    model.import_boundary(bnd.read_metno.NORA3(), pick.Area())
    model.import_forcing(wnd.read_metno.NORA3())
    model.import_spectra(spc.read.BoundaryToSpectra(model.boundary()), point_picker=pick.TrivialPicker())
    model.import_waveseries(wsr.read.SpectraToWaveSeries(model.spectra()), point_picker=pick.TrivialPicker())

    model.run_model(run.SWAN())
    print(model)
