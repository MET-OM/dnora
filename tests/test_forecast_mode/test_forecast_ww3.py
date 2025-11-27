import dnora as dn
from dnora.utils.io import read_ww3_nml
import shutil
import os
def cleanup():
    if os.path.isdir("forecast_WW3"):
        shutil.rmtree("forecast_WW3", onerror=handle_remove_readonly)

def test_ww3_restart_no_stride():
    """Should request restart only at end of forecast"""
    cleanup()
    grid = dn.grid.Grid(lon=(0,10), lat=(50,60), name= 'forecast')
    grid.set_spacing(dlon=1, dlat=1)
    grid.mesh_grid()

    t0 = "2020-01-01 00:00"
    t1 = "2020-01-01 10:00"
    model = dn.modelrun.Constant(grid)
    model.activate_forecast_mode(reference_time=t0, forecast_length=48)
    model.import_wind(u=1, v=2)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("forecast_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == "'3600'"
    cleanup()
    
def test_ww3_restart_stride_match_length():
    """Should request restart only at end of forecast"""
    cleanup()
    grid = dn.grid.Grid(lon=(0,10), lat=(50,60), name= 'forecast')
    grid.set_spacing(dlon=1, dlat=1)
    grid.mesh_grid()

    t0 = "2020-01-01 00:00"
    t1 = "2020-01-01 10:00"
    model = dn.modelrun.Constant(grid)
    model.activate_forecast_mode(reference_time=t0, forecast_length=48, stride=48)
    model.import_wind(u=1, v=2)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("forecast_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == "'3600'"
    cleanup()

def test_ww3_restart_stride_half_length():
    """Should request restart only at end of forecast"""
    cleanup()
    grid = dn.grid.Grid(lon=(0,10), lat=(50,60), name= 'forecast')
    grid.set_spacing(dlon=1, dlat=1)
    grid.mesh_grid()

    t0 = "2020-01-01 00:00"
    t1 = "2020-01-01 10:00"
    model = dn.modelrun.Constant(grid)
    model.activate_forecast_mode(reference_time=t0, forecast_length=48, stride=24)
    model.import_wind(u=1, v=2)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("forecast_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200101000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == f"'{str(3600*24)}'"
    cleanup()