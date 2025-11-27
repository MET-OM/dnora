import dnora as dn
from dnora.utils.io import read_ww3_nml
import shutil
import os
import stat
def handle_remove_readonly(func, path, exc_info):
    """Clear the read-only bit and reattempt the removal."""
    os.chmod(path, stat.S_IWRITE)  # Change to writable
    func(path)


def cleanup():
    if os.path.isdir("GridName_WW3"):
        shutil.rmtree("GridName_WW3", onerror=handle_remove_readonly)

def test_ww3_restart_no_stride():
    """Should request restart only at end of forecast"""
    cleanup()
    model = dn.modelrun.ModelRun()
    model.activate_forecast_mode(reference_time="2020-01-01 00:00", forecast_length=48)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("GridName_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == "'3600'"
    cleanup()

def test_ww3_restart_stride_match_length():
    """Should request restart only at end of forecast"""
    cleanup()
    model = dn.modelrun.ModelRun()
    model.activate_forecast_mode(reference_time="2020-01-01 00:00", forecast_length=48, stride=48)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("GridName_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == f"'{str(3600*48)}'"
    cleanup()

def test_ww3_restart_stride_half_length():
    """Should request restart only at end of forecast"""
    cleanup()
    model = dn.modelrun.ModelRun()
    model.activate_forecast_mode(reference_time="2020-01-01 00:00", forecast_length=48, stride=24)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("GridName_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20200102000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20200103000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == f"'{str(3600*24)}'"
    cleanup()

def test_ww3_restart_destine_type():
    """Should request first restart 24 hours from start """
    cleanup()
    model = dn.modelrun.ModelRun()
    model.activate_forecast_mode(reference_time='20251125 000000', forecast_length=47, stride=24)

    exe = dn.executer.WW3(model)
    exe.write_input_file()

    nml_dict = read_ww3_nml("GridName_WW3/ww3_shel.nml")
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['START'] == "'20251126000000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STOP'] == "'20251126230000'"
    assert nml_dict['OUTPUT_DATE_NML']['DATE']['RESTART']['STRIDE'] == f"'{str(3600*24)}'"
    cleanup()