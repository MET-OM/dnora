import pytest
import dnora as dn
from dnora.utils.io import read_ww3_nml
import os
import shutil
import stat
from pathlib import Path


def handle_remove_readonly(func, path, exc_info):
    """Clear the read-only bit and reattempt the removal."""
    os.chmod(path, stat.S_IWRITE)  # Change to writable
    func(path)


def cleanup():
    if os.path.isdir("TestGrid_WW3"):
        shutil.rmtree("TestGrid_WW3", onerror=handle_remove_readonly)


@pytest.fixture(scope="session")
def grid():
    grid = dn.grid.Constant(lon=(5, 6), lat=(70, 71), name="TestGrid")
    grid.set_spacing(dlon=0.2, dlat=0.1)
    grid.set_boundary_points(dn.grid.mask.Edges(["W", "N", "S"]))
    grid.import_topo(topo=55)
    grid.mesh_grid()
    return grid


def test_grid(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-31 00:00", end_time="2020-02-02 01:00"
    )
    exp = dn.export.WW3(model)
    exp.export_grid()

    exe = dn.executer.WW3(model)
    exe.write_grid_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_grid.nml")
    assert nml_dict["DEPTH_NML"]["DEPTH"]["FILENAME"] == "'TestGrid_bathy.txt'"
    assert nml_dict["MASK_NML"]["MASK"]["FILENAME"] == "'TestGrid_mapsta.txt'"

    exe.write_grid_file(folder_on_server="/lustre/folder/myfolder")
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_grid.nml")
    assert (
        Path(nml_dict["DEPTH_NML"]["DEPTH"]["FILENAME"]).as_posix()
        == "'/lustre/folder/myfolder/TestGrid_bathy.txt'"
    )
    assert (
        Path(nml_dict["MASK_NML"]["MASK"]["FILENAME"]).as_posix()
        == "'/lustre/folder/myfolder/TestGrid_mapsta.txt'"
    )
    cleanup()


def test_wind(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    exp = dn.export.WW3(model)
    exp.export_wind()

    exe = dn.executer.WW3(model)
    exe.write_wind_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_prnc_wind.nml")

    assert (
        nml_dict["FILE_NML"]["FILE"]["FILENAME"]
        == "'wind_ConstantData_TestGrid_20200101T0000_20200131T2300.nc'"
    )
    assert nml_dict["FILE_NML"]["FILE"]["LONGITUDE"] == "'lon'"
    assert nml_dict["FILE_NML"]["FILE"]["LATITUDE"] == "'lat'"
    assert nml_dict["FILE_NML"]["FILE"]["VAR(1)"] == "'u'"
    assert nml_dict["FILE_NML"]["FILE"]["VAR(2)"] == "'v'"

    exe.write_wind_file(folder_on_server="/lustre/folder/myfolder")
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_prnc_wind.nml")
    assert (
        Path(nml_dict["FILE_NML"]["FILE"]["FILENAME"]).as_posix()
        == "'/lustre/folder/myfolder/wind_ConstantData_TestGrid_20200101T0000_20200131T2300.nc'"
    )
    assert nml_dict["FILE_NML"]["FILE"]["LONGITUDE"] == "'lon'"
    assert nml_dict["FILE_NML"]["FILE"]["LATITUDE"] == "'lat'"
    assert nml_dict["FILE_NML"]["FILE"]["VAR(1)"] == "'u'"
    assert nml_dict["FILE_NML"]["FILE"]["VAR(2)"] == "'v'"
    cleanup()


def test_spectra(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_spectra()
    exp = dn.export.WW3(model)
    exp.export_spectra()

    exe = dn.executer.WW3(model)
    exe.write_spectra_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_bounc.nml")
    assert nml_dict["BOUND_NML"]["BOUND"]["MODE"] == "'WRITE'"
    assert nml_dict["BOUND_NML"]["BOUND"]["INTERP"] == "1"
    assert nml_dict["BOUND_NML"]["BOUND"]["VERBOSE"] == "2"
    assert nml_dict["BOUND_NML"]["BOUND"]["FILE"] == "'spectral_boundary_files.list'"
    with open("TestGrid_WW3/spectral_boundary_files.list", "r") as file:
        assert file.readline()[0:3] == "ww3"

    exe.write_spectra_file(
        folder_on_server="/lustre/folder/myfolder", verbose_level=1, method="linear"
    )
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_bounc.nml")
    assert nml_dict["BOUND_NML"]["BOUND"]["MODE"] == "'WRITE'"
    assert nml_dict["BOUND_NML"]["BOUND"]["INTERP"] == "2"
    assert nml_dict["BOUND_NML"]["BOUND"]["VERBOSE"] == "1"
    assert nml_dict["BOUND_NML"]["BOUND"]["FILE"] == "'spectral_boundary_files.list'"
    with open("TestGrid_WW3/spectral_boundary_files.list", "r") as file:
        assert Path(file.readline()[0:7]).as_posix() == "/lustre"

    cleanup()


def test_shel(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()

    exe = dn.executer.WW3(model)
    exe.write_input_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'F'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'F'"
    cleanup()


def test_shel_all_forcing(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    model.import_waterlevel()
    model.import_current()

    exe = dn.executer.WW3(model)
    exe.write_input_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'T'"
    assert (
        nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["POINT"]["FILE"] == "'spectral_points.list'"
    )
    assert len(nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["FIELD"]["LIST"]) > 0
    cleanup()


def test_shel_all_forcing(grid):
    cleanup()
    grid.set_output_points(dn.grid.mask.All())
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    model.import_waterlevel()
    model.import_current()

    exe = dn.executer.WW3(model)
    exe.write_input_file()
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'T'"
    assert (
        nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["POINT"]["FILE"] == "'spectral_points.list'"
    )
    assert len(nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["FIELD"]["LIST"]) > 0
    cleanup()


def test_shel_homog_wind(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_waterlevel()

    exe = dn.executer.WW3(model)
    exe.write_input_file(homog={"wind": (5, 6)})
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'H'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'F'"
    assert len(nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["FIELD"]["LIST"]) > 0

    assert nml_dict["HOMOG_COUNT_NML"]["HOMOG_COUNT"]["N_WND"] == "1"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["NAME"] == "'WND'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["DATE"] == "'20200130000000'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["VALUE1"] == "5"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["VALUE2"] == "6"
    cleanup()


def test_shel_homog_waterlevel(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    exe = dn.executer.WW3(model)
    exe.write_input_file(homog={"waterlevel": 1})
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'H'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'F'"
    assert len(nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["FIELD"]["LIST"]) > 0

    assert nml_dict["HOMOG_COUNT_NML"]["HOMOG_COUNT"]["N_LEV"] == "1"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["NAME"] == "'LEV'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["DATE"] == "'20200130000000'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["VALUE1"] == "1"
    cleanup()


def test_shel_homog_current(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    exe = dn.executer.WW3(model)
    exe.write_input_file(homog={"current": (0, 3)})
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_shel.nml")
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["START"] == "'20200130000000'"
    assert nml_dict["DOMAIN_NML"]["DOMAIN"]["STOP"] == "'20200131230000'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WINDS"] == "'T'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["WATER_LEVELS"] == "'F'"
    assert nml_dict["INPUT_NML"]["INPUT"]["FORCING"]["CURRENTS"] == "'H'"
    assert len(nml_dict["OUTPUT_TYPE_NML"]["TYPE"]["FIELD"]["LIST"]) > 0

    assert nml_dict["HOMOG_COUNT_NML"]["HOMOG_COUNT"]["N_CUR"] == "1"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["NAME"] == "'CUR'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["DATE"] == "'20200130000000'"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["VALUE1"] == "0"
    assert nml_dict["HOMOG_INPUT_NML"]["HOMOG_INPUT(1)"]["VALUE2"] == "3"
    cleanup()


def test_ounp(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    exe = dn.executer.WW3(model)
    exe.write_input_file(homog={"wind": (0, 3)})
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_ounp.nml")
    assert nml_dict["POINT_NML"]["POINT"]["TIMESTART"] == "'20200130000000'"
    assert nml_dict["POINT_NML"]["POINT"]["TIMECOUNT"] == "'48'"
    assert nml_dict["POINT_NML"]["POINT"]["TIMESTRIDE"] == "'3600'"
    assert nml_dict["POINT_NML"]["POINT"]["TIMESPLIT"] == "6"

    assert nml_dict["FILE_NML"]["FILE"]["NETCDF"] == "3"
    assert nml_dict["SPECTRA_NML"]["SPECTRA"]["OUTPUT"] == "3"
    assert nml_dict["PARAM_NML"]["PARAM"]["OUTPUT"] == "6"
    cleanup()


def test_ounf(grid):
    cleanup()
    model = dn.modelrun.Constant(
        grid, start_time="2020-01-30 00:00", end_time="2020-01-31 23:00"
    )
    model.import_wind()
    exe = dn.executer.WW3(model)
    exe.write_input_file(homog={"wind": (0, 3)})
    nml_dict = read_ww3_nml("TestGrid_WW3/ww3_ounf.nml")
    assert nml_dict["FIELD_NML"]["FIELD"]["TIMESTART"] == "'20200130000000'"
    assert nml_dict["FIELD_NML"]["FIELD"]["TIMECOUNT"] == "'48'"
    assert nml_dict["FIELD_NML"]["FIELD"]["TIMESTRIDE"] == "'3600'"
    assert nml_dict["FIELD_NML"]["FIELD"]["TIMESPLIT"] == "6"
    assert nml_dict["FIELD_NML"]["FIELD"][
        "LIST"
    ] == "'HS LM TP DIR SPR DP T02 T0M1 T01 UST CHA DPT WND USS TUS TAW TWO TOC FAW FOC PHS PTP PTM10 PT01 PT02 PDIR PDP MXE MXH MXHC SDMH SDMHC ABR UBR FBB TBB CGE WCC WBT'".replace(
        " ", ""
    )
    assert nml_dict["FILE_NML"]["FILE"]["NETCDF"] == "3"
    cleanup()
