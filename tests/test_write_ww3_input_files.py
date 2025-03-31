import pytest 
import dnora as dn
from dnora.utils.io import read_ww3_nml
import os
import shutil
def cleanup():
    if os.path.isdir("TestGrid_WW3"):
        shutil.rmtree("TestGrid_WW3")

@pytest.fixture(scope="session")
def grid():
    grid = dn.grid.Constant(lon=(5, 6), lat=(70, 71), name='TestGrid')
    grid.set_spacing(dlon=0.2, dlat=0.1)
    grid.set_boundary_points(dn.grid.mask.Edges(['W','N','S']))
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
    nml_dict = read_ww3_nml('TestGrid_WW3/ww3_grid.nml')
    assert nml_dict['DEPTH_NML']['DEPTH']['FILENAME'] == "'TestGrid_WW3/TestGrid_bathy.txt'"
    assert nml_dict['MASK_NML']['MASK']['FILENAME'] == "'TestGrid_WW3/TestGrid_mapsta.txt'"
    
    exe.write_grid_file(folder_on_server='/lustre/folder/myfolder')
    nml_dict = read_ww3_nml('TestGrid_WW3/ww3_grid.nml')
    assert nml_dict['DEPTH_NML']['DEPTH']['FILENAME'] == "'/lustre/folder/myfolder/TestGrid_bathy.txt'"
    assert nml_dict['MASK_NML']['MASK']['FILENAME'] == "'/lustre/folder/myfolder/TestGrid_mapsta.txt'"
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
    nml_dict = read_ww3_nml('TestGrid_WW3/ww3_prnc.nml')

    assert nml_dict['FILE_NML']['FILE']['FILENAME'] == "'TestGrid_WW3/wind_ConstantData_TestGrid_20200101T0000_20200131T2300.nc'"
    assert nml_dict['FILE_NML']['FILE']['LONGITUDE'] == "'lon'"
    assert nml_dict['FILE_NML']['FILE']['LATITUDE'] == "'lat'"
    assert nml_dict['FILE_NML']['FILE']['VAR(1)'] == "'u'"
    assert nml_dict['FILE_NML']['FILE']['VAR(2)'] == "'v'"
    
    
    exe.write_wind_file(folder_on_server='/lustre/folder/myfolder')
    nml_dict = read_ww3_nml('TestGrid_WW3/ww3_prnc.nml')
    assert nml_dict['FILE_NML']['FILE']['FILENAME'] == "'/lustre/folder/myfolder/wind_ConstantData_TestGrid_20200101T0000_20200131T2300.nc'"
    assert nml_dict['FILE_NML']['FILE']['LONGITUDE'] == "'lon'"
    assert nml_dict['FILE_NML']['FILE']['LATITUDE'] == "'lat'"
    assert nml_dict['FILE_NML']['FILE']['VAR(1)'] == "'u'"
    assert nml_dict['FILE_NML']['FILE']['VAR(2)'] == "'v'"
    cleanup()