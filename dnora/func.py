import numpy as np

# Import objects
from .grd.grd_mod import Grid

from .grd.read import ForceFeed
from .grd.mesh import TrivialMesher
from .grd.boundary import SetMatrix

# Import aux_funcsiliry functions
from .aux_funcs import read_ww3_info, add_folder_to_filename

def regenerate_ww3(gridname: str, folder: str='') -> Grid:
    """Recreate a WW3 grid object based on the _info, _bathy and _mapsta files"""

    filename = add_folder_to_filename(f'{gridname}_info.txt', folder)

    print(filename)
    lon_min, lon_max, lat_min, lat_max, dlon, dlat, NX, NY = read_ww3_info(filename)

    filename = add_folder_to_filename(f'{gridname}_bathy.txt', folder)
    topo=np.loadtxt(filename).reshape((NY,NX))
    filename = add_folder_to_filename(f'{gridname}_mapsta.txt', folder)
    mask=np.loadtxt(filename).reshape((NY,NX)) == 2 # Boundary points given as value 2

    grid = Grid(lon_min, lon_max, lat_min, lat_max, name=gridname)
    grid.set_spacing(nx = NX, ny = NY)
    grid.import_topo(topo_reader=ForceFeed(topo, grid.lon(), grid.lat()))
    grid.mesh_grid(TrivialMesher())
    grid.set_boundary(boundary_setter=SetMatrix(mask))

    return grid
