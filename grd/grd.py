import numpy as np
from scipy.interpolate import griddata
from .read import ForceFeed
from ..obj import Grid

# -----------------------------------------------------------------------------
# MISC STAND ALONE FUNCTIONS
# -----------------------------------------------------------------------------
def read_ww3_info(filename):
    """Read grid specification from the GridName_info.txt file"""
    with open(filename,'r') as f:
        lines = f.readlines()

    for n in range (len(lines)):
        line = lines[n].split()

        if len(line):
            if line[0] == 'lon:':
                lon_min = float(line[1])
                lon_max = float(line[3][0:-1])
                lat_min = float(line[5])
                lat_max = float(line[7])
            elif line[0] == 'dlon,':
                dlon = float(line[3][0:-1])
                dlat = float(line[4])
            elif line[0] == 'nx,':
                nx = int(line[3])
                ny = int(line[5])
    return lon_min, lon_max, lat_min, lat_max, dlon, dlat, nx, ny

def regenerate_ww3(gridname):
    """Recreate a WW3 grid object based on the _info, _bathy and _mapsta files"""
    lon_min, lon_max, lat_min, lat_max, dlon, dlat, NX, NY = read_ww3_info(f'{gridname}_info.txt')

    topo=-np.loadtxt(f'{gridname}_bathy.txt').reshape((NY,NX))
    mask=np.loadtxt(f'{gridname}_mapsta.txt').reshape((NY,NX)) == 2 # Boundary points given as value 2

    # Regenerate grid by force feeding data to the TopoFetcher
    grid = Grid(lon_min, lon_max, lat_min, lat_max, name = gridname)
    grid.set_spacing(nx = NX, ny = NY)
    topo_fetcher = ForceFeed(topo, grid.lon(), grid.lat())
    grid.import_topo(topo_fetcher)
    grid.mesh_grid(TrivialMesher())
    grid.set_boundary(given_array = mask)

    return grid
# -----------------------------------------------------------------------------








