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

