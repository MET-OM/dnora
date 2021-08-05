from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from scipy.interpolate import griddata
from dnora2.bnd import distance_2points
import dnora2.msg as msg
import matplotlib.pyplot as plt
import sys

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
    grid = WW3Grid(lon_min, lon_max, lat_min, lat_max, name = gridname)
    grid.set_spacing(nx = NX, ny = NY) 
    topo_fetcher = TopoForceFeed(topo, grid.lon, grid.lat)
    grid.import_topo(topo_fetcher, TrivialMesher())
    grid.set_boundary(given_bnd = mask)
    
    return grid
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# FILTER CLASSES RESPONSIBLE FOR FILTERING THE BATHYMETRIC INFORMATION
# -----------------------------------------------------------------------------
class Filter(ABC):
    def __init__(self):
        
        pass
      

class TrivialFilter(Filter):
    def __init__(self):
        pass
    
    def __call__(self, data, lon, lat):
        return data
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# MESHER CLASSES RESPONSIBLE FOR INTERPOLATING THE BATHYMETRIC INFORMATION
# -----------------------------------------------------------------------------
class Mesher(ABC):
    @abstractmethod
    def __init__(self):
        pass
    
    @abstractmethod
    def __call__(self, data, lon, lat, lonQ, latQ):
        pass
        
class InterpolationMesher(Mesher):
    def __init__(self, method = 'linear'):
        self.method = method
        #msg.info(f"Initializing mesher with method: {self.method}")

        return
    
    def __call__(self, data, lon, lat, lonQ, latQ):
        lon0, lat0 = np.meshgrid(lon, lat)
        lon1, lat1 = np.meshgrid(lonQ, latQ)
        data[np.isnan(data)] = 0 # Keeping land points as nan lets the shoreline creep out
        M = np.column_stack((data.ravel(), lon0.ravel(),lat0.ravel()))
        meshed_data = griddata(M[:,1:], M[:,0], (lon1, lat1), method=self.method)
        meshed_data[meshed_data>=0] = 32767
        
        return meshed_data

class TrivialMesher(Mesher):
    def __init__(self):
        pass
    
    def __call__(self, data, lon, lat, lonQ, latQ):
        return data
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# TOPOFETCH CLASSES RESPONSIBLE IMPORTING THE BATHYMETRIC INFORMATION
# -----------------------------------------------------------------------------
class TopoFetcher(ABC):
    @abstractmethod
    def __init__(self):
        pass
    
    @abstractmethod
    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        pass
 
    
class TopoEMODNET2018(TopoFetcher):
    """Reads data from EMODNET"""
    def __init__(self):
        self.source='/lustre/storeB/project/fou/om/WW3/bathy/emodnet_115m_x_115m/C5_2018.dtm' 
        print("*** Initialized Topography reader for EMODNET C5_2018 ***")
        return
    
    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # If we limit ourselves to exactly the grid, we will get nans at the edges in the interpolation. Add 10% tolerance around all edges.
        tolerance_lon = (lon_max-lon_min)*0.1
        tolerance_lat = (lat_max-lat_min)*0.1
        
        ds = xr.open_dataset(self.source).sel(COLUMNS=slice(lon_min-tolerance_lon, lon_max+tolerance_lon), LINES=slice(lat_min-tolerance_lat, lat_max+tolerance_lat))
        topo = ds.DEPTH.values
        topo_lon = ds.COLUMNS.values
        topo_lat = ds.LINES.values
        return topo, topo_lon, topo_lat


class TopoEmpty(TopoFetcher):
    """Creates an empty topography. Called when setting initial spacing."""
    def __init__(self):
        pass
    
    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # Creates a trivial topography with all water points
        topo = np.ones((2,2))*-9999
        topo_lon = np.array([lon_min, lon_max])
        topo_lat = np.array([lat_min, lat_max])
        return topo, topo_lon, topo_lat


class TopoForceFeed(TopoFetcher):
    """Simply passes on the data it was fed upon initialization"""
    def __init__(self, topo, topo_lon, topo_lat):
        self.topo = topo
        self.topo_lon = topo_lon
        self.topo_lat = topo_lat
        return

    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # Just use the values it was forcefed on initialization
        topo = self.topo.copy()
        topo_lon = self.topo_lon.copy()
        topo_lat = self.topo_lat.copy()
        return topo, topo_lon, topo_lat
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# GRID CLASSES RESPONSIBLE ALL THE GRID INFORMATION AND METHODS
# -----------------------------------------------------------------------------
class Grid(ABC):
    def __init__(self, lon_min, lon_max, lat_min, lat_max, name="AnonymousGrid"):
        self.lon_min = lon_min
        self.lat_min = lat_min
        self.lon_max = lon_max
        self.lat_max = lat_max
        self.name = name
        
        return        
    
    def import_topo(self, topo_fetcher, mesher = InterpolationMesher(method = 'linear'), filt = TrivialFilter()):
        """Main function responsible for initialising a proper grid. Reads to topography, which is filtered and meshed to the grid spacing. Finally the land-sea mask and an empty boundary mask are initialized."""
        if hasattr(self, 'lon') and hasattr(self, 'lon'):
            msg.info('Importing topography')
            topo, topo_lon, topo_lat = topo_fetcher(self.lon_min, self.lon_max, self.lat_min, self.lat_max)
    
            
            msg.info('Filtering bathmetry')
            topo_filt = filt(topo, topo_lon, topo_lat)
    
            msg.info('Generating grid bathymetry')
            self.topo = mesher(topo_filt, topo_lon, topo_lat, self.lon, self.lat)
            
            msg.info('Setting land-sea mask')
            self.mask = self.land_sea_mask()
            
            # This creates an empty one if it didn't exist. Keeps an old boundary if it existed. 
            self.set_boundary()
    
        else:
            msg.templates('no_spacing')

    def land_sea_mask(self):
        """Returns an array with the land-sea mask (1 = sea, 0 = land)."""
        mask_map = np.copy(self.topo)
        mask_map[mask_map >= 0] = 0 # land points
        mask_map[mask_map < 0] = 1 # sea points
        
        return mask_map

    def set_boundary(self, bounN = -1, edges = ['N', 'S', 'E', 'W'], given_array = None):
        """Define boundary points in grid either by setting every bounN point at edges, or providins a boolean mask (array)."""
        if given_array is not None: # Array where boundary points are True and rest False is given
            self.bnd = given_array
        else:
            # Initializing condition
            if not hasattr(self, 'bnd') and bounN <= 0:
                bounN = 0
                
            bounN = int(bounN)
    
            #self.bnd=np.full(self.mask.shape, False)
            #if bounN > 0:
            #    self.bounN = bounN # Update if given by user
            
            # In case of no input we still use a possible old saved value
            #if hasattr(self, 'bounN'):
            #    bounN = self.bounN
                   
            if bounN >0: 
                msg.info(f"Setting boundary points using bounN = {bounN}")
                # --------- North boundary ----------
                if 'N' in edges:
                    self.bnd[-1,::bounN] = True 
                ## --------- South boundary ----------
                if 'S' in edges:
                    self.bnd[0,::bounN] = True
                ## --------- East boundary ----------
                if 'E' in edges:
                    self.bnd[::bounN,-1] = True
                ## --------- West boundary ----------
                if 'W' in edges:
                    self.bnd[::bounN,0] = True
                msg.info(f'Set {sum(sum(self.bnd)):d} boundary points (some can be on land)')   
            elif bounN == 0:
                if hasattr(self, 'bnd'):
                    print("*** Removing all boundary points ***")
                self.bnd=np.full(self.mask.shape, False)
                
            elif bounN < 0:
                msg.info("Kept old boundary. Use bounN = 0 to clear boundary points and bounN > 0 to set them.")
           
        print(self)
        return
    
    def bnd_points(self):
        """Returns a lon, lat list of the set boundary points."""
        mask_flat = np.logical_and(self.bnd.ravel(), self.mask.ravel()) 
        meshlon, meshlat=np.meshgrid(self.lon,self.lat)
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        BOUND = np.column_stack((lonlat_flat,mask_flat))
        BOUND = BOUND[BOUND[:,2] == 1] # 1 = True
        BOUND = BOUND[:,0:2]
        
        return BOUND


    def set_spacing(self, dlon = 0, dlat = 0, dm = 0, nx = 0, ny = 0, floating_edge = False):
        """Defines longitude and latitude vectors based on desired spacing. Use either dlon and dlat (in deg), dm (in m), or nx and ny (in grid points)."""
        if dlon and dlat:
            if floating_edge:
                #Use exactly given dlon/dlat and change lon_max/lat_max accordingly

                self.dlon = dlon
                self.dlat = dlat
                
                
                msg.info(f"Setting spacing based on dlon = {self.dlon} and dlat = {self.dlat}")
                msg.info("floating_edge = True. Making sure dlon/dlat are keep exactly fixed")
                self.lon=np.arange(self.lon_min,self.lon_max+self.dlon/2,self.dlon)  
                self.lat=np.arange(self.lat_min,self.lat_max+self.dlon/2,self.dlat)        
                
                msg.info(f"Setting lon_max ({self.lon_max} >> {self.lon[-1]}), lat_max ({self.lat_max} >> {self.lat[-1]})")
                self.lon_max = self.lon[-1]
                self.lat_max = self.lat[-1]
                
                distance_x = distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
                distance_y = distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)
    
                # Number of points
                self.nx = int((self.lon_max-self.lon_min)/self.dlon + 1)
                self.ny = int((self.lat_max-self.lat_min)/self.dlat + 1)
                
                # dx, dy in metres
                self.dx = distance_x*1000/self.nx
                self.dy = distance_y*1000/self.ny
            
            else:
                # Keeping edges fixed and rounding dlon/dlat to something suitable                
            
                msg.info(f"Setting spacing based on dlon = {dlon} and dlat = {dlat}")
                
                # Number of points    
                self.nx = int((self.lon_max-self.lon_min)/dlon + 1)
                self.ny = int((self.lat_max-self.lat_min)/dlat + 1)
                
                # Define longitudes and latitudes
                self.lon = np.linspace(self.lon_min, self.lon_max, self.nx)
                self.lat = np.linspace(self.lat_min, self.lat_max, self.ny)
                self.dlon = (self.lon_max-self.lon_min)/(self.nx-1)
                self.dlat = (self.lat_max-self.lat_min)/(self.ny-1)
                
                distance_x = distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
                distance_y = distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)
            
                # dx, dy in metres
                self.dx = distance_x*1000/self.nx
                self.dy = distance_y*1000/self.ny
            if dm:
                msg.info(f"Ignoring value dm={dm} metres!")
        elif dm:
            msg.info(f"Setting spacing based on (approximately) dm={dm} metres")
            distance_x = distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
            distance_y = distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)

            # Number of points
            self.nx = int(np.round(distance_x*1000/dm))
            self.ny = int(np.round(distance_y*1000/dm))
            
            # dx, dy in metres
            self.dx = distance_x*1000/self.nx
            self.dy = distance_y*1000/self.ny


            # Define longitudes and latitudes
            self.lon = np.linspace(self.lon_min, self.lon_max, self.nx)
            self.lat = np.linspace(self.lat_min, self.lat_max, self.ny)
            self.dlon = (self.lon_max-self.lon_min)/(self.nx-1)
            self.dlat = (self.lat_max-self.lat_min)/(self.ny-1)

        elif nx and ny:
            msg.info(f"Setting spacing to have nx = {nx}, ny = {ny} points.")
            self.nx = nx
            self.ny = ny

            # Define longitudes and latitudes
            self.lon = np.linspace(self.lon_min, self.lon_max, self.nx)
            self.lat = np.linspace(self.lat_min, self.lat_max, self.ny)
            self.dlon = (self.lon_max-self.lon_min)/(self.nx-1)
            self.dlat = (self.lat_max-self.lat_min)/(self.ny-1)
            
            distance_x = distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
            distance_y = distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)
        
            # dx, dy in metres
            self.dx = distance_x*1000/self.nx
            self.dy = distance_y*1000/self.ny

        else:
            msg.advice("Doing nothing. Run set_spacing with either dlon AND dlat (in degrees), nx AND ny (in grid points), or dm (in metres).")

        
        # Initialize the grid with an empty topography
        msg.info("Initializing with an empty topography")
        topo_fetcher = TopoEmpty()
        self.import_topo(topo_fetcher)

    def plot_topo(self, save_fig = False, filename = ''):
        """Creates a plot of the topography"""
        
        if hasattr(self, 'topo'):
            if not filename:
                filename = f"{self.name}_topo.pdf"
            levels = np.linspace(np.min(self.topo), 0, 100, endpoint=True)
            
            plt.figure()
            plt.contourf(self.lon,self.lat,self.topo,levels)
            
            plt.colorbar()
            plt.title(f"{self.name} topograpy")
            plt.show()
            
            if save_fig:
                plt.savefig(filename, dpi=300)
                msg.to_file(filename)
        else:
            msg.templates('no_topo')
  
            
    def plot_mask(self, save_fig = False, filename = ''):
        """Creates a plot of the land-sea mask"""
        if hasattr(self, 'topo'):
            if not filename:
                filename = f"{self.name}_mask.pdf"
                
            plt.figure()
            plt.contourf(self.lon,self.lat,self.mask)
            BOUND = self.bnd_points()
            plt.plot(BOUND[:,0], BOUND[:,1],'r*')
            
            plt.colorbar()
            plt.title(f"{self.name} land-sea mask")
            plt.show()
            
            if save_fig:
                plt.savefig(filename, dpi=300)
                msg.to_file(filename)
        else:
            msg.templates('no_mask')

    def write_status(self, filename = ''):
        """Writes out the status of the grid to a file."""
        if not filename:
            filename = f"{self.name}_info.txt"
        
        msg.to_file(filename)
        
        stdout = sys.stdout
        sys.stdout = open(filename, 'w')
        print(self)
        sys.stdout.close()  
        sys.stdout = stdout


    def set_min_depth(self, min_depth = 0, to_land = False):
        """Sets a minimun depth in the grid. Use to_land = True to set shallower point to land points (isntead of setting to the minimum depth)."""
        
        if min_depth > 0:
            min_depth = -min_depth
        
        shallow_points = self.topo > min_depth
        
        if to_land:
            land_value = max([self.topo.max(), 0])
            msg.info(f"Setting points shallower than {min_depth} to land")
            self.topo[np.logical_and(shallow_points, self.mask)] = land_value # Don't touch land points by usign self.mask
            msg.info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, self.mask))} points")
            
            # Redfine the land-sea mask
            msg.info("Updating land-sea mask")
            self.mask = self.land_sea_mask()
        else:
            msg.info(f"Setting points shallower than {min_depth} to {min_depth}")
            # Set points to the limiter
            self.topo[np.logical_and(shallow_points, self.mask)] = min_depth # Don't touch land points by usign self.mask
            msg.info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, self.mask))} points")
        
    def __str__(self):
            msg.print_line()
            msg.plain(f'lon: {self.lon_min} - {self.lon_max}, lat: {self.lat_min} - {self.lat_max}')
            if hasattr(self, 'dlon') and hasattr(self, 'dlon'):
                msg.plain(f'dlon, dlat = {self.dlon}, {self.dlat} deg')
                msg.plain(f'Inverse: dlon, dlat = 1/{1/self.dlon}, 1/{1/self.dlat} deg')
            if hasattr(self, 'dx') and hasattr(self, 'dy'):
                msg.plain(f'dx, dy approximately {self.dx}, {self.dy} metres')
            if hasattr(self, 'nx') and hasattr(self, 'ny'):
                msg.plain(f'nx, ny = {self.nx} x {self.ny} grid points')
            if hasattr(self, 'topo'):
                msg.plain(f"Mean depth: {np.mean(self.topo[self.mask>0]):.1f} m")
                msg.plain(f"Max depth: {np.min(self.topo[self.mask>0]):.1f} m")
                msg.plain(f"Min depth: {np.max(self.topo[self.mask>0]):.1f} m")
            if hasattr(self, 'mask'):
                msg.print_line()
                msg.plain('Grid contains:')
                msg.plain(f'{sum(sum(self.mask==1)):d} sea points')
                msg.plain(f'{sum(sum(self.mask==0)):d} land points')
                msg.plain(f'{sum(sum(np.logical_and(self.bnd, self.mask==1))):d} boundary points')
            msg.print_line()
            
            return ''      

        
class WW3Grid(Grid):

    
    def write_topo(self, matrix = False):
        msg.info('Create files for regular grid')
        bnd_mask = np.logical_and(self.bnd, self.mask==1)
        mask_out = self.mask.copy()
        if np.logical_and(self.bnd, self.mask==1).any():
            msg.info(f'Setting {sum(sum(bnd_mask)):d} boundary points in grid...')   
            mask_out[bnd_mask] = 2
        
        if matrix:
            fn1 = 'mat_'+self.name+'_bathy.txt'
            msg.to_file(fn1)
            np.savetxt(fn1, -1*self.topo, delimiter=',',fmt='%1.6f')
            
            fn2 = 'mat_'+self.name+'_mapsta.txt'
            msg.to_file(fn2)
            np.savetxt(fn2, mask_out, delimiter=',',fmt='%1.0f')
        else:
            fn1 = self.name+'_bathy.txt'
            msg.to_file(fn1)
            np.savetxt(fn1, -1*self.topo.ravel(), delimiter=',',fmt='%1.6f')
            
            fn2 = self.name+'_mapsta.txt'
            msg.to_file(fn2)
            np.savetxt(fn2, mask_out.ravel(), delimiter=',',fmt='%1.0f')
            
        self.write_status()
    
    def __str__(self):
        msg.print_line()
        msg.plain(f'WW3 Grid: {self.name}')
        return super().__str__()

    