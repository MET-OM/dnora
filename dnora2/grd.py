from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from scipy.interpolate import griddata
import dnora2.bnd as bnd
import matplotlib.pyplot as plt
import sys

def msg_print_line(length = 75):
    print("-" * length)

def msg_to_file(filename):
    msg_info(f"Writing to file {filename}")
    
def msg_plain(msg):
    print(msg)

def msg_info(msg):
    print(f"*** {msg} ***")

def msg_advice(msg):
    print(f"!!! {msg} !!!")
    
    
def msg_templates(code):
    if code == 'no_spacing':
        msg_info("No information about grid spacing exists")
        msg_advice("First set spacing with .set_spacing(dlon, dlat) or .set_spacing(dm) (in metres)")
    elif code == 'no_topo':
        msg_info("No topography exists")
        msg_advice("First load topography with .import_topo(topography_fetcher)")
    elif code == 'no_mask':
        msg_info("No land-sea mask exists")
        msg_advice("First load topography with .import_topo(topography_fetcher)")


class Mesher(ABC):
    def __init__(self):
        pass
    
    
class BilinearMesher(Mesher):
    def __init__(self):
        pass
    
    def __call__(self, data, lon, lat, lonQ, latQ):
        lon0, lat0 = np.meshgrid(lon, lat)
        lon1, lat1 = np.meshgrid(lonQ, latQ)
        data[np.isnan(data)] = 0 # Keeping land points as nan lets the shoreline creep out
        M = np.column_stack((data.ravel(), lon0.ravel(),lat0.ravel()))
        meshed_data = griddata(M[:,1:], M[:,0], (lon1, lat1), method='linear')
        meshed_data[meshed_data>=0] = 32767
        
        return meshed_data

    
    
class TopoFetcher(ABC):
    @abstractmethod
    def __init__(self):
        pass
   
   
class TopoEMODNET2018(TopoFetcher):
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


class TopoTrivial(TopoFetcher):
    def __init__(self):
        pass
    
    def __call__(self, lon_min, lon_max, lat_min, lat_max):
        # Creates a trivial topography with all water points
        topo = np.ones((2,2))*-9999
        topo_lon = np.array([lon_min, lon_max])
        topo_lat = np.array([lat_min, lat_max])
        return topo, topo_lon, topo_lat





class Grid(ABC):
    def __init__(self, lon_min, lon_max, lat_min, lat_max, dlon = 0, dlat = 0, dm = 0, name="AnonymousGrid"):
        self.lon_min = lon_min
        self.lat_min = lat_min
        self.lon_max = lon_max
        self.lat_max = lat_max
        self.name = name
        
        if any([dlon, dlat, dm]):
            self.set_spacing(dlon, dlat, dm)
        else:
            msg_info("Grid initialized")
            msg_advice("Set spacing with .set_spacing(dlon, dlat) or .set_spacing(dm) (in metres)")
            print(self)
            
 
 
        return        
    
    def __str__(self):
        msg_print_line()
        msg_plain(f'lon: {self.lon_min}-{self.lon_max}, lat: {self.lat_min}-{self.lat_max}')
        if hasattr(self, 'dlon') and hasattr(self, 'dlon'):
            msg_plain(f'dlon, dlat={self.dlon}, {self.dlat} deg')
            msg_plain(f'dlon, dlat = 1/{1/self.dlon}, 1/{1/self.dlat} deg')
        if hasattr(self, 'dx') and hasattr(self, 'dy'):
            msg_plain(f'dx ,dy approximately {self.dx}, {self.dy} metres')
        if hasattr(self, 'nx') and hasattr(self, 'ny'):
            msg_plain(f'nx, ny = {self.nx} x {self.ny} grid points')
        if hasattr(self, 'topo'):
            msg_plain(f"Mean depth: {np.mean(self.topo[self.mask>0]):.1f} m")
            msg_plain(f"Max depth: {np.min(self.topo[self.mask>0]):.1f} m")
            msg_plain(f"Min depth: {np.max(self.topo[self.mask>0]):.1f} m")
        if hasattr(self, 'mask'):
            msg_print_line()
            msg_plain('Grid contains:')
            msg_plain(f'{sum(sum(self.mask==1)):d} sea points')
            msg_plain(f'{sum(sum(self.mask==0)):d} land points')
            msg_plain(f'{sum(sum(np.logical_and(self.bnd, self.mask==1))):d} boundary points')
        msg_print_line()
        
        return ''

    def import_topo(self, topo_fetcher, mesher = BilinearMesher()):
        if hasattr(self, 'lon') and hasattr(self, 'lon'):
            msg_info('Importing topography')
            topo, topo_lon, topo_lat = topo_fetcher(self.lon_min, self.lon_max, self.lat_min, self.lat_max)
    
            
            msg_info('Generating grid bathymetry')
            self.topo = mesher(topo, topo_lon, topo_lat, self.lon, self.lat)
            
            self.mask = self.set_land_sea_mask()
            
            # This creates an empty one if it didn't exist. Keeps an old boundary if it existed. 
            self.set_boundary()
    
        else:
            msg_templates('no_spacing')

    def set_land_sea_mask(self):
        mask_map = np.copy(self.topo)
        mask_map[mask_map >= 0] = 0 # land points
        mask_map[mask_map < 0] = 1 # sea points
        
        return mask_map

    def set_boundary(self, bounN = -1, edges = ['N', 'S', 'E', 'W']):
        
        
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
            msg_info(f"Setting boundary points using bounN = {bounN}")
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
            msg_info(f'Set {sum(sum(self.bnd)):d} boundary points (some can be on land)')   
        elif bounN == 0:
            if hasattr(self, 'bnd'):
                print("*** Removing all boundary points ***")
            self.bnd=np.full(self.mask.shape, False)
            
        elif bounN < 0:
            msg_info("Kept old boundary. Use bounN = 0 to clear boundary points and bounN > 0 to set them.")
            
        print(self)
        return
    
    def bnd_points(self):
        """Returns the set boundary points"""
        mask_flat = np.logical_and(self.bnd.ravel(), self.mask.ravel()) 
        meshlon, meshlat=np.meshgrid(self.lon,self.lat)
        lonlat_flat = np.column_stack((meshlon.ravel(),meshlat.ravel()))
        BOUND = np.column_stack((lonlat_flat,mask_flat))
        BOUND = BOUND[BOUND[:,2] == 1] # 1 = True
        BOUND = BOUND[:,0:2]
        
        return BOUND


    def set_spacing(self, dlon = 0, dlat = 0, dm = 0, floating_edge = False):
        """Defines longitude and latitude vectors based on desired spacing"""
        if dlon and dlat:
            if floating_edge:
                #Use exactly given dlon/dlat and change lon_max/lat_max accordingly

                self.dlon = dlon
                self.dlat = dlat
                
                
                msg_info(f"Setting spacing based on dlon = {self.dlon} and dlat = {self.dlat}")
                msg_info("floating_edge = True. Making sure dlon/dlat are keep exactly fixed")
                self.lon=np.arange(self.lon_min,self.lon_max+self.dlon/2,self.dlon)  
                self.lat=np.arange(self.lat_min,self.lat_max+self.dlon/2,self.dlat)        
                
                msg_info(f"Setting lon_max ({self.lon_max} >> {self.lon[-1]}), lat_max ({self.lat_max} >> {self.lat[-1]})")
                self.lon_max = self.lon[-1]
                self.lat_max = self.lat[-1]
                
                distance_x = bnd.distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
                distance_y = bnd.distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)
    
                # Number of points
                self.nx = int((self.lon_max-self.lon_min)/self.dlon + 1)
                self.ny = int((self.lat_max-self.lat_min)/self.dlat + 1)
                
                # dx, dy in metres
                self.dx = distance_x*1000/self.nx
                self.dy = distance_y*1000/self.ny
            
            else:
                # Keeping edges fixed and rounding dlon/dlat to something suitable                
            
                msg_info(f"Setting spacing based on dlon = {dlon} and dlat = {dlat}")
                
                # Number of points    
                self.nx = int((self.lon_max-self.lon_min)/dlon + 1)
                self.ny = int((self.lat_max-self.lat_min)/dlat + 1)
                
                # Define longitudes and latitudes
                self.lon = np.linspace(self.lon_min, self.lon_max, self.nx)
                self.lat = np.linspace(self.lat_min, self.lat_max, self.ny)
                self.dlon = (self.lon_max-self.lon_min)/(self.nx-1)
                self.dlat = (self.lat_max-self.lat_min)/(self.ny-1)
                
                distance_x = bnd.distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
                distance_y = bnd.distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)
            
                # dx, dy in metres
                self.dx = distance_x*1000/self.nx
                self.dy = distance_y*1000/self.ny
            if dm:
                msg_info(f"Ignoring value dm={dm} metres!")
        elif dm:
            msg_info(f"Setting spacing based on (approximately) dm={dm} metres")
            distance_x = bnd.distance_2points((self.lat_min+self.lat_max)/2, self.lon_min, (self.lat_min+self.lat_max)/2, self.lon_max)
            distance_y = bnd.distance_2points(self.lat_min, self.lon_min, self.lat_max, self.lon_min)

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

        else:
            msg_advice("Doing nothing. Run set_spacing with either dlon AND dlat (in degrees), or dm (in metres).")

        
        # Initialize the grid with an empty topography
        msg_info("Initializing with an empty topography")
        topo_fetcher = TopoTrivial()
        self.import_topo(topo_fetcher)




    def plot_topo(self, save_fig = False, filename = ''):
        if hasattr(self, 'topo'):
            if not filename:
                filename = f"{self.name}_topo.pdf"
            levels = np.linspace(np.min(self.topo), 0, 100, endpoint=True)
            plt.figure()
            plt.contourf(self.lon,self.lat,self.topo,levels)
            plt.colorbar()
            plt.title(f"{self.name} topograpy")
            #plt.show()
            if save_fig:
                plt.savefig(filename, dpi=300)
                msg_to_file(filename)
            #plt.clf()
        else:
            msg_templates('no_topo')
            
    def plot_mask(self, save_fig = False, filename = ''):
        
        if hasattr(self, 'topo'):

            if not filename:
                filename = f"{self.name}_mask.pdf"
                
            plt.figure()
            plt.contourf(self.lon,self.lat,self.mask)
            BOUND = self.bnd_points()
            plt.plot(BOUND[:,0], BOUND[:,1],'r*')
            plt.colorbar()
            plt.title(f"{self.name} land-sea mask")
            if save_fig:
                plt.savefig(filename, dpi=300)
                msg_to_file(filename)
            #plt.clf()
        
        else:
            msg_templates('no_mask')


    def write_status(self, filename = ''):
        if not filename:
            filename = f"{self.name}_info.txt"
        
        msg_to_file(filename)
        
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
            msg_info(f"Setting points shallower than {min_depth} to land")
            self.topo[np.logical_and(shallow_points, self.mask)] = land_value # Don't touch land points by usign self.mask
            msg_info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, self.mask))} points")
            
            # Redfine the land-sea mask
            msg_info("Updating land-sea mask")
            self.mask = self.set_land_sea_mask()
        else:
            msg_info(f"Setting points shallower than {min_depth} to {min_depth}")
            # Set points to the limiter
            self.topo[np.logical_and(shallow_points, self.mask)] = min_depth # Don't touch land points by usign self.mask
            msg_info(f"Affected {np.count_nonzero(np.logical_and(shallow_points, self.mask))} points")
        
        

        
class WW3Grid(Grid):

    
    def write_topo(self, matrix = False):
        msg_info('Create files for regular grid')
        bnd_mask = np.logical_and(self.bnd, self.mask==1)
        mask_out = self.mask
        if np.logical_and(self.bnd, self.mask==1).any():
            msg_info(f'Masking out {sum(sum(bnd_mask)):d} boundary points in grid...')   
            mask_out[bnd_mask] = 2
        
        if matrix:
            fn1 = 'mat_'+self.name+'_bathy.txt'
            msg_to_file(fn1)
            np.savetxt(fn1, -1*self.topo, delimiter=',',fmt='%1.6f')
            
            fn2 = 'mat_'+self.name+'_mapsta.txt'
            msg_to_file(fn2)
            np.savetxt(fn2, mask_out, delimiter=',',fmt='%1.0f')
        else:
            fn1 = self.name+'_bathy.txt'
            msg_to_file(fn1)
            np.savetxt(fn1, -1*self.topo.ravel(), delimiter=',',fmt='%1.6f')
            
            fn2 = self.name+'_mapsta.txt'
            msg_to_file(fn2)
            np.savetxt(fn2, mask_out.ravel(), delimiter=',',fmt='%1.0f')
    
    def __str__(self):
        msg_print_line()
        msg_plain(f'WW3 Grid: {self.name}')
        return super().__str__()



    