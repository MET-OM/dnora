from abc import ABC, abstractmethod
import xarray as xr
import numpy as np
from scipy.interpolate import griddata


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

class Grid(ABC):
    @abstractmethod
    def __init__(self):
        pass
    
    def __str__(self):
        print(f'lon: {self.lon_min}-{self.lon_max}, lat{self.lat_min}-{self.lat_max}, dlon/dlat={self.dlon}/{self.dlat}')
        if hasattr(self, 'mask'):
            print('Grid contains:')
            print(f'{sum(sum(self.mask==1)):d} sea points')
            print(f'{sum(sum(self.mask==0)):d} land points')
            print(f'{sum(sum(np.logical_and(self.bnd, self.mask==1))):d} boundary points')
        return ''


    def set_land_sea_mask(self):
        mask_map = np.copy(self.topo)
        mask_map[mask_map >= 0] = 0 # land points
        mask_map[mask_map < 0] = 1 # sea points
        
        return mask_map

    def set_boundary(self, bounN = 0):
        print('Setting boundary points...')
        bounN = int(bounN)
        
        self.bnd=np.full(self.mask.shape, False)
        
        if bounN > 0:
            self.bounN = bounN # Update if given by user
        
        # In case of no input we still use a possible old saved value
        if hasattr(self, 'bounN'):
            bounN = self.bounN
               
        if bounN: 
            # --------- North boundary ----------
            self.bnd[-1,::bounN] = True 
            ## --------- East boundary ----------
            self.bnd[::bounN,-1] = True
            ## --------- West boundary ----------
            self.bnd[::bounN,0] = True
            ## --------- South boundary ----------
            self.bnd[0,::bounN] = True
            print(f'Set {sum(sum(self.bnd)):d} boundary points (some can be on land).')   
      
        return
    
    



class WW3Grid(Grid):
    def __init__(self, lon_min, lon_max, lat_min, lat_max, dlon, dlat, name="AnonymousGrid"):
        self.lon_min = lon_min
        self.lat_min = lat_min
        self.lon_max = lon_max
        self.lat_max = lat_max
        self.dlon = dlon
        self.dlat = dlat
        self.name = name
    
        return        


    def import_topo(self, topo_fetcher, mesher = BilinearMesher()):
        topo, topo_lon, topo_lat = topo_fetcher(self.lon_min, self.lon_max, self.lat_min, self.lat_max)
        self.lon=np.arange(self.lon_min,self.lon_max+self.dlon/2,self.dlon)  
        self.lat=np.arange(self.lat_min,self.lat_max+self.dlon/2,self.dlat)  
        
        print('Generating Grid Bathymetry...')
        self.topo = mesher(topo, topo_lon, topo_lat, self.lon, self.lat)
        
        self.mask = self.set_land_sea_mask()
        
        # This resets the boundary if it exists and creates an empty one if it didn't exist 
        self.set_boundary()

        print(self)
    
    def write_topo(self, matrix = False):
        print('Create files for regular grid')
        bnd_mask = np.logical_and(self.bnd, self.mask==1)
        mask_out = self.mask
        if np.logical_and(self.bnd, self.mask==1).any():
            print(f'Masking out {sum(sum(bnd_mask)):d} boundary points in grid...')   
            mask_out[bnd_mask] = 2
        
        if matrix:
            np.savetxt('mat_'+self.name+'_bathy.txt', -1*self.topo, delimiter=',',fmt='%1.6f')
            np.savetxt('mat_'+self.name+'_mapsta.txt', mask_out, delimiter=',',fmt='%1.0f')
        else:
            np.savetxt(self.name+'_bathy.txt', -1*self.topo.ravel(), delimiter=',',fmt='%1.6f')
            np.savetxt(self.name+'_mapsta.txt', mask_out.ravel(), delimiter=',',fmt='%1.0f')
    
    def __str__(self):
        print(f'WW3 Grid {self.name}')
        return super().__str__()



    