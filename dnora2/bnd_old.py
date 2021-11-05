import numpy as np
import pandas as pd
import xarray as xr
from abc import ABC, abstractmethod
import netCDF4
import dnora2.msg as msg
from scipy import interpolate
from statistics import mode

# -----------------------------------------------------------------------------
# MISC STAND ALONE FUNCTIONS
# -----------------------------------------------------------------------------
def lon_lat_from_dataset(bnd_in):
    """Get longitude and latitude from xarray dataset"""
    lat = bnd_in["latitude"][0].values
    lon = bnd_in["longitude"][0].values 
    return lon, lat

def distance_2points(lat1,lon1,lat2,lon2):
    """Calculate distance between two points"""
    R = 6371.0
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    distance = R * c # in km
    return distance


def min_distance(lon, lat, lon_vec, lat_vec):
    dx = []
    for n in range(len(lat_vec)):
        dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))
        
        return np.array(dx).min(), np.array(dx).argmin()

def slice_and_gather_xr(inds, bnd_in):
    """Goes from a list of dataset for each day, to a list of datasets at each point (given by inds)"""
    bnd_out = []
    bnd_mask = [True] * len(inds)
    msg.info("Merging times")
    for n in range(len(inds)): # Loop over output points
        datasets = []
        for k in range(len(bnd_in)): # Loop over time (days)
            new_data=bnd_in[k].sel(x=(inds[n])) # X-indexing goes from 
            # This trivial dimension of one will cause probems since the latitudes are either then defined with only a y-dimension, or in two dimensions (time, y)
            # In the latter case new_data["latitude"].values[0] is not a number, but a one element array of a one element list
            # We solve this by squeezing out this unceccesary dimension
            datasets.append(new_data.squeeze('y')) 
         
        # This list contains one element for each requested output point
        # Each element is a one-spatial-point xr-dataset containing all times
        print(f"Point {n}/{len(inds)}")
        bnd_out.append(xr.concat(datasets, dim="time"))        
        
        if np.isnan(bnd_out[n].SPEC.values).any():
            msg.info(f"Point {n} contains NaN's. Masking as False.")
            bnd_mask[n] = False
    
  
    return bnd_out, bnd_mask


def flip_spec(spec,D):
    # This check enables us to flip directions with flip_spec(D,D)
    
    if len(spec.shape) == 1:
        flipping_dir = True
        spec = np.array([spec])
    else:
        flipping_dir = False
    spec_flip = np.zeros(spec.shape)

    ind = np.arange(0,len(D), dtype='int')
    dD = np.diff(D).mean()
    steps = D/dD # How many delta-D from 0
    
    ind_flip = ((ind - 2*steps).astype(int) + len(D)) % len(D)
    
    spec_flip=spec[:, list(ind_flip)]
    
    if flipping_dir:
        spec_flip = spec_flip[0]
    return spec_flip


def shift_spec(spec, D, shift = 0):
    # This check enables us to flip directions with flip_spec(D,D)
    if len(spec.shape) == 1:
        shifting_dir = True
        spec = np.array([spec])
    else:
        shifting_dir = False
    spec_shift = np.zeros(spec.shape)

    D = np.round(D)
    ind = np.arange(0,len(D), dtype='int')
    dD = mode(abs(np.diff(D)))
    
    if not (shift/dD).is_integer():
        print('aa')
        raise Exception ('Shift needs to be multiple of frequency resolution! Otherwise interpolation would be needed.')
      
    ind_flip = ((ind + int(shift/dD)).astype(int) + len(D)) % len(D)
    
    spec_shift=spec[:, list(ind_flip)]
    if shifting_dir:
        spec_shift = spec_shift[0]
    return spec_shift
    

def ocean_to_naut(oceanspec, D):
    """Convert spectrum in nautical convention (0 north, 90 east, direction from) to oceanic convention (0 north, 90 east, direction to)"""
    nautspec = shift_spec(oceanspec,D, 180)

    return nautspec


def naut_to_ocean(nautspec, D): # Just defined separately to not make for confusing code
    """Convert spectrum in oceanic convention (0 north, 90 east, direction to) to nautical convention (0 north, 90 east, direction from)"""    
    return ocean_to_naut(nautspec, D)


def ocean_to_math(oceanspec, D):
    """Convert spectrum in oceanic convention (0 north, 90 east, direction to) to mathematical convention (90 north, 0 east, direction to)"""
    
    # Flip direction
    spec_flip = flip_spec(oceanspec, D)
    D_flip = flip_spec(D,D)            

    # Shift 0 to be at 90    
    mathspec = shift_spec(spec_flip, D_flip, -90)

    return mathspec


def ocean_to_math_old(oceanspec, D):
    """Convert spectrum in oceanic convention (0 north, 90 east, direction to) to mathematical convention (90 north, 0 east, direction to)"""
    # This has not yeat been fully tested!

    mathspec=np.zeros(oceanspec.shape)
    
    D_math = ((90+360)-D) % 360
    
    ind_ocean = np.arange(0,len(D), dtype='int')
    dD = np.diff(D).mean()
    step_ocean = D/dD # How many delta-D from 0
    
    #           Original ind + 90 deg shift + flipping direction  + don't want negatives + modulo 360/dD
    ind_math = ((ind_ocean + int(90/dD) - 2*step_ocean).astype(int) + len(D)) % len(D)
        
    mathspec=oceanspec[:, list(ind_math)]
    #mathspec=oceanspec
        
    return mathspec, D_math


def naut_to_math_orig(nautspec, D):
    
    #f = interpolate.RectBivariateSpline(freq_obs,dir_obs, SPEC_obs[j,:,:], kx=1, ky=1, s=0)
    #SPEC_ww3[j,0,:,:] = f(freq_model,dir_model)/Delta_dir_model # Divide by direction to keep ww3-units 
    ######## Change direction axis from 0..355(obs.) to 90...0...95(WW3-format)
    # Step 1a and 1b: convert the meteor. convention of observed spectra (0..180..355) to ocean. convetion (180..360..175)
    
    dir_points = nautspec.shape[1] # Spectrum is freq x dirs
    SPEC_ww3 = np.zeros(nautspec.shape)
    SPEC_ww3_new  = np.zeros(nautspec.shape)
    D=np.round(D)
    SPEC_ww3 = nautspec
    #SPEC_ww3[:,0:dir_points//2] = nautspec[:,dir_points//2:] # Step 1a: 180..355 to start of array
    #SPEC_ww3[:,dir_points//2:] = nautspec[:,0:dir_points//2] # Step 1b: 0..175 to end of array
    # Step 2a,b,c: convert to mathematical convention (90...0...95)
    ind_95 = np.where(D == 90)[0][0]+1 # index for point 95 degrees(for 72 dir): 19
    ind_265 = np.where(D == 270)[0][0]-1 # index for point 265 degrees(for 72 dir): 53
    SPEC_ww3_new[:,ind_265:dir_points] = SPEC_ww3[:,0:ind_95] # Step 2a: transfer dir:(0..90) to the end of array 
    SPEC_ww3_new[:,0:ind_265] = SPEC_ww3[:,ind_95:dir_points] # Step 2b: tranfer dir:(95..355) to the start of array
    # After the Step 2a and 2b, we have dir:(95...355,0..90) so we need to change order to (90...0..95) 
    SPEC_ww3_new[:,:] = SPEC_ww3_new[:,::-1] # Step 2c:change(reverse dir. axis) order from (95...0...90) to (90...0...95)
    #SPEC_ww3_new[:,0] = SPEC_ww3[:,0]
    #SPEC_ww3_new[:,1:] = SPEC_ww3[:,:0:-1] # Step 2c:change(reverse dir. axis) order from (95...0...90) to (90...0...95)

    D_new = (90-D) % 360
    #D_new = D
    return SPEC_ww3_new, D_new


def interp_spec(f, D, S, fi, Di):
    Sleft = S
    Sright = S
    Dleft = -D[::-1] 
    Dright = D + 360
    
    bigS = np.concatenate((Sleft, S, Sright),axis=1)
    bigD = np.concatenate((Dleft, D, Dright))
        
    Finterpolator = interpolate.RectBivariateSpline(f, bigD, bigS, kx=1, ky=1, s=0)
    
    Si = Finterpolator(fi,Di)
    
    return Si
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# SPECTRAL PROCESSOR CLASSES FOR PROCESSING SPECTRA BEFORE OUTPUT
# -----------------------------------------------------------------------------
class SpectralProcessor(ABC):
    def __init__(self):
        pass
    
    @abstractmethod
    def __call__(self, bnd_in, bnd_mask):
        pass


class TrivialSpectralProcessor(SpectralProcessor):
    def __init__(self, calib_spec = 1):
        self.calib_spec = calib_spec
        return
    
    def __call__(self, bnd_in, bnd_mask):
       
        for n in range(len(bnd_in)):
            bnd_in[n].SPEC[:,:,:]=bnd_in[n].SPEC[:,:,:]*self.calib_spec
        return bnd_in, bnd_mask
    

class InterpSpectralProcessor(SpectralProcessor):
    def __init__(self, nbins = 36, start = 0, data_set = None):
        if data_set is not None:
            # Read the values from the example dataset provided
            nbins = len(data_set.direction.values)
        
        dD=int(360/nbins)
        if start > dD:
            msg.info("First bin {start} is defined as larger than the directional resolution {dD}. This might spell trouble!")
        
        self.Di = np.array(range(0,360,dD), dtype='float32') + start
      
        return
    
    def __call__(self, bnd_in, bnd_mask):
        for n in range(len(bnd_in)):
            for k in range(len(bnd_in[n].time)):
                S = bnd_in[n].SPEC[k,:,:].values
                f = bnd_in[n].freq.values
                D = bnd_in[n].direction.values
                bnd_in[n].SPEC[k,:,:] = interp_spec(f, D, S, f, self.Di)
            bnd_in[n] = bnd_in[n].assign_coords(direction=self.Di)
        return bnd_in, bnd_mask        
    
class NaNCleanerSpectralProcessor(SpectralProcessor):
    def __init__(self):
        pass
    
    def __call__(self, bnd_in, bnd_mask):
        bnd_out = []
        bnd_mask_out =[]
        
        for n in range(len(bnd_in)):
            if bnd_mask[n]:
                bnd_out.append(bnd_in[n])
                bnd_mask_out.append(True)
                
        return bnd_out, bnd_mask_out
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# POINT PICKER CLASSES FOR DIFFERENT METHODS TO CHOOSE SPECTRA FROM DATABASE
# -----------------------------------------------------------------------------
class PointPicker(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, grid, bnd_in):
        return bnd_in


class PPTrivialPicker(PointPicker):
    def __init__(self):
        pass

    def __call__(self, grid, bnd_in):
        return bnd_in


class PPNearestGridPoint(PointPicker):
    def __init__(self):
        pass
    
    def __call__(self, grid, bnd_in):
        bnd_points = grid.bnd_points()
        
        inds = self.find_nearest_points(bnd_points, bnd_in)
       
        # Gather all the data from one point into one Dataset
        bnd_out, bnd_mask = slice_and_gather_xr(inds, bnd_in)
        return bnd_out, bnd_mask


    def find_nearest_points(self, bnd_points, bnd_in):
        lon = bnd_points[:,0]
        lat = bnd_points[:,1]
        
        lon_vec, lat_vec = lon_lat_from_dataset(bnd_in[0])
        
        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = min_distance(lon[n], lat[n], lon_vec, lat_vec)
            print(f"Point {n}: lat: {lat[n]}, lon: {lon[n]} <<< ({lat_vec[ind]: .4f}, {lon_vec[ind]: .4f}). Distance: {dx:.1f} km")
            inds.append(ind) 
        return inds
    
    
class PPAreaPicker(PointPicker):
    def __init__(self, expansion_factor = 1.5):
        self.expansion_factor = expansion_factor
        return

    def __call__(self, grid, bnd_in):
        
        inds = self.find_inds_inside_area(grid, bnd_in)
        
        # Gather all the data from one point into one Dataset
        bnd_out, bnd_mask = slice_and_gather_xr(inds, bnd_in)
        
        return bnd_out, bnd_mask

    def find_inds_inside_area(self, grid, bnd_in):
         # Define area to search in
        expand_lon = (grid.lon_max - grid.lon_min)*self.expansion_factor*0.5
        expand_lat = (grid.lat_max - grid.lat_min)*self.expansion_factor*0.5
        
        # Get all the spectra in this area
        lon0=grid.lon_min - expand_lon
        lon1=grid.lon_max + expand_lon
        
        lat0=grid.lat_min - expand_lat
        lat1=grid.lat_max + expand_lat

        masklon = np.logical_and(bnd_in[0].longitude.values < lon1, bnd_in[0].longitude.values > lon0)
        masklat = np.logical_and(bnd_in[0].latitude.values > lat0, bnd_in[0].latitude.values < lat1)
        mask=np.logical_and(masklon, masklat)
        
        inds = bnd_in[0].x.values[mask[0]]
        
        return inds

class PPLegacyPicker(PointPicker):
    def __init__(self):
        pass
    
    def __call__(self, grid, bnd_in):
        self.bnd_in = bnd_in
        bnd_points = grid.bnd_points()
        
        inds = self.find_nearest_points(bnd_points)
        
        # Gather all the data from one point into one Dataset
        bnd_out = self.slice_and_gather_xr(inds)
        return bnd_out

    def slice_and_gather_xr(self, inds):
        bnd_out=[]
        for n in range(len(inds)): # Loop over output points
            datasets = []
            title = self.bnd_in[0].title # mean('x') looses all attributes so save this
            for k in range(len(self.bnd_in)): # Loop over time (days)
                new_data = self.bnd_in[k].sel(x=inds[n]).mean('x')
                new_data.attrs["title"] = title
                # This trivial dimension of one will cause probems since the latitudes are either then defined with only a y-dimension, or in two dimensions (time, y)
                # In the latter case new_data["latitude"].values[0] is not a number, but a one element array of a one element list
                # We solve this by squeezing out this unceccesary dimension
                datasets.append(new_data.squeeze('y')) 


            # This list contains one element for each requested output point
            # Each element is a one-spatial-point xr-dataset containing all times
            
            bnd_out.append(xr.concat(datasets, dim="time"))

        return bnd_out

    def find_nearest_points(self, bnd_points):
        nr_spec_interpolate = 1
        lat_in = self.lat_in()
        lon_in = self.lon_in()

        lat_out = bnd_points[0,:]
        lon_out = bnd_points[1,:]
        inds = []
        for n in range(len(lat_out)):
            diff_lat = np.abs(lat_out[n]-lat_in)
            diff_lon = np.abs(lon_out[n]-lon_in)
            add_lonlat = diff_lon + diff_lat
            inds.append(np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[nr_spec_interpolate]))[0])
            
        # The index in Xarray datasets starts from 1
        inds = inds + 1
        return inds
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# INPUT MODEL CLASSES RESPONSIBLE FOR ACTUALLY READING THE SPECTRA
# -----------------------------------------------------------------------------
# The output of an InputModel class should be a list of N xarray datasets
#
# Each dataset in the list should contain only one point.
# Each dataset in the list should contain all time steps for said point.
#
# Not mandatory to use, but might help:
#   The function slice_and_gather_xr goes from daily datasets containing
#   several points to datasets vontaining one point and all times.
#   A list of N indices must be provided to choose the points.
#
# Each dataset must contain the following Coordinates:
# -----------------------------------------------------------------------------
# time (N,) datetime64 array
# freq (k,) float array
# direction (d,) float array
# Each dataset must contain the following Data variable:
# -----------------------------------------------------------------------------
# SPEC (N,k,d) float array
# longitude 1-dimensional array (use .flatten() if 0-dimension)
# latitude 1-dimensional array (use .flatten() if 0-dimension)
# -----------------------------------------------------------------------------
# longitude/latitude should be constant.
# In case length > 1, only first value is ever used.
# -----------------------------------------------------------------------------

class InputModel(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def __call__(self, start_time, end_time):
        pass

    def days(self):
        """Determins a Pandas data range of all the days in the time span of the InputModel objext"""
        days = pd.date_range(start=self.start_time.split('T')[0], end=self.end_time.split('T')[0], freq='D')
        return days

    def get_time_limits(self, ind):
        """Determines star and end time for the day. First and last day doesn't start at 00:00 or end at 23:59"""
        if ind == 0:
            t0 = self.start_time
            t1 = self.days()[0].strftime('%Y-%m-%d') + "T23:59:59"
        elif ind == (len(self.days())-1):
            t0 = self.days()[-1].strftime('%Y-%m-%d') + "T00:00:00"				
            t1 = self.end_time
        else:
            t0 = self.days()[ind].strftime('%Y-%m-%d') + "T00:00:00"	
            t1 = self.days()[ind].strftime('%Y-%m-%d') + "T23:59:59"
        return t0, t1

    def __str__(self):
        return (f"{self.start_time} - {self.end_time}")


class InputWAM4(InputModel):
    
    def __call__(self, start_time, end_time, grid, point_picker = PPTrivialPicker()):
        self.start_time = start_time
        self.end_time = end_time
        
        msg.header(f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
            
        bnd_in, bnd_mask = point_picker(grid, bnd_in)
         
        # longitude and latitude are 0-dimensional. Force to one-dimensional trivial array
        for n in range(len(bnd_in)):
            bnd_in[n]["latitude"] = bnd_in[n]["latitude"].values.flatten()
            bnd_in[n]["longitude"] = bnd_in[n]["longitude"].values.flatten()
        
        return bnd_in, bnd_mask

    def get_url(self, ind):
        day = self.days()[ind]
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        print(url)
        
        return url
    

class InputNORA3(InputModel):
    
    def __call__(self, start_time, end_time, grid, point_picker = PPTrivialPicker()):
        self.start_time = start_time
        self.end_time = end_time
        msg.header(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
        
        bnd_in, bnd_mask = point_picker(grid, bnd_in)
           
        return bnd_in, bnd_mask

    def get_url(self, ind):
        day = self.days()[ind]
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        print(url)
        return url
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# OUTPUT MODEL CLASSES RESPONSIBLE FOR WRITING SPECTRA IN CORRECT FORMAT
# -----------------------------------------------------------------------------
class OutputModel(ABC):
    bnd_in: list
    bnd_points: np.array
    message: str
    @abstractmethod
    def __call__(self, bnd_out):
        pass

    def spec_info(self, bnd_out):
        dirs=bnd_out[0].direction.values
        freq=bnd_out[0].freq.values
        return freq, dirs
  
    
class OutputWW3nc(OutputModel):
    def __init__(self):
        pass
                
    def __call__(self, bnd_out, bnd_mask):
        # Convert from oceanic to mathematical convention
        #dirs = bnd_out[0].direction.values
        for n in range(len(bnd_out)):
            #dirs = bnd_out[n].direction.values
            D = bnd_out[n].direction.values
            for k in range(len(bnd_out[0].time.values)):
                S = naut_to_ocean(bnd_out[n].SPEC[k,:,:].values, D)
                bnd_out[n].SPEC[k,:,:] = S
                
            D = ocean_to_math(D,D)
            bnd_out[n] = bnd_out[n].assign_coords(direction=D)


        msg.header('Writing WAVEWATCH-III netcdf-output')

        for n in range(len(bnd_out)):
            if bnd_mask[n]:
                print(f"Point {n}: ", end="")
                #self.write_netcdf(bnd_out[n], output_points[n,:])
                self.write_netcdf(bnd_out[n])
            else:
                msg.info(f"Skipping point {n}. Masked as False.")
        return
    
    def write_netcdf(self, bnd_out):
        """Writes WW3 compatible netcdf spectral output from a list containing xarray datasets."""
        lat=bnd_out["latitude"].values[0]
        lon=bnd_out["longitude"].values[0] 
        
        #lat = output_points[0]
        #lon = output_points[1]
        output_file = f"ww3_spec_E{lon:09.6f}N{lat:09.6f}.nc"
        #output_file = 'ww3_spec_E'+str(lon)+'N'+str(lat)+'.nc'
        #output_file = 'Test_ww3.nc'
        print(output_file)
        root_grp = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
        #################### dimensions
        root_grp.createDimension('time', None)
        root_grp.createDimension('station', 1)
        root_grp.createDimension('string16', 16)
        root_grp.createDimension('frequency', len(bnd_out.freq))
        root_grp.createDimension('direction', len(bnd_out.direction))
        #######################################################
        ####################### variables
        time = root_grp.createVariable('time', np.float64, ('time',))
        station = root_grp.createVariable('station', np.int32, ('station',))
        frequency = root_grp.createVariable('frequency',np.float32 , ('frequency',))
        direction = root_grp.createVariable('direction', np.float32, ('direction',))
        efth = root_grp.createVariable('efth', np.float32, ('time','station','frequency','direction',))
        latitude = root_grp.createVariable('latitude',np.float32 , ('time','station',))
        longitude = root_grp.createVariable('longitude',np.float32 , ('time','station',))
        station_name = root_grp.createVariable('station_name', 'S1', ('station','string16',))
        string16 = root_grp.createVariable('string16',np.int32 , ('string16',))
               
        ########################## Attributes
        time.units = 'seconds since 1970-01-01 00:00:00 UTC'
        time.calendar = "standard"
        time.standard_name = "time" 
        time.axis = "T" 
        
        station.long_name = "station id" 
        station.axis = "X" 
        
        frequency.units = "s-1" 
        frequency.long_name = "frequency of center band" 
        frequency.standard_name = "sea_surface_wave_frequency" 
        frequency.globwave_name = "frequency" 
        frequency.valid_min = 0 
        frequency.valid_max = 10 
        frequency.axis = "Y" 
        
        direction.units = "degree" 
        direction.long_name = "sea surface wave to direction" 
        direction.standard_name = "sea_surface_wave_to_direction" 
        direction.globwave_name = "direction" 
        direction.valid_min = 0
        direction.valid_max = 360
        direction.axis = "Z" 
        
        longitude.units='degree_east'
        longitude.long_name = "longitude" 
        longitude.standard_name = "longitude" 
        longitude.valid_min = -180
        longitude.valid_max = 180
        	#longitude:_FillValue = 9.96921e+36f ;
        longitude.content = "TX" 
        longitude.associates = "time station" 
        
        latitude.units = "degree_north" 
        latitude.long_name = "latitude" 
        latitude.standard_name = "latitude" 
        latitude.valid_min = -90
        latitude.valid_max = 90
        	#latitude:_FillValue = 9.96921e+36f ;
        latitude.content = "TX" 
        latitude.associates = "time station"
        
        station_name.long_name = "station name" 
        station_name.content = "XW" 
        station_name.associates = "station string16" 
        
        station.long_name = "station id" 
        station.axis = "X" 
        
        string16.long_name = "station_name number of characters" 
        string16.axis = "W" 
        
        efth.long_name = "sea surface wave directional variance spectral density" 
        efth.standard_name = "sea_surface_wave_directional_variance_spectral_density" 
        efth.globwave_name = "directional_variance_spectral_density" 
        efth.units = "m2 s rad-1" 
        efth.scale_factor = 1 
        efth.add_offset = 0 
        efth.valid_min = 0 
        #efth.valid_max = 1.0e+20 
        #efth._FillValue = 9.96921e+36 
        efth.content = "TXYZ" 
        efth.associates = "time station frequency direction" 
        #######################################################
        ############## Pass data
        time[:] = bnd_out.time.values.astype('datetime64[s]').astype('float64')
        frequency[:] = bnd_out.freq.values
        direction[:] = bnd_out.direction.values 
        #dtheta=np.radians(np.diff(direction).mean())
        #print("NOT CONVERTING TO OCEANIC CONVENTION!")
        
        efth[:] =  bnd_out.SPEC.values
        station[:] = 1
        longitude[:] = np.full((len(bnd_out.time),1), lon,dtype=float)
        latitude[:] = np.full((len(bnd_out.time),1), lat,dtype=float)
        #longitude[:] = bnd_out.longitude.values
        #latitude[:] = bnd_out.latitude.values
        station_name[:] = 1
        
        root_grp.close() 
        return
   
   
class OutputSWANascii(OutputModel):
    def __init__(self, grid, factor = 1E-4):
        self.factor = factor
        self.grid = grid
                
    def __call__(self, bnd_in, bnd_mask):
        output_points = self.grid.bnd_points()
        # Initialize the boundary file by writing the header
        
        msg.header('Writing SWAN ASCII-output')
        
        freq, dirs = self.spec_info(bnd_in)
        with open('outfile', 'w') as file_out:
            file_out.write('SWAN   1\n')
            file_out.write('$ Data produced by '+bnd_in[0].title+'\n')
            file_out.write('TIME\n')
            file_out.write('          1\n')
            file_out.write('LONLAT\n')    
            file_out.write('          '+format(len(bnd_in))+'\n')     
            for k in range(len(bnd_in)):
                file_out.write('   '+format(output_points[k,1],'.4f')+'  '+format(output_points[k,0],'.4f')+'\n')
            file_out.write('AFREQ\n')
            file_out.write('          '+str(len(freq))+'\n')
            #breakpoint()
            for l in range(len(freq)):
                file_out.write('   '+format(freq[l],'.4f')+'\n')
            file_out.write('NDIR\n')
            file_out.write('          '+format(len(dirs))+'\n')
            for m in range(len(dirs)):
                file_out.write('   '+format(dirs[m],'.1f')+'\n') 
            file_out.write('QUANT\n')
            file_out.write('          1\n')
            file_out.write('VaDens\n')
            file_out.write('m2/Hz/degr \n')
            file_out.write('-32767\n')
                #first day
            msg.info('Generating 2d spectra at boundaries:')
    
            with open('outfile', 'w') as file_out:
                times = pd.DatetimeIndex(bnd_in[0]["time"].values) # All point have the same time vector so use first Dataset
                days = pd.date_range(start=times[0], end=times[-1], freq='D')
                
                for d in range(len(days)):
                    msg.plain(days[d].strftime('%Y-%m-%d'))
                    day_inds = np.where(times.day == days[d].day)[0]
                    
                    for time_step in day_inds:
                        time_stamp = str(times[time_step]).split('-')[0]+str(times[time_step]).split('-')[1]+\
                        str(times[time_step]).split('-')[2][:2]+'.'+str(times[time_step]).split('-')[2][3:5]+'0000\n'
                        file_out.write(time_stamp)
                        
                        for i in range(len(bnd_in)):
                            file_out.write('FACTOR\n')
                            file_out.write(format(self.factor,'1.0E')+'\n')
                            SPEC_ocean_convection = bnd_in[i].SPEC[time_step,:,:].values
                            SPEC_naut_convection = ocean_to_naut(SPEC_ocean_convection)
                            delth = 360/len(dirs)
                            np.savetxt(file_out,SPEC_naut_convection/(delth*self.factor), fmt='%-10.0f') #


        return
# -----------------------------------------------------------------------------