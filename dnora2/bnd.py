import numpy as np
import pandas as pd
import xarray as xr
from abc import ABC, abstractmethod

class SpectraToGrid(ABC):
    @abstractmethod
    def fit_spectra_to_grid(self):
        pass

    def lat_in(self):
        return self.bnd_points[:, 0]
    
    def lon_in(self):
        return self.bnd_points[:, 1]
    
    def lat_out(self):
        return self.bnd_in[0]["latitude"][0].values 
    
    def lon_out(self):
        return self.bnd_in[0]["longitude"][0].values 

    
class NearestSpectraToGrid(SpectraToGrid):
    def __init__(self, bnd_in, bnd_points):
        self.bnd_in = bnd_in
        self.bnd_points = bnd_points
        pass
    
    def fit_spectra_to_grid(self):
        inds = self.find_nearest_points()
       
        # Gather all the data from one point into one Dataset
        bnd_out = self.slice_and_gather_xr(inds)
        return bnd_out

    def slice_and_gather_xr(self, inds):
        bnd_out=[]
        for n in range(len(inds)): # Loop over output points
            datasets = []
            for k in range(len(self.bnd_in)): # Loop over time (days)
                datasets.append(self.bnd_in[k].sel(x=inds[n]))
             
            # This list contains one element for each requested output point
            # Each element is a one-spatial-point xr-dataset containing all times
            bnd_out.append(xr.concat(datasets, dim="time"))        
        
        return bnd_out    


    def find_nearest_points(self):
        lat = self.lat_in()
        lon = self.lon_in()
        
        lat_vec = self.lat_out()
        lon_vec = self.lon_out()
        
        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = self.min_distance(n)
            print(f"Point {n}: lat: {lat[n]}, lon: {lon[n]} <<< ({lat_vec[ind]: .4f}, {lon_vec[ind]: .4f}). Distance: {dx:.1f} km")
            inds.append(ind) 
        return inds
    
    
    def min_distance(self, ind):
        lat = self.lat_in()[ind]
        lon = self.lon_in()[ind]
        
        # Longitudes seems to be a trivial list of list, therefore the extre [0] before taking .values    
        lat_vec = self.lat_out()
        lon_vec = self.lon_out()
        

        dx = []
        for n in range(len(lat_vec)):
            dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))
        
        
        return np.array(dx).min(), np.array(dx).argmin()




class InputModel(ABC):
    start_time: str
    end_time: str
	
    @abstractmethod
    def read_bnd(self):
        pass

    def days(self):
        days = pd.date_range(start=self.start_time.split('T')[0], end=self.end_time.split('T')[0], freq='D')
        return days

    def get_time_limits(self, ind):
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

    def check_output(self, bnd_out):
        print ("Check that the xr dataset is of the right format here...")
        return bnd_out


class InputWAM4(InputModel):
    def __init__(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
    		
    def read_bnd(self):
        print(f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
        return bnd_in

    def get_url(self, ind):
        day = self.days()[ind]
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        
        return url
    

class InputNORA3(InputModel):
    def __init__(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
    		
    def read_bnd(self):
        print(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
        return bnd_in

    def get_url(self, ind):
        day = self.days()[ind]
        url = 'https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/SPC'+day.strftime('%Y%m%d')+'00.nc'
        
        return url


class OutputModel(ABC):
    bnd_in: list
    bnd_points: np.array
    message: str
    @abstractmethod
    def output_spec(self):
        pass

    def spec_dim(self):
        dirN=self.bnd_in[0].direction.shape[0]
        freqN=self.bnd_in[0].freq.shape[0]
        return freqN, dirN
    
    


    
class OutputWW3nc(OutputModel):
    def __init__(self, bnd_out):
        self.bnd_out = bnd_out
                
    def output_spec(self):
        # For WW3-netcdf's we need to write one netcdf for each month for each location
        # The input data is given one xr-dataset per day, containing all the points

            
        # Change the metadata so that it matches WW3 requirements (does nothing for now)
        ww3_bnd_out = self.to_ww3_metadata(self.bnd_out)
        
        # Write netcdf-output
        self.write_netcdf(ww3_bnd_out)
        
        return
    
    def write_netcdf(self, xrlist):
        for n in range(len(xrlist)):
            fn = 'Test' + str(n) + '.nc'
            xrlist[n].to_netcdf(fn)    
        return
    
    def to_ww3_metadata(self, xrlist):
        return xrlist
    
    
   
   
class OutputSWANascii(OutputModel):
    def __init__(self, bnd_in, bnd_points):
        self.factor = 1E-4
        self.nr_spec_interpolate = 0
        self.calib_spec = 1
        self.bnd_in = bnd_in
        self.bnd_points = bnd_points
	
    def output_spec(self):
        # Initialize the boundary file by writing the header

        freqN, dirN = self.spec_dim()
        
        with open('outfile', 'w') as file_out:
            file_out.write('SWAN   1\n')
            file_out.write('$ Data produced by '+self.bnd_in[0].title+'\n')
            file_out.write('TIME\n')
            file_out.write('          1\n')
            file_out.write('LONLAT\n')    
            file_out.write('          '+format(self.bnd_points.shape[0])+'\n')     
            for k in range((self.bnd_points.shape[0])):
                file_out.write('   '+format(self.bnd_points[k][1],'.4f')+'  '+format(self.bnd_points[k][0],'.4f')+'\n')
            file_out.write('AFREQ\n')
            file_out.write('          '+str(self.bnd_in[0].freq.shape[0])+'\n')
            for l in range(freqN):
                file_out.write('   '+format(self.bnd_in[0].freq[l].values,'.4f')+'\n')
            file_out.write('NDIR\n')
            file_out.write('          '+format(self.bnd_in[0].direction.shape[0])+'\n')
            for m in range(dirN):
                file_out.write('   '+format(self.bnd_in[0].direction[m].values,'.1f')+'\n') 
            file_out.write('QUANT\n')
            file_out.write('          1\n')
            file_out.write('VaDens\n')
            file_out.write('m2/Hz/degr \n')
            file_out.write('-32767\n')
                #first day
            print('Generating 2d spectra at boundaries:')
        
            for i in range(len(self.bnd_in)):
                first_time = pd.to_datetime(self.bnd_in[i]["time"].values[0])
                print(first_time.strftime('%Y-%m-%d'))
                for time_step in range(self.bnd_in[i].time.shape[0]):
                    time_stamp = str(self.bnd_in[i].time.time[time_step].values).split('-')[0]+str(self.bnd_in[i].time.time[time_step].values).split('-')[1]+\
                    str(self.bnd_in[i].time.time[time_step].values).split('-')[2][:2]+'.'+str(self.bnd_in[i].time.time[time_step].values).split('-')[2][3:5]+'0000\n'
                    file_out.write(time_stamp)
                    for p in range(self.bnd_points.shape[0]):
                        file_out.write('FACTOR\n')    
                        file_out.write(format(self.factor,'1.0E')+'\n')
                        # find the nearest grid point
                        diff_lat = np.abs(self.bnd_points[p][0]-self.bnd_in[i].latitude[0,:])
                        diff_lon = np.abs(self.bnd_points[p][1]-self.bnd_in[i].longitude[0,:])
                        add_lonlat = diff_lon + diff_lat
                        index_min_dinstance = np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[self.nr_spec_interpolate]))[0]
                        SPEC_ocean_convection = self.bnd_in[i].SPEC[time_step,0,index_min_dinstance,:,:].mean('x').values
                        SPEC_naut_convection = ocean_to_naut(SPEC_ocean_convection)
                        delth = 360/dirN
                        np.savetxt(file_out,self.calib_spec*SPEC_naut_convection/(delth*self.factor), fmt='%-10.0f') #     
        return


def ocean_to_naut(oceanspec):
    nautspec=np.zeros(oceanspec.shape)
    dirN=oceanspec.shape[1]
    nautspec[:,0:dirN//2] = oceanspec[:,dirN//2:] # Step 1a: 180..355 to start of array
    nautspec[:,dirN//2:]  = oceanspec[:,0:dirN//2] # Step 1b: 0..175 to end of array
    
    return nautspec

def distance_2points(lat1,lon1,lat2,lon2):
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
