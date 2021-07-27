import numpy as np
import pandas as pd
import xarray as xr
from abc import ABC, abstractmethod

class PointPicker(ABC):
    @abstractmethod
    def pick_points(self):
        pass


    def lat_out(self):
        return self.bnd_points[:, 0]
    
    def lon_out(self):
        return self.bnd_points[:, 1]
    
    def lat_in(self):
        return self.bnd_in[0]["latitude"][0].values 
    
    def lon_in(self):
        return self.bnd_in[0]["longitude"][0].values 

class PPTrivialPicker(PointPicker):
    def __init__(self):
        pass

    def pick_points(self, bnd_in):
        return bnd_in

class PPLegacyPicker(PointPicker):
    def __init__(self, bnd_points):
        self.bnd_points = bnd_points
        return
    
    def pick_points(self, bnd_in):
        self.bnd_in = bnd_in
        inds = self.find_nearest_points()

        # Gather all the data from one point into one Dataset
        bnd_out = self.slice_and_gather_xr(inds)
        return bnd_out

    def slice_and_gather_xr(self, inds):
        bnd_out=[]
        for n in range(len(inds)): # Loop over output points
            datasets = []
            title = self.bnd_in[0].title # mean('x') looses all attributes so save this
            for k in range(len(self.bnd_in)): # Loop over time (days)
                addition = self.bnd_in[k].sel(x=inds[n]).mean('x')
                addition.attrs["title"] = title
                datasets.append(addition)


            # This list contains one element for each requested output point
            # Each element is a one-spatial-point xr-dataset containing all times
            bnd_out.append(xr.concat(datasets, dim="time"))

        return bnd_out

    def find_nearest_points(self):
        nr_spec_interpolate = 1
        lat_in = self.lat_in()
        lon_in = self.lon_in()

        lat_out = self.lat_out()
        lon_out = self.lon_out()
        inds = []
        for n in range(len(lat_out)):
            diff_lat = np.abs(lat_out[n]-lat_in)
            diff_lon = np.abs(lon_out[n]-lon_in)
            add_lonlat = diff_lon + diff_lat
            inds.append(np.where((add_lonlat >= np.sort(add_lonlat)[0]) & (add_lonlat <= np.sort(add_lonlat)[nr_spec_interpolate]))[0])
        return inds


class PPNearestGridPoint(PointPicker):
    def __init__(self, bnd_points):
        self.bnd_points = bnd_points
        return
    
    def pick_points(self, bnd_in):
        self.bnd_in = bnd_in
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
        lat = self.lat_out()
        lon = self.lon_out()
        
        lat_vec = self.lat_in()
        lon_vec = self.lon_in()
        
        # Go through all points where we want output and find the nearest available point
        inds = []
        for n in range(len(lat)):
            dx, ind = self.min_distance(n)
            print(f"Point {n}: lat: {lat[n]}, lon: {lon[n]} <<< ({lat_vec[ind]: .4f}, {lon_vec[ind]: .4f}). Distance: {dx:.1f} km")
            inds.append(ind) 
        return inds
    
    
    def min_distance(self, ind):
        lat = self.lat_out()[ind]
        lon = self.lon_out()[ind]
        
        lat_vec = self.lat_in()
        lon_vec = self.lon_in()
        

        dx = []
        for n in range(len(lat_vec)):
            dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))
        
        
        return np.array(dx).min(), np.array(dx).argmin()




class InputModel(ABC):
    def __init__(self, point_picker = None):
        if point_picker is not None:
            self.point_picker = point_picker
        else:
            self.point_picker = PPTrivialPicker()

    @abstractmethod
    def __call__(self, start_time, end_time):
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
    
    def __call__(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
        print(f"Getting boundary spectra from WAM4 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
            
        bnd_in = self.point_picker.pick_points(bnd_in)
        return bnd_in

    def get_url(self, ind):
        day = self.days()[ind]
        url = 'https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/'+day.strftime('%Y') +'/'+day.strftime('%m')+'/'+day.strftime('%d')+'/MyWave_wam4_SPC_'+day.strftime('%Y%m%d')+'T00Z.nc'
        
        return url
    

class InputNORA3(InputModel):
    
    def __call__(self, start_time, end_time):
        self.start_time = start_time
        self.end_time = end_time
        print(f"Getting boundary spectra from NORA3 from {self.start_time} to {self.end_time}")
        bnd_in = []    
        for n in range(len(self.days())):
            url = self.get_url(n)
            t0, t1 = self.get_time_limits(n)
            bnd_in.append(xr.open_dataset(url).sel(time=slice(t0, t1)))
        
        bnd_in = self.point_picker.pick_points(bnd_in)
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
    def __call__(self, bnd_in):
        pass

    def spec_info(self, bnd_in):
        dirs=bnd_in[0].direction.values
        freq=bnd_in[0].freq.values
        return freq, dirs
    
class OutputWW3nc(OutputModel):
    def __init__(self, bnd_out):
        self.bnd_out = bnd_out
                
    def __call__(self):
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
    def __init__(self):
        self.factor = 1E-4
        self.nr_spec_interpolate = 0
        self.calib_spec = 1
	
    def __call__(self, bnd_in):
        # Initialize the boundary file by writing the header

        freq, dirs = self.spec_info(bnd_in)
        with open('outfile', 'w') as file_out:
            file_out.write('SWAN   1\n')
            file_out.write('$ Data produced by '+bnd_in[0].title+'\n')
            file_out.write('TIME\n')
            file_out.write('          1\n')
            file_out.write('LONLAT\n')    
            file_out.write('          '+format(len(bnd_in))+'\n')     
            for k in range(len(bnd_in)):
                file_out.write('   '+format(bnd_in[k]["latitude"].values[0][0],'.4f')+'  '+format(bnd_in[k]["longitude"].values[0][0],'.4f')+'\n')
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
            print('Generating 2d spectra at boundaries:')
    
            with open('outfile', 'w') as file_out:
                times = pd.DatetimeIndex(bnd_in[0]["time"].values) # All point have the same time vector so use first Dataset
                days = pd.date_range(start=times[0], end=times[-1], freq='D')
                
                for d in range(len(days)):
                    print (days[d].strftime('%Y-%m-%d'))
                    day_inds = np.where(times.day == days[d].day)[0]
                    
                    for time_step in day_inds:
                        time_stamp = str(times[time_step]).split('-')[0]+str(times[time_step]).split('-')[1]+\
                        str(times[time_step]).split('-')[2][:2]+'.'+str(times[time_step]).split('-')[2][3:5]+'0000\n'
                        file_out.write(time_stamp)
                        
                        for i in range(len(bnd_in)):
                            file_out.write('FACTOR\n')
                            file_out.write(format(self.factor,'1.0E')+'\n')
                            SPEC_ocean_convection = bnd_in[i].SPEC[time_step,0,:,:].values
                            SPEC_naut_convection = ocean_to_naut(SPEC_ocean_convection)
                            delth = 360/len(dirs)
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
