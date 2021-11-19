import numpy as np
from copy import copy
from .. import msg
from ..aux import check_if_folder, create_filename_obj, create_filename_time, create_filename_lonlat, add_file_extension, add_folder_to_filename
import netCDF4
import re

from .bnd_mod import BoundaryWriter # Abstract class
from .bnd_mod import Boundary # Boundary object

from .process import OceanToWW3
from ..defaults import dflt_bnd

class DumpToNc(BoundaryWriter):
    def __init__(self, folder: str='', filestring: str=dflt_bnd['fs']['General'], datestring: str=dflt_bnd['ds']['General']) -> None:
        self.folder = copy(folder)

        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, boundary: Boundary) -> None:
        msg.header(f'{type(self).__name__}: writing boundary spectra from {boundary.name()}')


        output_file = boundary.filename(filestring=self.filestring, datestring=self.datestring, extension='nc')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)

        # Dumping to a netcdf-file
        msg.to_file(output_path)
        boundary.data.to_netcdf(output_path)

        return output_file, self.folder


class NcFiles(BoundaryWriter):
    def __init__(self, folder: str='', filestring: str=dflt_bnd['fs']['General'], datestring: str=dflt_bnd['ds']['General']) -> None:
        self.folder = copy(folder)

        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, boundary: Boundary):
        msg.header(f'{type(self).__name__}: writing boundary spectra from {boundary.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_files = []
        for n in boundary.x():
            output_file = boundary.filename(filestring=self.filestring, datestring=self.datestring, n=n, extension='nc')
            output_files.append(output_file)

            output_path = add_folder_to_filename(output_files[n], folder=self.folder)
            msg.to_file(output_path)

            ds = boundary.slice_data(x=[n])
            ds.to_netcdf(output_path)

        return output_files, self.folder


class WW3(BoundaryWriter):
    def __init__(self, folder: str='', one_file: bool=True, filestring: str=dflt_bnd['fs']['WW3'], datestring: str=dflt_bnd['ds']['WW3']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        self.one_file = one_file
    def __call__(self, boundary: Boundary) -> None:

        boundary_in = copy(boundary)
        msg.header(f'{type(self).__name__}: writing boundary spectra from {boundary_in.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        # Convert from oceanic to mathematical convention
        boundary_in.process_spectra(OceanToWW3())

        msg.info('Writing WAVEWATCH-III netcdf-output')



        if self.one_file:
            if len(boundary_in.x()) == 1:
                # Uses $Lon $Lat
                output_files = boundary_in.filename(filestring=self.filestring, datestring=self.datestring, n=0, extension='nc')
            else:
                # Tries to remove $Lon $Lat in filestring
                output_files = boundary_in.filename(filestring=self.filestring, datestring=self.datestring, extension='nc')

            output_path = add_folder_to_filename(output_files, folder=self.folder)
            msg.plain(f"All points >> {output_path}")
            self.write_netcdf(boundary_in, output_path)

        else:
            output_files = []
            for n in boundary_in.x():
                if boundary_in.mask[n]: # This property is not really used and should always be true
                    output_file = boundary_in.filename(filestring=self.filestring, datestring=self.datestring, n=n, extension='nc')
                    output_files.append(output_file)
                    output_path = add_folder_to_filename(output_file, folder=self.folder)
                    msg.plain(f"Point {n} >> {output_path}")
                    self.write_netcdf(boundary_in, output_path, n)
                else:
                    msg.info(f"Skipping point {n} ({boundary_in.lon()[n]:10.7f}, {boundary_in.lat()[n]:10.7f}). Masked as False.")

        return output_files, self.folder

    def write_netcdf(self, boundary: Boundary, output_file: str, n: int=None) -> None:
        """Writes WW3 compatible netcdf spectral output from a list containing xarray datasets."""

        #if boundary.name == "AnonymousBoundary":
        #    output_file = f"ww3_spec_E{lon:09.6f}N{lat:09.6f}.nc"
        #else:
        #    output_file = f"ww3_{boundary.name}_E{lon:09.6f}N{lat:09.6f}.nc"
        #output_file = 'ww3_spec_E'+str(lon)+'N'+str(lat)+'.nc'
        #output_file = 'Test_ww3.nc'
        root_grp = netCDF4.Dataset(output_file, 'w', format='NETCDF4')
        #################### dimensions
        root_grp.createDimension('time', None)
        if self.one_file:
            root_grp.createDimension('station', len(boundary.x()))
        else:
            root_grp.createDimension('station', 1)
        root_grp.createDimension('string16', 16)
        root_grp.createDimension('frequency', len(boundary.freq()))
        root_grp.createDimension('direction', len(boundary.dirs()))

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
        time[:] = boundary.time().values.astype('datetime64[s]').astype('float64')
        frequency[:] =boundary.freq()
        direction[:] = boundary.dirs()

        if self.one_file:
            station[:] = boundary.x()
            efth[:] =  boundary.spec()
            longitude[:] = np.full((len(boundary.time()),len(boundary.lon())), boundary.lon(),dtype=float)
            latitude[:] = np.full((len(boundary.time()),len(boundary.lat())), boundary.lat(),dtype=float)
        else:
            station[:] = 1
            efth[:] =  boundary.spec(x=[n])
            longitude[:] = np.full((len(boundary.time()),1), boundary.lon()[n],dtype=float)
            latitude[:] = np.full((len(boundary.time()),1), boundary.lat()[n],dtype=float)
        #longitude[:] = bnd_out.longitude.values
        #latitude[:] = bnd_out.latitude.values
        station_name[:] = 1

        root_grp.close()
        return


class SWAN(BoundaryWriter):
    #def __init__(self, factor = 1E-4, folder: str='', boundary_in_filename: bool=True, time_in_filename: bool=True, grid_in_filename: bool=True) -> None:
    def __init__(self, factor = 1E-4, folder: str='', filestring: str=dflt_bnd['fs']['SWAN'], datestring: str=dflt_bnd['ds']['SWAN']) -> None:
        self.factor = factor

        self.folder = copy(folder)

        self.filestring = copy(filestring)
        self.datestring = copy(datestring)
        return

    def __call__(self, boundary: Boundary):
        msg.header(f'{type(self).__name__}: writing boundary spectra from {boundary.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = boundary.filename(filestring=self.filestring, datestring=self.datestring, extension='asc')
        output_path = add_folder_to_filename(output_file, folder=self.folder)

        swan_bnd_points = boundary.grid.boundary_points()
        days = boundary.days()

        with open(output_path, 'w') as file_out:
            file_out.write('SWAN   1\n')
            file_out.write('$ Data produced by '+boundary.data.source+'\n')
            file_out.write('TIME\n')
            file_out.write('          1\n')
            file_out.write('LONLAT\n')
            file_out.write('          '+format(len(boundary.x()))+'\n')
            for k in range(len(boundary.x())):
                file_out.write('   '+format(swan_bnd_points[k,0],'.4f')+'  '+format(swan_bnd_points[k,1],'.4f')+'\n')
            file_out.write('AFREQ\n')
            file_out.write('          '+str(len(boundary.freq()))+'\n')
            for l in range(len(boundary.freq())):
                file_out.write('   '+format(boundary.freq()[l],'.4f')+'\n')
            file_out.write('NDIR\n')
            file_out.write('          '+format(len(boundary.dirs()))+'\n')
            for m in range(len(boundary.dirs())):
                file_out.write('   '+format(boundary.dirs()[m],'.1f')+'\n')
            file_out.write('QUANT\n')
            file_out.write('          1\n')
            file_out.write('VaDens\n')
            file_out.write('m2/Hz/degr \n')
            file_out.write('-32767\n')
                #first day
            msg.info(f'Writing 2d spectra at boundaries to: {output_path}')

            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = boundary.times_in_day(day)
                for tim in times:
                    time_stamp = str(tim).split('-')[0]+str(tim).split('-')[1]+str(tim).split('-')[2][:2]+'.'+str(tim).split('-')[2][3:5]+'0000\n'
                    file_out.write(time_stamp)
                    for n in range(len(boundary.x())):
                        file_out.write('FACTOR\n')
                        file_out.write(format(self.factor,'1.0E')+'\n')
                        S = boundary.spec(start_time=tim, end_time=tim, x=[n]).squeeze()
                        delth = 360/len(boundary.dirs())

        return output_file, self.folder
