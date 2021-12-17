from __future__ import annotations # For TYPE_CHECKING

import numpy as np
from copy import copy
from abc import ABC, abstractmethod
import netCDF4
import re

# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .bnd_mod import Boundary

# Import default values and auxiliry functions
from .. import msg
from ..defaults import dflt_bnd, list_of_placeholders
from ..aux import add_folder_to_filename, clean_filename, create_filename_lonlat

class BoundaryWriter(ABC):
    """Writes the boundary spectra to a certain file format.

    This object is provided to the .export_boundary() method.
    """

    def _preferred_format(self) -> str:
        """For the file format using defauts.py, e.g. SWAN or WW3"""
        return 'General'

    def _preferred_extension(self):
        return 'nc'


    @abstractmethod
    def convention(self) -> str:
        """Defines in which format the incoming spectra should be.

        The conventions to choose from are predetermined:

        'Ocean':    Oceanic convention
                    Directional vector monotonically increasing.
                    Direction to. North = 0, East = 90.

        'Met':      Meteorological convention
                    Directional vector monotonically increasing.
                    Direction from. North = 0, East = 90.

        'Math':     Mathematical convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 90, East = 0.

        'WW3':      WAVEWATCH III output convention
                    Directional vector of type: [90 80 ... 10 0 350 ... 100]
                    Direction to. North = 0, East = 90.
        """
        return

    @abstractmethod
    def __call__(self, boundar: Boundary, filename: str, folder: str) -> Tuple[str, str]:
        """Writed the data from the Boundary object and returns the file and
        folder where data were written."""

        return

class DumpToNc(BoundaryWriter):
    def __init__(self, convention: str='Ocean', out_format: str='General') -> None:
        self._convention = convention
        self.out_format = out_format
        return

    def convention(self) -> str:
        """Convention of spectra"""
        return self._convention

    def _preferred_format(self) -> str:
        """Preferred format of file name"""
        return self.out_format

    def __call__(self, boundary: Boundary, filename: str, folder: str) -> Tuple[str, str]:

        output_file = clean_filename(filename, list_of_placeholders)

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=folder)

        # Dumping to a netcdf-file
        msg.to_file(output_path)
        boundary.data.to_netcdf(output_path)

        return output_file, folder


class NcFiles(BoundaryWriter):
    def __init__(self, convention: str='Ocean') -> None:
        self._convention = convention
        return

    def convention(self) -> str:
        """Convention of spectra"""
        return self._convention


    def __call__(self, boundary: Boundary, filename: str, folder: str) -> Tuple[str, str]:

        output_files = []
        for n in boundary.x():
            output_file = create_filename_lonlat(filename, lon=boundary.lon()[n], lat=boundary.lat()[n])
            output_file = clean_filename(output_file, list_of_placeholders)
            output_files.append(output_file)

            output_path = add_folder_to_filename(output_file, folder=folder)
            msg.to_file(output_path)

            ds = boundary.slice_data(x=[n])
            ds.to_netcdf(output_path)

        return output_files, folder


class WW3(BoundaryWriter):
    def __init__(self, one_file: bool=True, convention='WW3', out_format = 'WW3') -> None:
        self.one_file = one_file
        self._convention=convention
        self.out_format = out_format
        return

    def convention(self) -> str:
        """Convention of spectra"""
        return self._convention

    def _preferred_format(self):
        """Format of filenames"""
        return self.out_format

    def __call__(self, boundary: Boundary, filename: str, folder: str) -> Tuple[str, str]:


        #boundary_in.change_convention(wanted_convention='WW3')

        msg.info('Writing WAVEWATCH-III netcdf-output')



        if self.one_file:
            if len(boundary.x()) == 1:
                # Uses $Lon $Lat
                filename = create_filename_lonlat(filename, lon=boundary.lon()[0], lat=boundary.lat()[0])

            output_files = clean_filename(filename, list_of_placeholders)

            output_path = add_folder_to_filename(output_files, folder=folder)
            msg.plain(f"All points >> {output_path}")
            self.write_netcdf(boundary, output_path)

        else:
            output_files = []
            for n in boundary.x():
                if boundary.mask[n]: # This property is not really used and should always be true
                    output_file = create_filename_lonlat(filename, lon=boundary.lon()[n], lat=boundary.lat()[n])
                    output_file = clean_filename(output_file, list_of_placeholders)

                    output_files.append(output_file)
                    output_path = add_folder_to_filename(output_file, folder=folder)
                    msg.plain(f"Point {n} >> {output_path}")
                    self.write_netcdf(boundary, output_path, n)
                else:
                    msg.info(f"Skipping point {n} ({boundary.lon()[n]:10.7f}, {boundary.lat()[n]:10.7f}). Masked as False.")

        return output_files, folder

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
    def __init__(self, factor = 1E-4, out_format = 'SWAN') -> None:
        self.factor = factor
        self.out_format = out_format
        return

    def convention(self) -> str:
        """Convention of spectra"""
        return 'Met'

    def _preferred_format(self):
        """Format of filenames"""
        return self.out_format

    def _preferred_extension(self):
        return 'asc'


    def __call__(self, boundary: Boundary, filename: str, folder: str) -> Tuple[str, str]:
        #boundary_in = copy(boundary)

        output_file = clean_filename(filename, list_of_placeholders)
        output_path = add_folder_to_filename(filename, folder=folder)

        #boundary_in.change_convention(wanted_convention='Ocean')

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
            msg.to_file(f"{output_path}")

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
                        # SWAN uses m*m/Hz/deg normalization
                        np.savetxt(file_out,S*np.pi/(180*self.factor), fmt='%-10.0f')
        return output_file, folder
