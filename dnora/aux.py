from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy import interpolate
import os, re
from typing import TYPE_CHECKING, Tuple, List, Union

if TYPE_CHECKING:
    from .grd.grd_mod import Grid
    from .bnd.bnd_mod import Boundary
    from .wnd.wnd_mod import Forcing

def distance_2points(lat1,lon1,lat2,lon2) -> float:
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

def min_distance(lon, lat, lon_vec , lat_vec) -> Tuple[float, int]:
    """Calculates minimum distance between a given point and a list of
    point.

    Also returns index of the found minimum.
    """

    dx = []
    for n in range(len(lat_vec)):
        dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))

    return np.array(dx).min(), np.array(dx).argmin()

def lon_in_km(lat: float) -> float:
    """Converts one longitude degree to km for a given latitude."""

    return distance_2points(lat, 0, lat, 1)


def day_list(start_time, end_time):
    """Determins a Pandas data range of all the days in the time span of the InputModel objext"""
    days = pd.date_range(start=start_time.split('T')[0], end=end_time.split('T')[0], freq='D')
    return days

def month_list(start_time, end_time):
    """Determins a Pandas data range of all the months in the time span of the InputModel objext"""
    months = pd.date_range(start=start_time[:7], end=end_time[:7], freq='MS')
    return months

def add_extension(filename: str, extension: str):
    """Adds a file extension to the file name.

    If the file already has an extension, then no extension is added. An
    extension is defined as a . in the last five characters of the string.
    """

    if ('.' in filename[-5:]) or (not extension):
        return filename
    else:
        if not extension[0] == '.':
            extension = f".{extension}"
        return f"{filename}{extension}"


def add_prefix(filename: str, prefix: str) -> str:
    """Adds a prefix to a filename, e.g. FileName.txt -> new_FileName.txt"""
    if (not prefix == '') and (not prefix[-1] == '_'):
        return f"{prefix}_{filename}"
    else:
        return f"{prefix}{filename}"

def add_suffix(filename: str, suffix: str) -> str:
    """Adds a suffix to a filename, e.g. FileName.txt -> FileName_new.txt"""
    if (not suffix == '') and (not suffix[0] == '_'):
        suffix = f"_{suffix}"

    filename_list = filename.split('.')

    if len(filename_list) == 1:
        return filename+suffix
    else:
        return '.'.join(filename_list[0:-1]) + suffix + '.' + filename_list[-1]


def add_folder_to_filename(filename: str, folder: str) -> str:
    """Adds a folder to filename to create a path:

    e.g. FileName.txt, Folder -> Folder/FileName.txt
    """

    if (not folder == '') and (not folder[-1] == '/'):
        return f"{folder}/{filename}"
    else:
        return f"{folder}{filename}"

def create_filename_time(filestring: str, times: List, datestring: str='%Y%m%d%H%M') -> str:
    """Substitutes the strings #T0, #T1, #T2... etc. in filestring with time
    stamps from a list of times, using the format in datestring.

    e.g. #T0_file.txt, ['2020-05-04 18:00'], %Y%m%d%H%M -> 202005041800_file.txt
    """

    ct = 0
    for t in times:
        filestring = re.sub(f"#T{ct}", pd.Timestamp(t).strftime(datestring), filestring)
        ct = ct + 1

    return filestring

def create_filename_lonlat(filestring: str, lon: float, lat: float) -> str:
    """Substitutes the strings #Lon, #Lat in filestring with values of lon and
    lat.

    e.g. #Lon_#Lat_file.txt, 8.0, 60.05 -> 08.0000000_60.05000000_file.txt
    """

    filestring = re.sub("#Lon", f"{lon:010.7f}", filestring)
    filestring = re.sub("#Lat", f"{lat:010.7f}", filestring)

    return filestring

def create_filename_obj(filestring: str, objects: List[Union[Grid, Forcing, Boundary]]) -> str:
    """Substitutes the strings #{Object} in filestring with the name given to
    the object.

    e.g. #Grid_#Forcing_file.txt, [Grid(..., name="Sula"), Forcing(..., name='NORA3')]
        -> Sula_NORA3_file.txt
    """

    for object in objects:
        if object is not None:
            obj_str = type(object).__name__
            obj_name = object.name()
            filestring = re.sub(f"#{obj_str}", obj_name, filestring)

    return filestring

def clean_filename(filename: str, list_of_placeholders: List[str]) -> str:
    """ Cleans out the file name from possible used placeholders, e.g. #Grid
    as given in the list.

    Also removes multiple underscores '___'
    """

    for s in list_of_placeholders:
            filename = re.sub(s, '', filename)

    filename = re.sub("_{2,10}", '_', filename)

    return filename

def create_time_stamps(start_time: str, end_time: str, stride: int, hours_per_file: int=0, last_file: str='', lead_time: int=0) -> Tuple:
    """Create time stamps to read in blocks of wind forcing from files.

    Options

    start_time:     Wanted start times

    end_time:       Wanted end time

    stride:         Time between files (in hours). This many hours read from
                    each file.

    lead_time:      E.g. 12 means the time 12:00 is read from file 00:00, not
                    from time 12:00 (in hours; negative values accepted).

    last_file:      Don't try to read past a file with this time stamp.

    hours_per_file: Try to read this many hours from the last file. Only used
                    if last_file is given, and only meaningful if hours_per_file
                    is different from stride.

    Returns

    start_times:    Pandas DatetimeIndex with the start times.
    end_times:      Pandas DatetimeIndex with the end times.
    file_times:     Pandas DatetimeIndex with the time stamps for the files.

    I.e. loop through the objects and read data for the time
            start_times[n] - end_times[n]
    from a file with a time stamp
            file_times[n]
    """

    if hours_per_file == 0:
        hours_per_file = stride

    # FIND FILE STAMPS
    start_stamp = pd.Timestamp(start_time) - pd.DateOffset(hours=lead_time)
    if last_file is not '':
        end_stamp = pd.Timestamp(last_file)

        # E.g. we want to start a forecast at 06:00 but the last (and only) file is 00:00
        if start_stamp > end_stamp:
            start_stamp = end_stamp
    else:
        end_stamp = pd.Timestamp(end_time) - pd.DateOffset(hours=lead_time)

    # How many ours to remove if files are e.g. 00, 06, 12 and we request output from 01-08
    h0 = int(start_stamp.hour) % stride
    h1 = int(end_stamp.hour) % stride
    file_times = pd.date_range(start = start_stamp - pd.DateOffset(hours=h0), end = end_stamp - pd.DateOffset(hours=h1), freq=f'{stride}H')

    # FIND START AND END TIMES
    start_times = file_times + pd.DateOffset(hours=lead_time)
    end_times = start_times + pd.DateOffset(hours=stride-1)

    # First time might not coincide with first step in first file
    start_times.values[0] = pd.Timestamp(start_time)

    if last_file:
        # In operational systems we might want to read a longer segment from the last file
        end_times.values[-1] = min([pd.Timestamp(last_file) + pd.DateOffset(hours=(hours_per_file-1)), pd.Timestamp(end_time)])
    else:
        # Last time might not coincide with last step in last file
        end_times.values[-1] = pd.Timestamp(end_time)

    return start_times, end_times, file_times


def expand_area(lon_min: float, lon_max:float , lat_min:float , lat_max:float, expansion_factor: float) -> Tuple[float, float, float, float]:
    """
    Expands a lon-lat bounding box with an expansion factor.
    expansion_factor = 1 does nothing, and 1.5 expands 50% both directions.
    """

    expand_lon = (lon_max - lon_min)*(expansion_factor-1)*0.5
    expand_lat = (lat_max - lat_min)*(expansion_factor-1)*0.5

    new_lon_min = lon_min - expand_lon
    new_lon_max = lon_max + expand_lon

    new_lat_min = lat_min - expand_lat
    new_lat_max = lat_max + expand_lat

    return new_lon_min, new_lon_max, new_lat_min, new_lat_max



def check_if_folder(folder: str, create: bool=True) -> bool:
    """Creates a folder if it does not exist, and returns True if it
    already existed."""

    if folder == '':
        existed = True
    else:
        existed = os.path.isdir(folder)

    if not existed:
        os.mkdir(folder)

    return existed


# -----------------------------------------------------------------------------
# MISC STAND ALONE FUNCTIONS
# -----------------------------------------------------------------------------
def read_ww3_info(filename) -> Tuple[float, float, float, float, float, float, int, int]:
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


def u_v_from_dir(ws, wdir) -> Tuple[float, float]:
    """Converts wind speed and direction (from) to u and v components."""

    # see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))

    return u, v

def interp_spec(f, D, S, fi, Di):
    """Interpolates a spectrum to new frequncy and directions.

    Spectrum is two dimensional (len(f), len(D))

    """
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

def flip_spec(spec, D):
    """Flips the directionality of the spectrum (clock/anticlockwise).

    To flip the directional vector, use flip_spec(D,D)
    """

    # This check enables us to flip directions with flip_spec(D,D)
    if len(spec.shape) == 1:
        flipping_dir = True
        spec = np.array([spec])
    else:
        flipping_dir = False
    spec_flip = np.zeros(spec.shape)

    ind = np.arange(0,len(D), dtype='int')
    #dD = np.diff(D).mean()
    dD = 360/len(D)
    steps = D/dD # How many delta-D from 0

    # Need to move indeces the other way if the vector is decreasing than if it is increasing
    direction = np.sign(np.median(np.diff(D)))
    ind_flip = ((ind - 2*steps*direction).astype(int) + len(D)) % len(D)

    spec_flip=spec[..., list(ind_flip)]

    if flipping_dir:
        spec_flip = spec_flip[0]

    return spec_flip

def shift_spec(spec, D, shift=0):
    """Shifts the spectrum D degree. To shift the directional vector, use

    shift_spec(D, D, shift)
    """

    # This check enables us to flip directions with shift_spec(D, D, shift)
    if len(spec.shape) == 1:
        shifting_dir = True
        spec = np.array([spec])
    else:
        shifting_dir = False
    spec_shift = np.zeros(spec.shape)

    ind = np.arange(0,len(D), dtype='int')
    dD = 360/len(D)

    if abs(np.floor(shift/dD)-shift/dD) > 0:
        raise Exception (f'Shift {shift} needs to be multiple of frequency resolution {dD}, but shift/dD={shift/dD}! Otherwise interpolation would be needed.')

    ind_flip = ((ind + int(shift/dD)).astype(int) + len(D)) % len(D)

    spec_shift=spec[..., list(ind_flip)]
    if shifting_dir:
        spec_shift = spec_shift[0]

    return spec_shift
