import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from scipy import interpolate
from statistics import mode
import os, re

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

def min_distance(lon, lat, lon_vec , lat_vec):
    dx = []
    for n in range(len(lat_vec)):
        dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))

    return np.array(dx).min(), np.array(dx).argmin()

def lon_in_km(lat: float) -> float:
    return distance_2points(lat, 0, lat, 1)


def day_list(start_time, end_time):
    """Determins a Pandas data range of all the days in the time span of the InputModel objext"""
    days = pd.date_range(start=start_time.split('T')[0], end=end_time.split('T')[0], freq='D')
    return days

def month_list(start_time, end_time):
    """Determins a Pandas data range of all the months in the time span of the InputModel objext"""
    months = pd.date_range(start=start_time[:7], end=end_time[:7], freq='MS')
    return months

def add_file_extension(filename: str, extension: str):
    if not extension[0] == '.':
        extension = f".{extension}"

    if not filename[-(len(extension)):] == extension:
        return f"{filename}{extension}"
    else:
        return filename

def add_folder_to_filename(filename: str, folder: str):
    if (not folder == '') and (not folder[-1] == '/'):
        return f"{folder}/{filename}"
    else:
        return f"{folder}{filename}"

def create_filename_time(filestring: str, times, datestring: str='%Y%m%d%H%M'):
    ct = 0
    for t in times:
        filestring = re.sub(f"\$T{ct}", pd.Timestamp(t).strftime(datestring), filestring)
        ct = ct + 1

    return filestring

def create_filename_lonlat(filestring: str, lon: float, lat: float):
    filestring = re.sub(f"\$Lon", str(lon), filestring)
    filestring = re.sub(f"\$Lat", str(lat), filestring)

    return filestring

def create_filename_obj(filestring: str, objects):
    for object in objects:
        obj_str = type(object).__name__
        obj_name = object.name()
        filestring = re.sub(f"\${obj_str}", obj_name, filestring)

    return filestring

def create_time_stamps(start_time: str, end_time: str, stride: int, hours_per_file: int = 0, last_file: str = '', lead_time: int = 0):
    """Create time stamps to read in blocks of wind forcing from files"""
    if hours_per_file == 0:
        hours_per_file = stride

    # FIND FILE STAMPS
    #t0 = np.datetime64(start_time) - np.timedelta64(lead_time,'h')
    start_stamp = pd.Timestamp(start_time) - pd.DateOffset(hours=lead_time)
    if last_file is not '':
        end_stamp = pd.Timestamp(last_file)

        # E.g. we want to start a forecast at 06:00 but the last (and only) file is 00:00
        if start_stamp > end_stamp:
            start_stamp = end_stamp
    else:
        #t1 = np.datetime64(end_time) - np.timedelta64(lead_time,'h')
        end_stamp = pd.Timestamp(end_time) - pd.DateOffset(hours=lead_time)



    # How many ours to remove if files are e.g. 00, 06, 12 and we request output from 01-08
    h0 = int(start_stamp.hour) % stride
    h1 = int(end_stamp.hour) % stride
    file_times = pd.date_range(start = start_stamp - pd.DateOffset(hours=h0), end = end_stamp - pd.DateOffset(hours=h1), freq=f'{stride}H')


    # FIND START AND END TIMES
    start_times = file_times + pd.DateOffset(hours=lead_time)
    end_times = start_times + pd.DateOffset(hours=stride-1)

    # First time might not coincide with first step in first file
    #start_times.values[0] = np.datetime64(start_time)
    start_times.values[0] = pd.Timestamp(start_time)

    if last_file is not '':
        # In operational systems we might want to read a longer segment from the last file
        end_times.values[-1] = min([pd.Timestamp(last_file) + pd.DateOffset(hours=(hours_per_file-1)), pd.Timestamp(end_time)])
    else:
        # Last time might not coincide with last step in last file
        end_times.values[-1] = pd.Timestamp(end_time)
    return start_times, end_times, file_times


def expand_area(lon_min: float, lon_max:float , lat_min:float , lat_max:float, expansion_factor: float):
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


def u_v_from_dir(ws, wdir):
    # see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))
    return u, v

    return grid

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
