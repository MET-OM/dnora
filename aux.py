import numpy as np
import pandas as pd
from scipy.interpolate import griddata

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

def day_list(start_time, end_time):
    """Determins a Pandas data range of all the days in the time span of the InputModel objext"""
    days = pd.date_range(start=start_time.split('T')[0], end=end_time.split('T')[0], freq='D')
    return days

def month_list(start_time, end_time):
    """Determins a Pandas data range of all the months in the time span of the InputModel objext"""
    months = pd.date_range(start=start_time[:7], end=end_time[:7], freq='MS')
    return months

def create_time_stamps(start_time, end_time, stride, hours_per_file = None, last_file = None, lead_time = 0):
    """Create time stamps to read in blocks of wind forcing from files"""
    if hours_per_file is None:
        hours_per_file = stride

    # FIND FILE STAMPS
    t0 = np.datetime64(start_time) - np.timedelta64(lead_time,'h')
    if last_file is not None:
        t1 = np.datetime64(last_file)

        # E.g. we want to start a forecast at 06:00 but the last (and only) file is 00:00
        if t0 > t1:
            t0 = t1
    else:
        t1 = np.datetime64(end_time) - np.timedelta64(lead_time,'h')



    # How many ours to remove if files are e.g. 00, 06, 12 and we request output from 01-08
    h0 = int(pd.Timestamp(t0).hour) % stride
    h1 = int(pd.Timestamp(t1).hour) % stride
    file_times = pd.date_range(start = t0 - np.timedelta64(h0,'h'), end = t1 - np.timedelta64(h1,'h'), freq=f'{stride}H')


    # FIND START AND END TIMES
    start_times = file_times + np.timedelta64(lead_time,'h')
    end_times = start_times + np.timedelta64(stride-1, 'h')

    # First time might not coincide with first step in first file
    start_times.values[0] = np.datetime64(start_time)

    if last_file is not None and hours_per_file is not None:
        # In operational systems we might want to read a longer segment from the last file
        end_times.values[-1] = min([np.datetime64(last_file) + np.timedelta64((hours_per_file-1),'h'), np.datetime64(end_time)])
    else:
        # Last time might not coincide with last step in last file
        end_times.values[-1] = np.datetime64(end_time)
    return start_times, end_times, file_times



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
# -----------------------------------------------------------------------------

