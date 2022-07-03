import numpy as np
import pandas as pd
def is_gridded(data: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> bool:
    if data.shape == (len(lat), len(lon)):
        return True

    if len(data.shape) == 1 and len(lat) == data.shape[0] and len(lon) == data.shape[0]:
        return False

    raise Exception(f"Size of data is {data.shape} but len(lat) = {len(lat)} and len(lon) = {len(lon)}. I don't know what is going on!")

def day_list(start_time, end_time):
    """Determins a Pandas data range of all the days in the time span"""
    t0 = pd.Timestamp(start_time).strftime('%Y-%m-%d')
    t1 = pd.Timestamp(end_time).strftime('%Y-%m-%d')
    days = pd.date_range(start=t0, end=t1, freq='D')
    return days.strftime('%Y-%m-%d').to_list()

def month_list(start_time, end_time):
    """Determins a Pandas data range of all the months in the time span"""
    t0 = pd.Timestamp(start_time).strftime('%Y-%m')
    t1 = pd.Timestamp(end_time).strftime('%Y-%m')
    months = pd.date_range(start=t0, end=t1, freq='MS')
    return months.strftime('%Y-%m').to_list()

def year_list(start_time, end_time):
    """Determins a Pandas data range of all the years in the time span"""
    t0 = pd.Timestamp(start_time).strftime('%Y-%m-%d')
    t1 = pd.Timestamp(end_time).strftime('%Y-%m-%d')
    years = pd.date_range(start=t0, end=t1, freq='YS')
    return years.strftime('%Y').to_list()
