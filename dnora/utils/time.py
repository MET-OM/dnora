import pandas as pd
import calendar
import numpy as np


def day_list(start_time, end_time):
    """Determins a Pandas data range of all the days in the time span"""
    t0 = pd.Timestamp(start_time).strftime("%Y-%m-%d")
    t1 = pd.Timestamp(end_time).strftime("%Y-%m-%d")
    days = pd.date_range(start=t0, end=t1, freq="D")
    return days.strftime("%Y-%m-%d").to_list()


#
def month_list(start_time, end_time, fmt="%Y-%m"):
    """Determins a Pandas data range of all the months in the time span"""
    t0 = pd.Timestamp(start_time).strftime("%Y-%m")
    t1 = pd.Timestamp(end_time).strftime("%Y-%m")
    months = pd.date_range(start=t0, end=t1, freq="MS")
    return months.strftime(fmt).to_list()


#
def year_list(start_time, end_time):
    """Determins a Pandas data range of all the years in the time span"""
    t0 = pd.Timestamp(start_time).strftime("%Y-%m-%d")
    t1 = pd.Timestamp(end_time).strftime("%Y-%m-%d")
    years = pd.date_range(start=t0, end=t1, freq="YS")
    return years.strftime("%Y").to_list()


def int_list_of_years(start_time, end_time):
    year0 = min(pd.Series(pd.to_datetime(day_list(start_time, end_time))).dt.year)
    year1 = max(pd.Series(pd.to_datetime(day_list(start_time, end_time))).dt.year)
    return np.linspace(year0, year1, year1 - year0 + 1).astype(int)


def int_list_of_months(start_time, end_time):
    if len(int_list_of_years(start_time, end_time)) > 1:
        raise Exception("Only use this function for times within a single year!")
    month0 = min(pd.to_datetime(pd.Series(day_list(start_time, end_time))).dt.month)
    month1 = max(pd.to_datetime(pd.Series(day_list(start_time, end_time))).dt.month)
    return np.linspace(month0, month1, month1 - month0 + 1).astype(int)


def int_list_of_days(start_time, end_time):
    if len(int_list_of_months(start_time, end_time)) > 1:
        raise Exception("Only use this function for times within a single month!")
    day0 = min(pd.to_datetime(pd.Series(day_list(start_time, end_time))).dt.day)
    day1 = max(pd.to_datetime(pd.Series(day_list(start_time, end_time))).dt.day)
    return np.linspace(day0, day1, day1 - day0 + 1).astype(int)


def create_monthly_stamps(start_time: str, end_time: str) -> tuple:
    months = month_list(start_time, end_time)
    start_times = []
    end_times = []
    for month in months:
        n_of_days = calendar.monthrange(
            pd.to_datetime(month).year, pd.to_datetime(month).month
        )[1]

        start_times.append(pd.to_datetime(month))
        end_times.append(
            pd.to_datetime(month) + pd.Timedelta(hours=(n_of_days * 24 - 1))
        )
    start_times[0] = np.max([start_times[0], pd.to_datetime(start_time)])
    end_times[-1] = np.min([end_times[-1], pd.to_datetime(end_time)])
    return pd.to_datetime(start_times), pd.to_datetime(end_times)

def create_yearly_stamps(start_time: str, end_time: str) -> tuple:
    years = year_list(start_time, end_time)
    start_times = []
    end_times = []
    for year in years:
        t0 = pd.to_datetime(f"{year}-01-01 00:00:00")
        t1 = pd.to_datetime(f"{year}-12-31 23:59:59")
        start_times.append(t0)
        end_times.append(t1)
    start_times[0] = np.max([start_times[0], pd.to_datetime(start_time)])
    end_times[-1] = np.min([end_times[-1], pd.to_datetime(end_time)])
    return pd.to_datetime(start_times), pd.to_datetime(end_times)

def create_time_stamps(
    start_time: str,
    end_time: str,
    stride: int,
    hours_per_file: int = 0,
    last_file: str = "",
    lead_time: int = 0,
    offset: int = 0,
) -> tuple:
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
    offset:         if runs are e.g. 06 and 12, utse offset=6, stride=12
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
    if last_file != "":
        end_stamp = pd.Timestamp(last_file)

        # E.g. we want to start a forecast at 06:00 but the last (and only) file is 00:00
        if start_stamp > end_stamp:
            start_stamp = end_stamp
    else:
        end_stamp = pd.Timestamp(end_time) - pd.DateOffset(hours=lead_time)

    # How many ours to remove if files are e.g. 00, 06, 12 and we request output from 01-08
    h0 = int(start_stamp.hour - offset) % stride
    h1 = int(end_stamp.hour - offset) % stride
    file_times = pd.date_range(
        start=start_stamp - pd.DateOffset(hours=h0),
        end=end_stamp - pd.DateOffset(hours=h1),
        freq=f"{stride}h",
    )

    # FIND START AND END TIMES
    start_times = file_times + pd.DateOffset(hours=lead_time)
    end_times = start_times + pd.DateOffset(hours=stride - 1)

    # First time might not coincide with first step in first file
    start_times.values[0] = pd.Timestamp(start_time)

    if last_file:
        # In operational systems we might want to read a longer segment from the last file
        end_times.values[-1] = min(
            [
                pd.Timestamp(last_file) + pd.DateOffset(hours=(hours_per_file - 1)),
                pd.Timestamp(end_time),
            ]
        )
    else:
        # Last time might not coincide with last step in last file
        end_times.values[-1] = pd.Timestamp(end_time)

    return start_times, end_times, file_times

def get_first_file(
    start_time: str,
    stride: int,
    lead_time: int = 0,
    offset: int = 0,
):

    _, _, file_times = create_time_stamps(
        start_time, start_time, stride, lead_time=lead_time, offset=offset
    )
    return file_times[0]