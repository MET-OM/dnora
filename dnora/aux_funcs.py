import numpy as np
import pandas as pd
from scipy import interpolate
import geopy.distance
from pathlib import Path
import calendar
import re
import os


def distance_2points(lat1, lon1, lat2, lon2) -> float:
    """Calculate distance between two points"""
    return geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).km


def min_distance(lon, lat, lon_vec, lat_vec) -> tuple[float, int]:
    """Calculates minimum distance [km] between a given point and a list of
    points given in spherical coordinates (lon/lat degrees).

    Also returns index of the found minimum.
    """
    dx = []
    for n, __ in enumerate(lat_vec):
        dx.append(distance_2points(lat, lon, lat_vec[n], lon_vec[n]))

    return np.array(dx).min(), np.array(dx).argmin()


def min_cartesian_distance(x, y, x_vec, y_vec) -> tuple[float, int]:
    """ "Calculates minimum distance [km] between a given point and list of points given
    in cartesian coordinates [m].

    Also returns incex of found minimum"""
    dx = ((y - y_vec) ** 2 + (x - x_vec) ** 2) ** 0.5

    return dx.min() / 1000, dx.argmin()


def lon_in_km(lat: float) -> float:
    """Converts one longitude degree to km for a given latitude."""

    return distance_2points(lat, 0, lat, 1)


def domain_size_in_km(
    lon: tuple[float, float], lat: tuple[float, float]
) -> tuple[float, float]:
    """Calculates approximate size of grid in km."""

    km_x = distance_2points(
        (lat[0] + lat[1]) / 2, lon[0], (lat[0] + lat[1]) / 2, lon[1]
    )
    km_y = distance_2points(lat[0], lon[0], lat[1], lon[0])

    return km_x, km_y


def get_coordinates_from_ds(ds, return_dict: bool = False) -> tuple:
    """Determins if an xarray dataset is cartesian (x,y) or spherical (lon,lat)
    and returns the vecotrs (None for the ones that are not defined).

    If lon, lat is defined over time, the firs instance is grabbed."""

    if "time" in ds.dims:
        ds = ds.isel(time=0)

    if hasattr(ds, "lon") and hasattr(ds, "lat"):
        lon, lat = np.squeeze(ds.lon.values), np.squeeze(ds.lat.values)
        x, y = None, None

    if hasattr(ds, "longitude") and hasattr(ds, "latitude"):
        lon, lat = np.squeeze(ds.longitude.values), np.squeeze(ds.latitude.values)
        x, y = None, None
    if hasattr(ds, "x") and hasattr(ds, "y"):
        x, y = np.squeeze(ds.x.values), np.squeeze(ds.y.values)
        lon, lat = None, None

    if all_none([x, y, lon, lat]):
        raise AttributeError(
            "Dataset doesn't have a combination of lon(gitude)/lat(itude) or x/y!"
        )

    if return_dict:
        return {"lon": lon, "lat": lat, "x": x, "y": y}
    else:
        return lon, lat, x, y


def all_none(val) -> bool:
    return not [a for a in val if a is not None]


def is_gridded(data: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> bool:
    if lon is None or lat is None:
        return False

    if data.shape == (len(lat), len(lon)):
        return True

    if len(data.shape) == 1 and len(lat) == data.shape[0] and len(lon) == data.shape[0]:
        return False

    raise Exception(
        f"Size of data is {data.shape} but len(lat) = {len(lat)} and len(lon) = {len(lon)}. I don't know what is going on!"
    )


#
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
    h0 = int(start_stamp.hour) % stride + offset
    h1 = int(end_stamp.hour) % stride + offset
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


def expand_area(
    lon: tuple[float, float], lat: tuple[float, float], expansion_factor: float
) -> tuple[float, float, float, float]:
    """
    Expands a lon-lat bounding box with an expansion factor.
    expansion_factor = 1 does nothing, and 1.5 expands 50% both directions.
    """

    expand_lon = (lon[1] - lon[0]) * (expansion_factor - 1) * 0.5
    expand_lat = (lat[1] - lat[0]) * (expansion_factor - 1) * 0.5

    new_lon = lon[0] - expand_lon, lon[1] + expand_lon
    new_lat = lat[0] - expand_lat, lat[1] + expand_lat

    return new_lon, new_lat


def check_if_file(filename: str, halt=False) -> bool:
    """Checks if a file exists and halts with an error if halt=True"""
    exists = os.path.isfile(filename)
    if not exists and halt:
        raise FileNotFoundError(f"DNORA cannot find file {filename}")
    return exists


def check_if_folder(folder: str, create: bool = True) -> bool:
    """Creates a folder if it does not exist, and returns True if it
    already existed."""

    if folder == "":
        existed = True
    else:
        existed = os.path.isdir(folder)

    if not existed:
        os.mkdir(folder)

    return existed


# -----------------------------------------------------------------------------
# MISC STAND ALONE FUNCTIONS
# -----------------------------------------------------------------------------
def read_ww3_info(
    filename,
) -> tuple[float, float, float, float, float, float, int, int]:
    """Read grid specification from the GridName_info.txt file"""
    with open(filename, "r") as f:
        lines = f.readlines()

    for n in range(len(lines)):
        line = lines[n].split()

        if len(line):
            if line[0] == "lon:":
                lon_min = float(line[1])
                lon_max = float(line[3][0:-1])
                lat_min = float(line[5])
                lat_max = float(line[7])
            elif line[0] == "dlon,":
                dlon = float(line[3][0:-1])
                dlat = float(line[4])
            elif line[0] == "nx,":
                nx = int(line[3])
                ny = int(line[5])
    return lon_min, lon_max, lat_min, lat_max, dlon, dlat, nx, ny


def u_v_from_speed_dir(ws, wdir) -> tuple[float, float]:
    """Converts wind speed and direction (from) to u and v components."""

    # see http://tornado.sfsu.edu/geosciences/classes/m430/Wind/WindDirection.html
    u = -ws * (np.sin(np.deg2rad(wdir)))
    v = -ws * (np.cos(np.deg2rad(wdir)))

    return u, v


def speed_dir_from_u_v(u, v) -> tuple[float, float]:
    """Convert component to speed and direction (from)"""
    ws = (u**2 + v**2) ** 0.5
    wdir = np.mod((90 - np.rad2deg(np.arctan2(v, u))) + 180, 360)
    return ws, wdir


def interp_spec(f, D, S, fi, Di):
    """Interpolates a spectrum to new frequncy and directions.
    Spectrum is two dimensional (len(f), len(D))
    """
    Sleft = S
    Sright = S
    Dleft = -D[::-1]
    Dright = D + 360

    bigS = np.concatenate((Sleft, S, Sright), axis=1)
    bigD = np.concatenate((Dleft, D, Dright))

    Finterpolator = interpolate.RectBivariateSpline(f, bigD, bigS, kx=1, ky=1, s=0)
    Si = Finterpolator(fi, Di)

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

    ind = np.arange(0, len(D), dtype="int")
    # dD = np.diff(D).mean()
    dD = 360 / len(D)
    steps = D / dD  # How many delta-D from 0

    # Need to move indeces the other way if the vector is decreasing than if it is increasing
    direction = np.sign(np.median(np.diff(D)))
    ind_flip = ((ind - 2 * steps * direction).astype(int) + len(D)) % len(D)

    spec_flip = spec[..., list(ind_flip)]

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

    ind = np.arange(0, len(D), dtype="int")
    dD = 360 / len(D)

    if abs(np.floor(shift / dD) - shift / dD) > 0:
        raise Exception(
            f"Shift {shift} needs to be multiple of frequency resolution {dD}, but shift/dD={shift/dD}! Otherwise interpolation would be needed."
        )

    ind_flip = ((ind + int(shift / dD)).astype(int) + len(D)) % len(D)

    spec_shift = spec[..., list(ind_flip)]
    if shifting_dir:
        spec_shift = spec_shift[0]

    return spec_shift


def check_that_spectra_are_consistent(
    spec, dirs, freq, expected_dim: int = None
) -> int:
    if spec.shape[-1] == len(dirs) and spec.shape[-2] == len(freq):
        spec_dim = 2
    elif spec.shape[-1] == len(freq) and spec.shape[-1] == len(dirs):
        spec_dim = 1
    else:
        spec_dim = -1

    if expected_dim is None:
        return spec_dim

    if spec_dim != expected_dim:
        if spec_dim not in [1, 2]:
            ValueError("Provided array does not contain valid 1D or 2D spectra!")
        else:
            ValueError(
                f"Expected {expected_dim} dimensional spectra, but they seem to be {spec_dim} dimensional!"
            )


def identify_boundary_edges(boundary_mask: np.ndarray) -> list[str]:
    """Identifies which edges has some boundary points
    North = [-1,:]
    South = [0,:]
    East = [:,-1]
    West = [:,0]
    """
    edges = []

    if boundary_mask.shape[1] > 2:
        n0 = 1
        n1 = -1
    else:
        n0 = 0
        n1 = None

    if np.any(boundary_mask[-1, n0:n1]):
        edges.append("N")

    if np.any(boundary_mask[0, n0:n1]):
        edges.append("S")

    if boundary_mask.shape[0] > 2:
        n0 = 1
        n1 = -1
    else:
        n0 = 0
        n1 = None

    if np.any(boundary_mask[n0:n1, -1]):
        edges.append("E")

    if np.any(boundary_mask[n0:n1, 0]):
        edges.append("W")

    return edges


def create_ordered_boundary_list(edge_list):
    """Gets all edges, but given in a continuous clockwise direction.
    If this is not possible (e.g. ['N','S']), then ampty list is returned."""
    full_list = ["N", "E", "S", "W"]
    for ind, edge in enumerate(full_list):
        if edge not in edge_list:
            full_list[ind] = ""
    full_array = np.array(full_list)

    if len(np.where(full_array == "")[0]) == 0:
        return full_array.tolist()

    ct = 0
    while (
        np.where(full_array == "")[0][-1] != len(np.where(full_array == "")[0]) - 1
    ) and ct < 5:
        full_array = np.roll(full_array, 1)
        ct += 1

    if ct > 4:
        print(
            f"No continuous boundary can be found for edges {edge_list}. Returning empy list."
        )
        return []

    full_array = full_array[full_array != ""]

    return full_array.tolist()


def get_coords_for_boundary_edges(
    edges: list, lon_edges: tuple[float, float], lat_edges: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """Create coordinate vectors for clockwise running edges.
    Assumes that edges are clockwise and continuous, which is imposed by the
    function create_ordered_boundary_list.
    Empty list return empty arrays.
    """
    lon = []
    lat = []
    for edge in edges:
        if edge == "N":
            lon.append(lon_edges[0])
            lat.append(lat_edges[1])
        if edge == "S":
            lon.append(lon_edges[1])
            lat.append(lat_edges[0])
        if edge == "W":
            lon.append(lon_edges[0])
            lat.append(lat_edges[0])
        if edge == "E":
            lon.append(lon_edges[1])
            lat.append(lat_edges[1])

    if edges:
        edge = edges[-1]  # Close the loop
        if edge == "N":
            lon.append(lon_edges[1])
            lat.append(lat_edges[1])
        if edge == "S":
            lon.append(lon_edges[0])
            lat.append(lat_edges[0])
        if edge == "W":
            lon.append(lon_edges[0])
            lat.append(lat_edges[1])
        if edge == "E":
            lon.append(lon_edges[1])
            lat.append(lat_edges[0])

    return np.array(lon), np.array(lat)


def create_swan_segment_coords(boundary_mask, lon_edges, lat_edges):
    """Createsa longitude and latitude arrays for the SWAN BOUND SEGEMENT
    command based on boundary mask.
    Identifies edges (north, south etc.) and sets all the points on the edges
    as boundary points.
    If no continuous boundary can be identified, it returns empty list.
    """

    edge_list = identify_boundary_edges(boundary_mask)
    clean_edge_list = create_ordered_boundary_list(edge_list)
    lon, lat = get_coords_for_boundary_edges(clean_edge_list, lon_edges, lat_edges)
    return lon, lat


def pyfimex(
    input_file,
    output_file,
    projString,
    xAxisValues,
    yAxisValues,
    selectVariables,
    reduceTime_start,
    reduceTime_end,
    ensemble_member=False,
):
    import pyfimex0 as pyfi

    r = pyfi.createFileReader("netcdf", input_file)
    inter_ll = pyfi.createInterpolator(r)
    inter_ll.changeProjection(
        pyfi.InterpolationMethod.BILINEAR,
        projString,
        xAxisValues,
        yAxisValues,
        "degree",
        "degree",
    )
    extra = pyfi.createExtractor(inter_ll)
    extra.selectVariables(selectVariables)
    extra.reduceTimeStartEnd(reduceTime_start, reduceTime_end)
    if ensemble_member == True:
        extra.reduceDimensionStartEnd("ensemble_member", 1, 1)
    pyfi.createFileWriter(extra, "netcdf", output_file)


def get_url(
    folder: str,
    filename: str,
    time_stamp: pd.DatetimeIndex | str = None,
    get_list: bool = False,
) -> str:
    if not isinstance(filename, list):
        filename = [filename]
    if time_stamp is not None:
        time_stamp = pd.to_datetime(time_stamp)
        filename = [time_stamp.strftime(fn) for fn in filename]
        folder = time_stamp.strftime(folder)
    url = []
    for fn in filename:
        url_temp = Path(folder).joinpath(fn)
        url_temp = re.sub(f"https:/", "https://", str(url_temp), 1)
        url.append(url_temp)
    if len(url) == 1 and not get_list:
        return os.path.expanduser(url[0])
    else:
        return [os.path.expanduser(u) for u in url]
