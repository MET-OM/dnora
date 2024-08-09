import xarray as xr
import numpy as np
import pandas as pd
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.aux_funcs import lon_in_km
from subprocess import call


def create_fimex_xy_strings(
    lon: tuple[float, float], lat: tuple[float, float], resolution_in_km: float
) -> tuple[str, str]:
    # Set resolution
    dlat = resolution_in_km / 111
    mean_lon_in_km = (lon_in_km(lat[0]) + lon_in_km(lat[-1])) * 0.5
    dlon = resolution_in_km / mean_lon_in_km

    if len(np.unique(lon)) == 1:
        x_str = str(lon[0])
    else:
        if lon[1] < lon[0] + dlon:
            lon = (lon[0] - dlon, lon[0] + dlon)
        x_str = str(lon[0]) + "," + str(lon[0] + dlon) + ",...," + str(lon[-1])
    if len(np.unique(lat)) == 1:
        y_str = str(lat[0])
    else:
        if lat[1] < lat[0] + dlat:
            lat = (lat[0] - dlat, lat[0] + dlat)
        y_str = str(lat[0]) + "," + str(lat[0] + dlat) + ",...," + str(lat[-1])

    return x_str, y_str


def create_pyfimex_xy_vectors(
    lon: tuple[float, float], lat: tuple[float, float], resolution_in_km: float
) -> tuple[np.ndarray, np.ndarray]:
    # Set resolution
    dlat = resolution_in_km / 111
    mean_lon_in_km = (lon_in_km(lat[0]) + lon_in_km(lat[-1])) * 0.5
    dlon = resolution_in_km / mean_lon_in_km

    lon_vec = np.arange(lon[0], lon[1] + dlon, dlon)
    lat_vec = np.arange(lat[0], lat[1] + dlat, dlat)
    return lon_vec, lat_vec


def create_fimex_command(
    nc_fimex: str,
    url: str,
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    lon: np.ndarray,
    lat: np.ndarray,
    resolution_in_km: float,
    variables: list[str],
) -> list[str]:
    x_str, y_str = create_fimex_xy_strings(lon, lat, resolution_in_km=resolution_in_km)
    fimex_command = [
        "fimex",
        "--input.file=" + url,
        "--interpolate.method=bilinear",
        "--interpolate.projString=+proj=latlong +ellps=sphere +a=6371000 +e=0",
        "--interpolate.xAxisValues=" + x_str + "",
        "--interpolate.yAxisValues=" + y_str + "",
        "--interpolate.xAxisUnit=degree",
        "--interpolate.yAxisUnit=degree",
        "--process.rotateVector.all",
    ]
    for var in variables:
        fimex_command.append("--extract.selectVariables=" + var)

    fimex_command += [
        "--extract.reduceTime.start=" + start_time.strftime("%Y-%m-%dT%H:%M:%S"),
        "--extract.reduceTime.end=" + end_time.strftime("%Y-%m-%dT%H:%M:%S"),
        "--process.rotateVector.direction=latlon",
        "--output.file=" + nc_fimex,
    ]

    return fimex_command


def ds_fimex_read(
    start_time,
    end_time,
    url,
    lon,
    lat,
    resolution_in_km: float,
    variables: list[str],
    data_type: DnoraDataType,
    name: str,
    program: str = "pyfimex",  #'fimex' or 'pyfimex'
) -> xr.Dataset:
    nc_fimex = f"dnora_{data_type.name.lower()}_temp/{name}_{start_time.strftime('%Y%m%d%H%M')}_fimex.nc"
    if program not in ["fimex", "pyfimex"]:
        raise ValueError("program need to be 'fimex' or 'pyfimex'!")

    if program == "pyfimex":
        x_vec, y_vec = create_pyfimex_xy_vectors(lon, lat, resolution_in_km)
        pyfimex(
            input_file=url,
            output_file=nc_fimex,
            projString="+proj=latlong +ellps=sphere +a=6371000 +e=0",
            xAxisValues=x_vec,
            yAxisValues=y_vec,
            selectVariables=variables,
            reduceTime_start=start_time.strftime("%Y-%m-%dT%H:%M:%S"),
            reduceTime_end=end_time.strftime("%Y-%m-%dT%H:%M:%S"),
        )
    elif program == "fimex":
        fimex_command = create_fimex_command(
            nc_fimex,
            url,
            start_time,
            end_time,
            lon,
            lat,
            resolution_in_km=resolution_in_km,
            variables=variables,
        )
        call(fimex_command)
    ds = xr.open_dataset(nc_fimex).squeeze()
    return ds


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
