import os
import xarray as xr
import glob
import numpy as np
from typing import Tuple
import pandas as pd
from .process import RemoveEmpty

# Import abstract classes and needed instances of them
from dnora.read.abstract_readers import SpectralDataReader
from dnora.grid import Grid
import cdsapi
from dnora.type_manager.data_sources import DataSource

# Import aux_funcsiliry functions
from dnora import msg
from dnora.aux_funcs import (
    int_list_of_days,
    int_list_of_months,
    int_list_of_years,
)
from dnora.spectral_conventions import SpectralConvention


def renormalize_era5_spec(bnd_spec):
    bnd_spec = bnd_spec.assign_coords(direction=np.arange(7.5, 352.5 + 15, 15))
    bnd_spec = bnd_spec.assign_coords(
        frequency=np.full(30, 0.03453) * (1.1 ** np.arange(0, 30))
    )
    bnd_spec = 10**bnd_spec
    bnd_spec = bnd_spec.fillna(0)
    return bnd_spec


def reshape_bnd_spec(bnd_spec):
    pass


def download_era5_from_cds(
    start_time, end_time, lon, lat, dlon, dlat, folder="dnora_wnd_temp"
) -> str:
    """Downloads ERA5 wave spectra from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    c = cdsapi.Client()

    filename = f"{folder}/EC_ERA5.nc"

    years = [f"{y:4.0f}" for y in int_list_of_years(start_time, end_time)]
    months = [f"{m:02.0f}" for m in int_list_of_months(start_time, end_time)]
    days = [f"{d:02.0f}" for d in int_list_of_days(start_time, end_time)]

    # Create string for dates
    dates = []
    for y in years:
        for m in months:
            for d in days:
                dates.append(f"{y}-{m}-{d}")
    dates = "/".join(dates)

    cds_command = {
        "class": "ea",
        "date": dates,
        "direction": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24",
        "domain": "g",
        "expver": "1",
        "frequency": "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30",
        "param": "251.140",
        "stream": "wave",
        "time": "00:00:00/03:00:00/06:00:00/09:00:00/12:00:00/15:00:00/18:00:00/21:00:00",
        "area": f"{lat[1]}/{lon[0]}/{lat[0]}/{lon[1]}",  # north, west, south, east
        "grid": f"{dlat}/{dlon}",
        "type": "an",
        "format": "netcdf",
    }

    c.retrieve("reanalysis-era5-complete", cds_command, filename)
    return filename


class ERA5(SpectralDataReader):
    def __init__(self, dlon: float = None, dlat: float = None):
        """Default dlon=dlat=0.5. If only dlon given, dlon=dlat and vice versa."""
        if dlon is None and dlat is None:
            dlon = 0.5
            dlat = 0.5

        if dlon is None:
            dlon = dlat

        if dlat is None:
            dlat = dlon

        self.dlon = dlon
        self.dlat = dlat

    def convention(self) -> str:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def post_processing(self):
        return RemoveEmpty()

    def get_coordinates(
        self, grid, start_time, source: DataSource, folder: str, **kwargs
    ) -> dict:
        """Reads first time instance of first file to get longitudes and latitudes for the PointPicker"""
        # point_list = self.get_restricted_area()._point_list()
        # lon_all = point_list[:,0]
        # lat_all = point_list[:,1]

        # return lon_all, lat_all
        lon = np.floor(np.array(grid.edges("lon")) / self.dlon) * self.dlon
        lat = np.floor(np.array(grid.edges("lat")) / self.dlat) * self.dlat

        self._given_grid = Grid(
            lon=(lon[0] - self.dlon, lon[-1] + self.dlon),
            lat=(lat[0] - self.dlat, lat[-1] + self.dlat),
        )
        self._given_grid.set_spacing(dlon=self.dlon, dlat=self.dlat)
        lon_all, lat_all = self._given_grid.lonlat()

        return {"lat": lat_all, "lon": lon_all}

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        inds,
        **kwargs,
    ) -> tuple[dict]:
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        msg.info(f"Getting ERA5 boundary spectra from {start_time} to {end_time}")

        temp_folder = "dnora_bnd_temp"
        if not os.path.isdir(temp_folder):
            os.mkdir(temp_folder)
            print("Creating folder %s..." % temp_folder)

        msg.plain("Removing old files from temporary folder...")
        for f in glob.glob(f"{temp_folder}/EC_ERA5.nc"):
            os.remove(f)

        nc_file = download_era5_from_cds(
            start_time,
            end_time,
            lon=self._given_grid.edges("lon"),
            lat=self._given_grid.edges("lat"),
            dlon=self.dlon,
            dlat=self.dlat,
            folder=temp_folder,
        )

        bnd_spec = xr.open_dataset(nc_file)

        bnd_spec = renormalize_era5_spec(bnd_spec)

        lon, lat = np.meshgrid(
            bnd_spec.longitude.values, bnd_spec.latitude.values[::-1]
        )
        lon = lon.ravel()
        lat = lat.ravel()

        # This spec is time, freq, dir, lat, lon
        spec = bnd_spec.d2fd.values
        # Latitude was flipped to be ascending, so flip that dimension
        spec = np.flip(spec, 3)

        freq = bnd_spec.frequency.values
        dirs = bnd_spec.direction.values
        try:
            time = bnd_spec.time.values
        except AttributeError:  # Name changes in new beta
            time = bnd_spec.valid_time.values

        # This is time, freq, dir, station
        spec = np.reshape(
            spec,
            (
                len(time),
                len(freq),
                len(dirs),
                len(lon),
            ),
        )
        # This is time, station, freq, dir (as we want it)
        spec = np.moveaxis(spec, 3, 1)

        # Inds given by point picker
        lon = lon[inds]
        lat = lat[inds]
        spec = spec[:, inds, :, :]

        coord_dict = {"lon": lon, "lat": lat, "time": time, "freq": freq, "dirs": dirs}
        data_dict = {"spec": spec}
        meta_dict = {"source": "ECMWF-ERA5 from Copernicus Climate Data Store"}

        return coord_dict, data_dict, meta_dict
