import numpy as np
import xarray as xr
import subprocess

# from subprocess import call, run
import os, glob

# Import objects
from dnora.grid import Grid

from scipy.interpolate import griddata

# Import abstract classes
from dnora.read.abstract_readers import DataReader

# Import aux_funcsiliry functions
from dnora import msg
from dnora import utils
import pandas as pd
from dnora.type_manager.data_sources import DataSource
from dnora.read.ds_read_functions import setup_temp_dir
from dnora.type_manager.dnora_types import DnoraDataType

from dnora.read.depreciation_decorator import deprecated_class_call


def download_GTSM_from_cds(start_time, end_time, folder="dnora_wlv_temp") -> str:
    """Downloads GTSM model water level data from the Copernicus Climate Data Store for a
    given area and time period"""
    start_time = pd.Timestamp(start_time)
    end_time = pd.Timestamp(end_time)
    try:
        import cdsapi
    except ImportError as e:
        msg.advice(
            "The cdsapi package is required to use ECWMF products! Install by e.g. 'conda install cdsapi'"
        )
        raise e

    c = cdsapi.Client()

    filename = f"{folder}/EC_GTSM_ERA5.tar.gz"
    # cds_command_test = {
    #     'product_type': 'reanalysis',
    #     'format': 'netcdf',
    #     'variable': [
    #         '10m_u_component_of_wind', '10m_v_component_of_wind',
    #     ],
    #     'year': '2008',
    #     'month': '01',
    #     'day': [
    #         '01', '02',
    #     ],
    #     'time': [
    #         '00:00', '01:00', '02:00',
    #         '03:00', '04:00', '05:00',
    #         '06:00', '07:00', '08:00',
    #         '09:00', '10:00', '11:00',
    #         '12:00', '13:00', '14:00',
    #         '15:00', '16:00', '17:00',
    #         '18:00', '19:00', '20:00',
    #         '21:00', '22:00', '23:00',
    #     ],
    #     'area': [
    #         61.25, 4, 60.53,
    #         5.73,
    #     ],
    # }

    years = [f"{y:4.0f}" for y in utils.time.int_list_of_years(start_time, end_time)]
    if len(years) == 1:
        years = years[0]
    months = [f"{m:02.0f}" for m in utils.time.int_list_of_months(start_time, end_time)]
    if len(months) == 1:
        months = months[0]

    cds_command = {
        "format": "tgz",
        "variable": ["total_water_level"],
        "experiment": "reanalysis",
        "temporal_aggregation": "hourly",
        "year": years,  # 1979-2018
        "month": months,
    }

    c.retrieve("sis-water-level-change-timeseries-cmip6", cds_command, filename)
    return filename


@deprecated_class_call("ERA5", "era5", "waterlevel")
class GTSM_ERA5(DataReader):
    """Reads GTSM_ERA5 waterlevel data"""

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __call__(
        self,
        obj_type,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        expansion_factor: float = 1.1,
        **kwargs,
    ):
        """Reads hourly water level from GTSM_ERA5 database"""

        msg.info(f"Getting GTSM/ERA5 water level from {start_time} to {end_time}")
        temp_folder = setup_temp_dir(DnoraDataType.WATERLEVEL, self.name())

        # Define area to search in
        lon, lat = utils.grid.expand_area(
            lon=grid.edges("lon"),
            lat=grid.edges("lat"),
            expansion_factor=expansion_factor,
        )

        out_file = download_GTSM_from_cds(start_time, end_time, folder=temp_folder)

        temppath = os.path.dirname(out_file)
        # first unpack the tar.gz file.
        nc_file = (
            subprocess.run(["tar", "-ztf", out_file], stdout=subprocess.PIPE)
            .stdout.decode("utf-8")
            .split("\n")[0:-1]
        )
        nc_file = sorted([ff.strip("\r") for ff in nc_file])
        # print(nc_file)
        subprocess.run(
            ["tar", "-xzvf", out_file, "--directory", temppath], stdout=subprocess.PIPE
        )  # Extract tar file

        lon_local = np.arange(lon[0], lon[1], 0.1)
        lat_local = np.arange(lat[0], lat[1], 0.1)
        grid_x, grid_y = np.meshgrid(lon_local, lat_local, indexing="xy")

        print(nc_file)
        grid_tot = []
        time_tot = []
        for ncfile in nc_file:
            # print(os.path.join(temppath,nc_file))
            waterlevel = xr.open_dataset(
                os.path.join(temppath, ncfile), engine="netcdf4"
            )
            waterlevel = waterlevel.rename_vars(
                {"station_x_coordinate": "lon", "station_y_coordinate": "lat"}
            )
            waterlevel = waterlevel.sel(stations=waterlevel.lon >= lon[0])
            waterlevel = waterlevel.sel(stations=waterlevel.lon <= lon[1])
            waterlevel = waterlevel.sel(stations=waterlevel.lat >= lat[0])
            waterlevel = waterlevel.sel(stations=waterlevel.lat <= lat[1])

            waterlevel = waterlevel.sel(time=slice(start_time, end_time))
            # grid_x, grid_y = np.mgrid[lon_min:lon_max:100j, lat_min:lat_max:110j]

            points = np.array([waterlevel.lon, waterlevel.lat]).T
            time = waterlevel.time
            grid_z = np.zeros([len(time), len(lat_local), len(lon_local)])
            msg.info(
                "Interpolating data to a 0.1x0.1 degree resolution with cubic interpolation..."
            )
            for i_t, t in enumerate(time):
                values = waterlevel.waterlevel[i_t, :]
                grid_z[i_t, :, :] = griddata(
                    points, values, (grid_x, grid_y), method="cubic", fill_value=0.0
                )
                # plt.imshow(grid_z[i_t,:,:], extent=(lon_min, lon_max, lat_min, lat_max), origin='lower')
                # plt.scatter(waterlevel.lon, waterlevel.lat)
                # plt.colorbar()
                # plt.show()
            grid_tot.append(grid_z)
            time_tot.append(time)

        grid_tot = np.concatenate(grid_tot, axis=0)
        time_tot = np.concatenate(time_tot, axis=0)

        # Finally we put the new gridded data into a dataset
        # waterlevel_gridded = xr.Dataset(
        #     data_vars=dict(
        #         waterlevel=(["time", "lat", "lon"], grid_tot)
        #     ),
        #     coords=dict(
        #         time=(["time"], time_tot),
        #         lat=(["lat"], lat_local),
        #         lon=(["lon"], lon_local),
        #     ),
        #     attrs=dict(description="waterlevel"),
        # )

        # print(waterlevel_gridded)
        coord_dict = {"lon": lon_local, "lat": lat_local, "time": time}
        data_dict = {"eta": grid_tot}
        meta_dict = {"description": "Waterlevel from GTSM/ERA5"}

        return coord_dict, data_dict, meta_dict
        # return (
        #     time,
        #     grid_tot,
        #     lon_local,
        #     lat_local,
        #     None,
        #     None,
        #     dict(description="Waterlevel from GTSM/ERA5"),
        # )
