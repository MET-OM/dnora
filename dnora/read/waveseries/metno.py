from dnora.type_manager.data_sources import DataSource
from dnora.read.waveseries import WW3Unstruct
import geo_parameters as gp
from dnora.read.abstract_readers import PointDataReader


class NORAC(WW3Unstruct):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/norac_wave/field"
    }
    _default_filename = "ww3.%Y%m.nc"
    _decode_cf = False
    _keep_gp_names = True
    _data_vars = [
        gp.wave.Hs,
        gp.wave.Tm01("t01"),
        gp.wave.Tm02("t02"),
        gp.wave.Tm_10("t0m1"),
        gp.wave.Dirm("dir"),
        gp.wave.Dirp("dp"),
    ]

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE


class E39(PointDataReader):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/%Y/%m"
    }

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __init__(self, loc: str = "D", mode="wave"):
        self._loc = loc  # Given as "D", or "D_Breisundet"
        self._mode = mode  # 'wind' or 'wave'
        self._default_filename = f"%Y%m_E39_{self.loc()}_{self._mode}.nc"

    def _buoy_dict(self) -> dict:
        return {
            "A": "A_Sulafjorden",
            "B": "B_Sulafjorden",
            "B1": "B1_Sulafjorden",
            "C": "C_Sulafjorden",
            "C1": "C1_Sulafjorden",
            "D": "D_Breisundet",
            "F": "F_Vartdalsfjorden",
            "G": "G_Halsafjorden",
        }

    def loc(self) -> list[str]:
        if len(self._loc) > 2:
            return self._loc
        else:
            return self._buoy_dict()[self._loc]

    # def _folder_filename(
    #     self, source: DataSource, folder: str, filename: str
    # ) -> tuple[str]:
    #     if source == DataSource.REMOTE:
    #         folder = "https://thredds.met.no/thredds/dodsC/obs/buoy-svv-e39/%Y/%m"
    #     if filename is None:
    #         filename = f"%Y%m_E39_{self.loc()}_{self._mode}.nc"
    #     return folder, filename

    def get_coordinates(
        self,
        grid,
        start_time,
        source: DataSource,
        folder: str,
        filename: str,
        **kwargs,
    ) -> dict:
        # start_time = pd.to_datetime(start_time)
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)
        ds = read_first_ds(folder, filename, start_time)
        return {"lon": ds.longitude.values, "lat": ds.latitude.values}

    def __call__(
        self,
        grid,
        start_time,
        end_time,
        source: DataSource,
        folder: str,
        filename: str,
        inds,
        **kwargs,
    ) -> tuple:
        # loc = np.array(self._buoys())[inds][0]

        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times = create_monthly_stamps(start_time, end_time)
        file_times = start_times

        ds_creator_function = ds_xarray_read
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )
        # breakpoint()
        # ds_list = []
        # for month in months:
        #     url = get_url(folder, filename, month)
        #     ds_list.append(xr.open_dataset(url))
        ds = xr.concat(ds_list, dim="time")
        ds["lon"] = np.nanmedian(ds.longitude.values)
        ds["lat"] = np.nanmedian(ds.latitude.values)

        lon, lat, x, y = utils.get_coordinates_from_ds(ds)
        data_dict = {}
        for var in ds.data_vars:
            if var not in ["lon", "lat", "longitude", "latitude", "x", "y"]:
                if hasattr(ds[var], "standard_name"):
                    meta_param = gp.get(
                        ds[var].standard_name
                    )  # Find geo-parameter based on standard name
                else:
                    meta_param = None
                if meta_param is not None:
                    data_dict[meta_param] = np.swapaxes(
                        np.expand_dims(ds.get(var).values, axis=1), 0, 1
                    )

        coord_dict = {"time": ds.time.values, "lon": lon, "lat": lat}
        meta_dict = ds.attrs

        return coord_dict, data_dict, meta_dict
