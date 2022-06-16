from abc import ABC,  abstractmethod
import xarray as xr
# Import objects
from ..grd.grd_mod import Grid
from .. import file_module
from .. import aux_funcs
from .. import msg

class ForcingReader(ABC):
    """Reads forcing data from some source and provide it to the object.

    The area is defined from the Grid object that is passed.
    """
    @abstractmethod
    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        pass

    def name(self) -> str:
        return type(self).__name__

class DnoraNc(ForcingReader):
    def __init__(self, files: str) -> None:
        self.files = files

    def __call__(self, grid, start_time, end_time, expansion_factor):
        def _crop(ds):
            return ds.sel(time=slice(start_time, end_time), lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
        lon_min, lon_max, lat_min, lat_max = aux_funcs.expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)
        msg.info(f"Getting wind forcing from cached netcdf (e.g. {self.files[0]}) from {start_time} to {end_time}")

        # These files might get deleted, so we don't want to use dask for a lazy load
        ds = xr.open_mfdataset(self.files, preprocess=_crop)
        return ds


class File_WW3Nc(ForcingReader):
    def __init__(self, folder: str='', filename: str='ww3_wind_T0', dateformat: str='%Y%m%dT%H%M', stride: int=None, hours_per_file: int=24, last_file: str='', lead_time: int=0) -> None:
        self.stride = stride
        self.hours_per_file = hours_per_file
        self.lead_time = lead_time
        self.last_file = last_file

        if (not folder == '') and (not folder[-1] == '/'):
            self.folder = folder + '/'
        else:
            self.folder = folder

        self.filename = filename
        self.dateformat = dateformat

        return

    def __call__(self, grid: Grid, start_time: str, end_time: str, expansion_factor: float):
        """Reads in all wind data from a WW3 style wind input"""

        if self.stride is None:  # Read everything from one file
            start_times = [start_time]
            end_times = [end_time]
            file_times = [start_time]
        else:
            start_times, end_times, file_times = aux_funcs.create_time_stamps(start_time, end_time, stride = self.stride, hours_per_file = self.hours_per_file, last_file = self.last_file, lead_time = self.lead_time)

        lon_min, lon_max, lat_min, lat_max = aux_funcs.expand_area(min(grid.lon()), max(grid.lon()), min(grid.lat()), max(grid.lat()), expansion_factor)

        msg.info(f"Getting wind data from {self.filename} from {start_time} to {end_time}")
        wnd_list = []

        for time0, time1, file_time in zip(start_times, end_times, file_times):
            filename = self.get_filename(file_time)
            msg.from_file(filename)
            msg.plain(f"Reading wind data: {time0}-{time1}")

            wnd_list.append(xr.open_dataset(filename).sel(time=slice(time0, time1),
                            lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max)))

        wind_forcing = xr.concat(wnd_list, dim="time")

        return wind_forcing

    def get_filename(self, time) -> str:
        filename = self.folder + file_module.replace_times(self.filename,
                                                        self.dateformat,
                                                        [time]) + '.nc'
        return filename
