from dnora.dnora_type_manager.data_sources import DataSource
from dnora.dnora_type_manager.dnora_types import DnoraDataType
from dnora.readers.abstract_readers import DataReader
from dnora.readers.fimex_functions import ds_fimex_read
from dnora.readers.ds_read_functions import read_ds_list, setup_temp_dir
from dnora.readers.file_structure import FileStructure
from dnora.grid import Grid
from dnora import msg
from dnora.aux_funcs import expand_area
from functools import partial
import xarray as xr


class NORA3(DataReader):
    """Reads ice data of the NORA3 wave hindcast directly from MET Norways servers.

    The NORA3 high-resolution (ca 3 km) hindcast for the Nordic Seas and the Arctic Ocean.

    Breivik, Ø., Carrasco, A., Haakenstad, H., Aarnes, O. J., Behrens, A., Bidlot, J.-R., et al. (2022).
    The impact of a reduced high-wind Charnock parameter on wave growth with application to the North Sea,
    the Norwegian Sea, and the Arctic Ocean.
    Journal of Geophysical Research: Oceans, 127, e2021JC018196. https://doi.org/10.1029/2021JC018196
    """

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m"
    }
    _default_filename = "%Y%m%d_MyWam3km_hindcast.nc"

    def __init__(
        self,
        stride: int = 24,
        hours_per_file: int = 24,
        last_file: str = "",
        lead_time: int = 0,
    ):
        """The data is currently in daily files. Do not change the default
        setting unless you have a good reason to do so.
        """

        self.file_structure = FileStructure(
            stride=stride,
            hours_per_file=hours_per_file,
            last_file=last_file,
            lead_time=lead_time,
        )
        return

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def __call__(
        self,
        grid: Grid,
        start_time: str,
        end_time: str,
        source: DataSource,
        folder: str,
        filename: str,
        expansion_factor: float = 1.2,
        program: str = "pyfimex",
        **kwargs,
    ):
        """Reads in all boundary spectra between the given times and at for the given indeces"""
        folder = self._folder(folder, source)
        filename = self._filename(filename, source)

        start_times, end_times, file_times = self.file_structure.create_time_stamps(
            start_time, end_time
        )
        msg.info(f"Getting ice forcing from NORA3 from {start_time} to {end_time}")
        setup_temp_dir(DnoraDataType.ICE, self.name())
        # Define area to search in
        msg.info(f"Using expansion_factor = {expansion_factor:.2f}")
        lon, lat = expand_area(grid.edges("lon"), grid.edges("lat"), expansion_factor)

        msg.process(f"Applying {program}")
        ds_creator_function = partial(
            ds_fimex_read,
            lon=lon,
            lat=lat,
            resolution_in_km=3,
            data_vars=["SIC", "SIT"],
            data_type=DnoraDataType.ICE,
            name=self.name(),
            program=program,
        )
        ds_list = read_ds_list(
            start_times,
            end_times,
            file_times,
            folder,
            filename,
            ds_creator_function,
        )

        ice_forcing = xr.concat(ds_list, dim="time")

        data_dict = {
            "sic": ice_forcing.SIC.data,
            "sit": ice_forcing.SIT.data,
        }
        coord_dict = {
            "time": ice_forcing.time.data,
            "lon": ice_forcing.rlon.data,
            "lat": ice_forcing.rlat.data,
        }
        meta_dict = ice_forcing.attrs

        return coord_dict, data_dict, meta_dict
