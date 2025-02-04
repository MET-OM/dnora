# Import abstract classes and needed instances of them
from dnora.process.spectra import RemoveEmpty
from dnora.type_manager.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora.type_manager.data_sources import DataSource
from dnora.read.spectra import SpectralProductReader

from dnora.read.product_configuration import ProductConfiguration
from dnora.read.file_structure import FileStructure
from functools import partial
from dnora.read.ds_read_functions import basic_xarray_read
import geo_parameters as gp


class WAM4km(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="MyWave_wam4_SPC_%Y%m%dT%HZ.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m/%d",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="x"),
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.REMOTE,
        ds_aliases={"SPEC": gp.wave.Efth},
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=73,
    )

    def post_processing(self):
        return RemoveEmpty()


class NORA3(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="SPC%Y%m%d00.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m",
            DataSource.INTERNAL: "WINDSURFER/mw3hindcast/spectra/%Y/%m",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="x"),
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.REMOTE,
        ds_aliases={"SPEC": gp.wave.Efth},
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )


class WW3_4km(SpectralProductReader):
    """Operational WW3 for Norwegian waters"""

    product_configuration = ProductConfiguration(
        filename="ww3_4km_#TILE_SPC_%Y%m%dT%HZ.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/%Y/%m/%d",
            DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
        },
        tile="POI",
        ds_creator_function=partial(basic_xarray_read, inds_var="x"),
        convention=SpectralConvention.WW3,
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=73,
    )

    def post_processing(self):
        return RemoveEmpty()


class WAM800(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="MyWave_wam800_#TILESPC_%Y%m%dT%HZ.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800#TILENAME",
            DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
        },
        tile="c3",
        tile_names={
            "c0": "Finnmark",
            "c1": "NordNorge",
            "c2": "MidtNorge",
            "c3": "Vestlandet",
            "c4": "Skagerrak",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="x"),
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.IMMUTABLE,
        ds_aliases={"SPEC": gp.wave.Efth},
    )

    file_structure = FileStructure(
        stride=12,
        hours_per_file=73,
    )


class NORAC(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="ww3_spec.%Y%m.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/norac_wave/spec/",
            DataSource.INTERNAL: "sfiblues/wave_hindcast/hindcast_v2/spec",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="station"),
        convention=SpectralConvention.WW3,
        default_data_source=DataSource.REMOTE,
    )

    # This defines how the file structure of the model output is set up
    # file_structure = FileStructure(stride="month")

    # _default_filename = "ww3_spec.%Y%m.nc"

    # stride = "month"  # int (for hourly), or 'month'
    # hours_per_file = None  # int (if not monthly files)
    # offset = 0  # int

    # def set_up_for_ds_read(self, obj_type) -> tuple:
    #     if obj_type == DnoraDataType.WAVESERIES:
    #         dynamic = True
    #         ignore_vars = ["station_name"]
    #         ds_aliases = {"dpt": gp.ocean.WaterDepth}
    #         aliases = {}
    #     else:
    #         ignore_vars = []
    #         ds_aliases = {}
    #         aliases = {}
    #         dynamic = False
    #     return dynamic, aliases, ds_aliases, ignore_vars

    # def _ds_creator_function(self, inds):
    #     return partial(basic_xarray_read, inds=inds, inds_var="station")


# class NORAC(WW3):
#     _default_folders = {
#         DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/norac_wave/spec/",
#         DataSource.INTERNAL: "sfiblues/wave_hindcast/hindcast_v2/spec",
#     }

#     def default_data_source(self) -> DataSource:
#         return DataSource.REMOTE


# class WAM3(WAM):
#     """covers Nordic Seas and the Arctic"""

#     _default_folders = {
#         DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam3_latest/",
#         DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
#     }
#     _default_filenames = {DataSource.REMOTE: "MyWave_wam3_WAVE_%Y%m%dT%HZ.nc"}
#     _default_filename = "MyWave_wam3_SPC_%Y%m%dT%HZ.nc"

#     stride: int = 12
#     hours_per_file: int = 121
#     offset: int = 6

#     def default_data_source(self) -> DataSource:
#         return DataSource.IMMUTABLE


class WAM3(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="MyWave_wam3_SPC_%Y%m%dT%HZ.nc",
        default_filenames={DataSource.REMOTE: "MyWave_wam3_WAVE_%Y%m%dT%HZ.nc"},
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam3_latest/",
            DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
        },
        ds_creator_function=partial(basic_xarray_read, inds_var="x"),
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.IMMUTABLE,
        ds_aliases={"SPEC": gp.wave.Efth},
    )

    file_structure = FileStructure(
        stride=12,
        hours_per_file=121,
        offset=6,
    )
