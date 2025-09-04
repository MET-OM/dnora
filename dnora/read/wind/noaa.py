from dnora.type_manager.data_sources import DataSource
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.file_structure import FileStructure
from dnora.read.ds_read_functions import basic_xarray_read
from functools import partial
import geo_parameters as gp
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("NOAA", "noaa", "wind")
class PacIOOS(ProductReader):
    """Reads GFS wind data from the best estimate files. 3h data for 5 days.

    Creator: NOAA National Centers for Environmental Prediction (NCEP)
    Publisher: Pacific Islands Ocean Observing System (PacIOOS)

    U.S. National Oceanic and Atmospheric Administration (NOAA) National Centers for Environmental Prediction (NCEP) Global Forecast System (GFS) numerical weather prediction model 8-day, 3-hourly global forecast at approximately 50-km or 0.5-deg
    """

    product_configuration = ProductConfiguration(
        filename="NCEP_Global_Atmospheric_Model_best.ncd",
        default_folders={
            DataSource.REMOTE: "https://pae-paha.pacioos.hawaii.edu/thredds/dodsC/ncep_global/",
        },
        data_vars=["ugrd10m", "vgrd10m"],
        ds_aliases={"ugrd10m": gp.wind.XWind, "vgrd10m": gp.wind.YWind},
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure()


@deprecated_class_call("NOAA", "noaa", "wind")
class NCEP(ProductReader):
    """Reads GFS wind data from the NCEP thredds server

    3h data for 16 days
    """

    product_configuration = ProductConfiguration(
        filename="gfs_0p25_%Hz",
        default_folders={
            DataSource.REMOTE: "https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs%Y%m%d/",
        },
        ds_creator_function=partial(
            basic_xarray_read, time_reference_str="days since 0001-01-01 00:00:0.0"
        ),
        data_vars=["ugrd10m", "vgrd10m"],
        ds_aliases={"ugrd10m": gp.wind.XWind, "vgrd10m": gp.wind.YWind},
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(stride=6, hours_per_file=385)


@deprecated_class_call("NOAA", "noaa", "wind")
class NCEP1h(ProductReader):
    """Reads GFS wind data from the NCEP thredds server

    1h data for 5 days
    """

    product_configuration = ProductConfiguration(
        filename="gfs_0p25_1hr_%Hz",
        default_folders={
            DataSource.REMOTE: "https://nomads.ncep.noaa.gov/dods/gfs_0p25_1hr/gfs%Y%m%d/",
        },
        ds_creator_function=partial(
            basic_xarray_read, time_reference_str="days since 0001-01-01 00:00:0.0"
        ),
        data_vars=["ugrd10m", "vgrd10m"],
        ds_aliases={"ugrd10m": gp.wind.XWind, "vgrd10m": gp.wind.YWind},
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(stride=6, hours_per_file=121)
