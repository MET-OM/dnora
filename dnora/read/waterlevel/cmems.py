from dnora.type_manager.data_sources import DataSource
from dnora.read.file_structure import FileStructure

from functools import partial
from dnora.read.cmems_functions import ds_cmems_read
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.process.gridded import FillNaNs
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("CMEMS", "cmems", "waterlevel")
class Global(ProductReader):
    """The Operational Mercator global ocean analysis and forecast system at 1/12 degree is providing 10 days of 3D global
    ocean forecasts updated daily. The time series is aggregated in time in order to reach a two full year's time series
    sliding window. This product includes daily and monthly mean files of temperature, salinity, currents, sea level,
    mixed layer depth and ice parameters from the top to the bottom over the global ocean. It also includes hourly
    mean surface fields for sea level height, temperature and currents. The global ocean output files are displayed
    with a 1/12 degree horizontal resolution with regular longitude/latitude equirectangular projection.
    50 vertical levels are ranging from 0 to 5500 meters. This product also delivers a special dataset for
    surface current which also includes wave and tidal drift called SMOC (Surface merged Ocean Current).

    DOI (product): https://doi.org/10.48670/moi-00016
    https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
    """

    product_configuration = ProductConfiguration(
        ds_creator_function=partial(
            ds_cmems_read,
            dataset_id="cmems_mod_glo_phy_anfc_0.083deg_PT1H-m",
            variables=["zos"],
            minimum_depth=0.49402499198913574,
            maximum_depth=0.49402499198913574,
        ),
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure()

    def post_processing(self):
        return FillNaNs(0)
