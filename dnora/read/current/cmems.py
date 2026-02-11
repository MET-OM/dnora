from dnora.type_manager.data_sources import DataSource
from dnora.read.file_structure import FileStructure

from functools import partial
from dnora.read.cmems_functions import ds_cmems_read
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.process.gridded import FillNaNs
from dnora.read.depreciation_decorator import deprecated_class_call
import geo_parameters as gp

@deprecated_class_call("CMEMS", "cmems", "current")
class Global(ProductReader):
    """This product is a L4 REP and NRT global total velocity field at 0m and 15m together wiht its individual components 
    (geostrophy and Ekman) and related uncertainties. It consists of the zonal and meridional velocity at a 1h frequency 
    and at 1/4 degree regular grid. The total velocity fields are obtained by combining CMEMS satellite Geostrophic surface 
    currents and modelled Ekman currents at the surface and 15m depth (using ERA5 wind stress in REP and ERA5* in NRT). 
    1 hourly product, daily and monthly means are available. This product has been initiated in the frame of CNES/CLS projects. 
    Then it has been consolidated during the Globcurrent project (funded by the ESA User Element Program)..

    DOI (product): https://doi.org/10.48670/mds-00327
    https://https://data.marine.copernicus.eu/product/MULTIOBS_GLO_PHY_MYNRT_015_003/description
    """

    product_configuration = ProductConfiguration(
        ds_creator_function=partial(
            ds_cmems_read,
            dataset_id="cmems_obs-mob_glo_phy-cur_my_0.25deg_PT1H-i",
            variables=["uo","vo"],
            minimum_depth=0,
            maximum_depth=0,
        ),
        ds_aliases={"uo": gp.ocean.XCurrent, "vo": gp.ocean.YCurrent},
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure()

    def post_processing(self):
        return FillNaNs(0)