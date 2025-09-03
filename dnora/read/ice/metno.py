from dnora.type_manager.data_sources import DataSource
from dnora.read.fimex_functions import ds_fimex_read
from dnora.read.file_structure import FileStructure
from dnora.grid import Grid
from functools import partial
import re
from dnora.read.product_readers import ProductReader
from dnora.read.product_configuration import ProductConfiguration
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("MET Norway's", "metno", "ice")
class NORA3(ProductReader):
    """Reads wind data (from monthly files 'arome3km_1hr_YYYYMM.nc') of the NORA3 hindcast directly from MET Norways servers.

    The NORA3 HARMONIE-AROME high-resolution (ca 3 km) hindcast for the
    North Sea, the Norwegian Sea, and the Barents Sea.

    Haakenstad, H., Breivik, Ø., Furevik, B. R., Reistad, M., Bohlinger, P., &
    Aarnes, O. J. (2021). NORA3: A Nonhydrostatic High-Resolution Hindcast of
    the North Sea, the Norwegian Sea, and the Barents Sea,
    Journal of Applied Meteorology and Climatology, 60(10), 1443-1464,
    DOI: 10.1175/JAMC-D-21-0029.1
    """

    product_configuration = ProductConfiguration(
        filename="%Y%m%d_MyWam3km_hindcast.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_files/%Y/%m",
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3,
        ),
        data_vars=["SIC", "SIT"],
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=24,
        hours_per_file=24,
    )


@deprecated_class_call("MET Norway's", "metno", "ice")
class Barents25(ProductReader):
    """Reads sea ice data of the Barents 2.5 km operational ocean model directly from MET Norways servers.

    Barents-2.5km is an operational data-assimilative coupled ocean and sea ice ensemble prediction model for
    the Barents Sea and Svalbard

    Röhrs, J., Gusdal, Y., Rikardsen, E. S. U., Durán Moro, M., Brændshøi, J., Kristensen, N. M.,
    Fritzner, S., Wang, K., Sperrevik, A. K., Idžanović, M., Lavergne, T., Debernard, J. B., and Christensen, K. H.:
    Barents-2.5km v2.0: an operational data-assimilative coupled ocean and sea ice ensemble prediction model
    for the Barents Sea and Svalbard, Geosci. Model Dev., 16, 5401-5426,
    https://doi.org/10.5194/gmd-16-5401-2023, 2023.
    """

    product_configuration = ProductConfiguration(
        filename="barents_sfc_%Y%m%dT[6]Zm00.nc",
        default_folders={
            DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/barents_eps_surface/%Y/%m/%d/T[6]Z"
        },
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=2.5,
        ),
        data_vars=["ice_concentration", "ice_thickness"],
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=97,
    )


@deprecated_class_call("MET Norway's", "metno", "ice")
class CLIMAREST(ProductReader):
    product_configuration = ProductConfiguration(
        filename="sic_HCLIM43_MPIESM12LR_3hr_%Y_%m.nc",
        ds_creator_function=partial(
            ds_fimex_read,
            resolution_in_km=3.0,
        ),
        data_vars=["sic"],
        default_data_source=DataSource.LOCAL,
    )

    file_structure = FileStructure(
        stride="month",
    )
