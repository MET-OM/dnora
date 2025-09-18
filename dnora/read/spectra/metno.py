# Import abstract classes and needed instances of them
from dnora.process.spectra import RemoveEmpty, RemoveNanTimes
from dnora.type_manager.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora.type_manager.data_sources import DataSource
from dnora.read.product_readers import SpectralProductReader

from dnora.read.product_configuration import ProductConfiguration
from dnora.read.file_structure import FileStructure
from functools import partial
from dnora.read.ds_read_functions import basic_xarray_read
import geo_parameters as gp
from dnora.read.aliases import WW3_DS_ALIASES
from dnora.read.depreciation_decorator import deprecated_class_call


@deprecated_class_call("MET Norway's", "metno", "spectra")
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


@deprecated_class_call("MET Norway's", "metno", "spectra")
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


@deprecated_class_call("MET Norway's", "metno", "spectra")
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
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.REMOTE,
    )

    file_structure = FileStructure(
        stride=6,
        hours_per_file=73,
    )

    def post_processing(self):
        return RemoveEmpty()


@deprecated_class_call("MET Norway's", "metno", "spectra")
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
        default_data_source=DataSource.REMOTE,
        ds_aliases={"SPEC": gp.wave.Efth},
    )

    file_structure = FileStructure(
        stride=12,
        hours_per_file=73,
    )


@deprecated_class_call("MET Norway's", "metno", "spectra")
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
        ds_aliases=WW3_DS_ALIASES,
    )


@deprecated_class_call("MET Norway's", "metno", "spectra")
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


@deprecated_class_call("MET Norway's", "metno", "spectra")
class CLIMAREST(SpectralProductReader):
    product_configuration = ProductConfiguration(
        filename="CLIMAREST_*_spec_2040_2070.nc",
        convention=SpectralConvention.OCEAN,
        default_data_source=DataSource.LOCAL,
        ds_aliases={
            "efth": gp.wave.Efth,
            "dir": gp.wave.Dirs,
            "hs": gp.wave.Hs,
            "tp": gp.wave.Tp,
            "Pdir": gp.wave.DirpTo,
        },
    )

    # stride = None means one point per file and all times in one file
    file_structure = FileStructure(
        stride=None,
    )

    def post_processing(self):
        return RemoveNanTimes()
