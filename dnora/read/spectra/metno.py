# Import abstract classes and needed instances of them
from dnora.process.spectra import RemoveEmpty
from dnora.type_manager.spectral_conventions import SpectralConvention

# Import aux_funcsiliry functions
from dnora.type_manager.data_sources import DataSource
from dnora.read.spectra import WW3, WAM


class WAM4km(WAM):
    """Operational WAM for Norwegian waters"""

    WAM_OTHER_VARS = [
        "time",
        "longitude",
        "latitude",
        "ff",
        "dd",
        "Pdir",
        "hs",
        "tp",
        "thq_sea",
        "thq_swell",
        "depth",
    ]
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam4archive/%Y/%m/%d",
    }
    _default_filename = "MyWave_wam4_SPC_%Y%m%dT%HZ.nc"
    stride: int = 6
    hours_per_file: int = 73

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE

    def post_processing(self):
        return RemoveEmpty()


class NORA3(WAM):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/windsurfer/mywavewam3km_spectra/%Y/%m",
        DataSource.INTERNAL: "WINDSURFER/mw3hindcast/spectra/%Y/%m",
    }
    _default_filename = "SPC%Y%m%d00.nc"
    stride: int = 24
    hours_per_file: int = 24

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE


# This WW3 product uses WAM-conventions in the spectral files with respect to variable names etc
class WW3_4km(WAM):
    """Operational WW3 for Norwegian waters"""

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/ww3_4km_archive_files/%Y/%m/%d",
        DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
    }
    stride: int = 6
    hours_per_file: int = 73

    def __init__(self, tile="POI") -> None:
        super().__init__()
        self._default_filename = f"ww3_4km_{tile}_SPC_%Y%m%dT%HZ.nc"

    def post_processing(self):
        return RemoveEmpty()

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE


class WAM3(WAM):
    """covers Nordic Seas and the Arctic"""

    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam3_latest/",
        DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
    }
    _default_filenames = {DataSource.REMOTE: "MyWave_wam3_WAVE_%Y%m%dT%HZ.nc"}
    _default_filename = "MyWave_wam3_SPC_%Y%m%dT%HZ.nc"

    stride: int = 12
    hours_per_file: int = 121
    offset: int = 6

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE


class WAM800(WAM):
    """c0 covers Finnmark, c1 covers NordNorge, c2 covers MidtNorge, c3 covers Vestlandet and c4 covers Skagerrak (these are the names of the different domains as used below)."""

    tile_names = {
        "c0": "Finnmark",
        "c1": "NordNorge",
        "c2": "MidtNorge",
        "c3": "Vestlandet",
        "c4": "Skagerrak",
    }
    WAM_OTHER_VARS = [
        "time",
        "longitude",
        "latitude",
        "ff",
        "dd",
        "Pdir",
        "hs",
        "tp",
        "thq_sea",
        "thq_swell",
        "depth",
    ]
    stride: int = 12
    hours_per_file: int = 73

    def __init__(
        self,
        tile="c3",
    ) -> None:
        super().__init__()
        self._default_folders = {
            DataSource.REMOTE: f"https://thredds.met.no/thredds/dodsC/fou-hi/mywavewam800{self.tile_names[tile][0].lower()}",
            DataSource.IMMUTABLE: "DNMI_WAVE/%Y/%m/%d",
        }
        self._default_filenames = {DataSource.REMOTE: f"MyWave_wam800_{tile}SPC%H.nc"}
        self._default_filename = f"MyWave_wam800_{tile}SPC_%Y%m%dT%HZ.nc"

    def convention(self) -> SpectralConvention:
        return SpectralConvention.OCEAN

    def default_data_source(self) -> DataSource:
        return DataSource.IMMUTABLE


class NORAC(WW3):
    _default_folders = {
        DataSource.REMOTE: "https://thredds.met.no/thredds/dodsC/norac_wave/spec/",
        DataSource.INTERNAL: "sfiblues/wave_hindcast/hindcast_v2/spec",
    }

    def default_data_source(self) -> DataSource:
        return DataSource.REMOTE
