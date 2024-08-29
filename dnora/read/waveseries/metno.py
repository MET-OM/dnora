from dnora.type_manager.data_sources import DataSource
from dnora.waveseries.read import WW3Unstruct
import geo_parameters as gp


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
