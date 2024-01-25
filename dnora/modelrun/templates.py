# from ..run import model_executers
from dnora.modelrun import ModelRun

from dnora import spectra, wind, waterlevel

from dnora.dnora_types import DnoraDataType


class NORA3(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: spectra.read_metno.NORA3(),
        DnoraDataType.WIND: wind.read_metno.NORA3(),
    }


class ERA5(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: spectra.read_ec.ERA5(),
        DnoraDataType.WIND: wind.read_ec.ERA5(),
        DnoraDataType.WATERLEVEL: waterlevel.read_ec.GTSM_ERA5(),
    }


class WAM4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: spectra.read_metno.WAM4km(),
        DnoraDataType.WIND: wind.read_metno.MEPS(),
    }


class WW3_4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: spectra.read_metno.WW3_4km(),
        DnoraDataType.WIND: wind.read_metno.MEPS(),
    }
