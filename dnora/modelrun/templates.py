# from ..run import model_executers
from .modelrun import ModelRun

from .. import spectra, wind, waterlevel, current, ice

from ..dnora_object_type import DnoraObjectType


class NORA3(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: spectra.read_metno.NORA3(),
        DnoraObjectType.Forcing: wind.read_metno.NORA3(),
    }


class ERA5(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: spectra.read_ec.ERA5(),
        DnoraObjectType.Forcing: wind.read_ec.ERA5(),
        DnoraObjectType.WaterLevel: waterlevel.read_ec.GTSM_ERA5(),
    }


class WAM4km(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: spectra.read_metno.WAM4km(),
        DnoraObjectType.Forcing: wind.read_metno.MEPS(),
    }


class WW3_4km(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: spectra.read_metno.WW3_4km(),
        DnoraObjectType.Forcing: wind.read_metno.MEPS(),
    }
