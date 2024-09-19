# from ..run import model_executers
from dnora.modelrun.modelrun import ModelRun
from dnora.read.generic import ConstantData
from dnora.type_manager.dnora_types import DnoraDataType
import dnora


class Constant(ModelRun):
    _reader_dict = {
        DnoraDataType.WIND: ConstantData(),
        DnoraDataType.SPECTRA: ConstantData(),
        DnoraDataType.SPECTRA1D: ConstantData(),
        DnoraDataType.ICE: ConstantData(),
        DnoraDataType.WATERLEVEL: ConstantData(),
        DnoraDataType.CURRENT: ConstantData(),
        DnoraDataType.WAVESERIES: ConstantData(),
    }


class NORA3(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.NORA3(),
        DnoraDataType.WIND: dnora.read.wind.metno.NORA3(),
        DnoraDataType.ICE: dnora.read.ice.metno.NORA3(),
    }


class ERA5(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.ec.ERA5(),
        DnoraDataType.WIND: dnora.read.wind.ec.ERA5(),
        DnoraDataType.WATERLEVEL: dnora.read.waterlevel.ec.GTSM_ERA5(),
    }


class NOAA(ModelRun):
    _reader_dict = {
        DnoraDataType.WIND: dnora.read.wind.noaa.GFS(),
    }


class WAM4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.WAM4km(),
        DnoraDataType.WIND: dnora.read.wind.metno.MEPS(),
    }


class WW3_4km(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.metno.WW3_4km(),
        DnoraDataType.WIND: dnora.read.wind.metno.MEPS(),
    }

class NCHMF(ModelRun):
    _reader_dict = {
        DnoraDataType.SPECTRA: dnora.read.spectra.nchmf.SWAN(),
        DnoraDataType.WIND: dnora.read.wind.nchmf.Oper(),
        DnoraDataType.WATERLEVEL: dnora.read.waterlevel.cmems.Global(),
    }
