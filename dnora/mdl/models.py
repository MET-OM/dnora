# from ..run import model_executers
from .mdl_mod import ModelRun

from .. import bnd, wnd, wlv, ocr, ice

from ..dnora_object_type import DnoraObjectType


class NORA3(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: bnd.read_metno.NORA3(),
        DnoraObjectType.Forcing: wnd.read_metno.NORA3(),
    }


class ERA5(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: bnd.read_ec.ERA5(),
        DnoraObjectType.Forcing: wnd.read_ec.ERA5(),
        DnoraObjectType.WaterLevel: wlv.read_ec.GTSM_ERA5(),
    }


class WAM4km(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: bnd.read_metno.WAM4km(),
        DnoraObjectType.Forcing: wnd.read_metno.MEPS(),
    }


class WW3_4km(ModelRun):
    _reader_dict = {
        DnoraObjectType.Boundary: bnd.read_metno.WW3_4km(),
        DnoraObjectType.Forcing: wnd.read_metno.MEPS(),
    }
