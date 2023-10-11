# from ..run import model_executers
from .mdl_mod import ModelRun

from .. import bnd, wnd, wlv, ocr, ice


class NORA3(ModelRun):
    _reader_dict = {
        "Boundary": bnd.read_metno.NORA3(),
        "Forcing": wnd.read_metno.NORA3(),
    }


class ERA5(ModelRun):
    _reader_dict = {
        "Boundary": bnd.read_ec.ERA5(),
        "Forcing": wnd.read_ec.ERA5(),
        "WaterLevel": wlv.read_ec.GTSM_ERA5(),
    }


class WAM4km(ModelRun):
    _reader_dict = {
        "Boundary": bnd.read_metno.WAM4km(),
        "Forcing": wnd.read_metno.MEPS(),
    }


class WW3_4km(ModelRun):
    _reader_dict = {
        "Boundary": bnd.read_metno.WW3_4km(),
        "Forcing": wnd.read_metno.MEPS(),
    }
