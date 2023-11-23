from .exp_mod import DataExporter
from . import grd, bnd, wnd, spc, wsr, wlv, inp
from .general.general_writing_functions import DnoraNc, Null, DumpToNc

from ..dnora_object_type import DnoraObjectType
from .exp_mod import WriterFunction
from ..model_formats import ModelFormat


class NullExporter(DataExporter):
    _writer_dict = {}

    def _get_default_writer(self) -> WriterFunction:
        return Null()


class Cacher(DataExporter):
    def _get_default_writer(self) -> WriterFunction:
        return DnoraNc()

    def _get_default_format(self) -> str:
        return ModelFormat.CACHE


class SWAN(DataExporter):
    _writer_dict = {
        DnoraObjectType.Boundary: bnd.SWAN(),
        DnoraObjectType.Forcing: wnd.SWAN(),
        DnoraObjectType.Grid: grd.SWAN(),
        DnoraObjectType.InputFile: inp.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(DataExporter):
    _writer_dict = {
        DnoraObjectType.Boundary: bnd.SWAN(),
        DnoraObjectType.Forcing: wnd.SWAN(),
        DnoraObjectType.Grid: grd.SWAN(),
        DnoraObjectType.InputFile: inp.SWASH(),
    }

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(DataExporter):
    _writer_dict = {
        DnoraObjectType.Boundary: bnd.WW3(),
        DnoraObjectType.Forcing: wnd.WW3(),
        DnoraObjectType.Grid: grd.WW3(),
        DnoraObjectType.TriGrid: grd.WW3Triangular(),
        DnoraObjectType.InputFile: inp.WW3(),
    }

    def _get_default_format(self):
        return ModelFormat.WW3


class HOS_ocean(DataExporter):
    _writer_dict = {
        DnoraObjectType.InputFile: inp.HOS_ocean(),
    }

    def _get_default_format(self):
        return ModelFormat.HOS_Ocean


class REEF3D(DataExporter):
    _writer_dict = {
        DnoraObjectType.Spectra: spc.REEF3D(),
        DnoraObjectType.InputFile: inp.REEF3D(),
    }

    def _get_default_format(self):
        return ModelFormat.REEF3D
