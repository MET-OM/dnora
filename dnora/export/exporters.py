from .exp_mod import DataExporter
from . import grid, spectra, wind, spectra1d, waveseries, waterlevel, inputfile
from .generic.generic_writers import DnoraNc, Null, DumpToNc

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
        DnoraObjectType.Boundary: spectra.SWAN(),
        DnoraObjectType.Forcing: wind.SWAN(),
        DnoraObjectType.Grid: grid.SWAN(),
        DnoraObjectType.InputFile: inputfile.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(DataExporter):
    _writer_dict = {
        DnoraObjectType.Boundary: spectra.SWAN(),
        DnoraObjectType.Forcing: wind.SWAN(),
        DnoraObjectType.Grid: grid.SWAN(),
        DnoraObjectType.InputFile: inputfile.SWASH(),
    }

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(DataExporter):
    _writer_dict = {
        DnoraObjectType.Boundary: spectra.WW3(),
        DnoraObjectType.Forcing: wind.WW3(),
        DnoraObjectType.Grid: grid.WW3(),
        DnoraObjectType.TriGrid: grid.WW3Triangular(),
        DnoraObjectType.InputFile: inputfile.WW3(),
    }

    def _get_default_format(self):
        return ModelFormat.WW3


class HOS_ocean(DataExporter):
    _writer_dict = {
        DnoraObjectType.InputFile: inputfile.HOS_ocean(),
    }

    def _get_default_format(self):
        return ModelFormat.HOS_Ocean


class REEF3D(DataExporter):
    _writer_dict = {
        DnoraObjectType.Spectra: spectra1d.REEF3D(),
        DnoraObjectType.InputFile: inputfile.REEF3D(),
    }

    def _get_default_format(self):
        return ModelFormat.REEF3D
