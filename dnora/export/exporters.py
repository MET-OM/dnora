from .exporter import DataExporter
from . import grid, spectra, wind, spectra1d, waveseries, waterlevel, inputfile
from .generic.generic_writers import DnoraNc, Null, DumpToNc

from ..dnora_types import DnoraDataType
from .exporter import WriterFunction
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
        DnoraDataType.SPECTRA: spectra.SWAN(),
        DnoraDataType.WIND: wind.SWAN(),
        DnoraDataType.GRID: grid.SWAN(),
        # DnoraDataType.InputFile: inputfile.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra.SWAN(),
        DnoraDataType.WIND: wind.SWAN(),
        DnoraDataType.GRID: grid.SWAN(),
        # DnoraDataType.InputFile: inputfile.SWASH(),
    }

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra.WW3(),
        DnoraDataType.WIND: wind.WW3(),
        DnoraDataType.GRID: grid.WW3(),
        DnoraDataType.TRIGRID: grid.WW3Triangular(),
        # DnoraDataType.InputFile: inputfile.WW3(),
    }

    def _get_default_format(self):
        return ModelFormat.WW3


class HOS_ocean(DataExporter):
    _writer_dict = {
        # DnoraDataType.InputFile: inputfile.HOS_ocean(),
    }

    def _get_default_format(self):
        return ModelFormat.HOS_Ocean


class REEF3D(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA1D: spectra1d.REEF3D(),
        # DnoraDataType.InputFile: inputfile.REEF3D(),
    }

    def _get_default_format(self):
        return ModelFormat.REEF3D
