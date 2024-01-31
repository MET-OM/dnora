from .exporter import DataExporter
from . import grid, spectra, wind, spectra1d, waveseries, waterlevel
from .generic import generic_writers

from dnora.dnora_types import DnoraDataType, DnoraFileType
from .exporter import WriterFunction
from dnora.model_formats import ModelFormat


class NullExporter(DataExporter):
    _writer_dict = {}

    def _get_default_writer(self) -> WriterFunction:
        return generic_writers.Null()


class Cacher(DataExporter):
    def _get_default_writer(self) -> WriterFunction:
        return generic_writers.Netcdf(monthly_files=True)

    def _get_default_format(self) -> str:
        return ModelFormat.CACHE


class SWAN(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra.SWAN(),
        DnoraDataType.WIND: generic_writers.SWAN(),
        DnoraDataType.WATERLEVEL: generic_writers.SWAN(),
        DnoraDataType.CURRENT: generic_writers.SWAN(),
        DnoraDataType.ICE: generic_writers.SWAN(),
        DnoraDataType.GRID: grid.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra.SWAN(),
        DnoraDataType.WIND: wind.SWAN(),
        DnoraDataType.GRID: grid.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra.WW3(),
        DnoraDataType.WIND: generic_writers.Netcdf(monthly_files=True),
        DnoraDataType.GRID: grid.WW3(),
        DnoraDataType.TRIGRID: grid.WW3Triangular(),
    }

    def _get_default_format(self):
        return ModelFormat.WW3


class HOS_Ocean(DataExporter):
    def _get_default_format(self):
        return ModelFormat.HOS_OCEAN


class REEF3D(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA1D: spectra1d.REEF3D(),
    }

    def _get_default_format(self):
        return ModelFormat.REEF3D
