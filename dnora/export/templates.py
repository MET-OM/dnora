from .exporter import DataExporter
from . import data_writers, grid_writers, spectra_writers, spectra1d_writers

from dnora.type_manager.dnora_types import DnoraDataType
from .exporter import WriterFunction
from dnora.type_manager.model_formats import ModelFormat


class NullExporter(DataExporter):
    def _get_default_writer(self) -> WriterFunction:
        return data_writers.Null()


class Netcdf(DataExporter):
    def _get_default_writer(self) -> WriterFunction:
        return data_writers.Netcdf()


class Cacher(DataExporter):
    _silent = True

    def _get_default_writer(self) -> WriterFunction:
        return data_writers.Netcdf(daily_files=True)

    def _get_default_format(self) -> str:
        return ModelFormat.CACHE


class SWAN(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra_writers.SWAN(),
        DnoraDataType.WIND: data_writers.SWAN(),
        DnoraDataType.WATERLEVEL: data_writers.SWAN(),
        DnoraDataType.CURRENT: data_writers.SWAN(),
        DnoraDataType.ICE: data_writers.SWAN(),
        DnoraDataType.GRID: grid_writers.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra_writers.SWAN(),
        DnoraDataType.WIND: data_writers.SWAN(),
        DnoraDataType.GRID: grid_writers.SWAN(),
    }

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA: spectra_writers.WW3(),
        DnoraDataType.WIND: data_writers.Netcdf(monthly_files=True),
        DnoraDataType.GRID: grid_writers.WW3(),
        DnoraDataType.TRIGRID: grid_writers.WW3Triangular(),
    }

    def _get_default_format(self):
        return ModelFormat.WW3


class HOS_Ocean(DataExporter):
    def _get_default_format(self):
        return ModelFormat.HOS_OCEAN


class REEF3D(DataExporter):
    _writer_dict = {
        DnoraDataType.SPECTRA1D: spectra1d_writers.REEF3D(),
    }

    def _get_default_format(self):
        return ModelFormat.REEF3D
