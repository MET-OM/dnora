from .executer import ModelExecuter
from . import model_runners
from . import inputfile
from dnora.type_manager.model_formats import ModelFormat
from dnora.type_manager.dnora_types import DnoraFileType


class SWAN(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWAN()}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWAN()}

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWASH()}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWASH()}

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(ModelExecuter):

    _input_file_writers = {
        DnoraFileType.INPUT: inputfile.WW3(),
        DnoraFileType.GRID: inputfile.WW3Grid(),
        DnoraFileType.WIND: inputfile.WW3Forcing(DnoraFileType.WIND),
        DnoraFileType.CURRENT: inputfile.WW3Forcing(DnoraFileType.CURRENT),
        DnoraFileType.WATERLEVEL: inputfile.WW3Forcing(DnoraFileType.WATERLEVEL),
        DnoraFileType.ICE: inputfile.WW3Forcing(DnoraFileType.ICE),
        DnoraFileType.SPECTRA: inputfile.WW3Spectra(),
    }
    _model_runners = {DnoraFileType.GRID: model_runners.WW3('grid'), DnoraFileType.WIND: model_runners.WW3('prnc'), DnoraFileType.SPECTRA: model_runners.WW3('bounc'), DnoraFileType.INPUT: model_runners.WW3('shel'), DnoraFileType.CURRENT: model_runners.WW3('prnc'), DnoraFileType.WATERLEVEL: model_runners.WW3('prnc'),DnoraFileType.ICE: model_runners.WW3('prnc')}

    def _get_default_format(self):
        return ModelFormat.WW3
