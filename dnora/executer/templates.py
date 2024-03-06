from .executer import ModelExecuter
from . import model_runners
from . import inputfile
from dnora.model_formats import ModelFormat
from dnora.dnora_type_manager.dnora_types import DnoraFileType


class SWAN(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWAN}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWAN}

    def _get_default_format(self):
        return ModelFormat.SWAN


class WW3(ModelExecuter):

    _input_file_writers = {
        DnoraFileType.INPUT: inputfile.WW3,
        DnoraFileType.GRID: inputfile.WW3Grid,
        DnoraFileType.WIND: inputfile.WW3Wind,
        DnoraFileType.SPECTRA: inputfile.WW3Spectra,
    }
    _model_runners = {}

    def _get_default_format(self):
        return ModelFormat.WW3
