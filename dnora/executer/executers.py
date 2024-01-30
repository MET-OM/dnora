from .executer import ModelExecuter
from . import run
from . import inputfile
from dnora.model_formats import ModelFormat


class SWAN(ModelExecuter):
    def _get_input_file_writer(self):
        return inputfile.SWAN

    def _get_model_runner(self):
        return run.SWAN

    def _get_default_format(self):
        return ModelFormat.SWAN
