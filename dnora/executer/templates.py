from .executer import ModelExecuter
from . import model_runners
from . import inputfile
from dnora.model_formats import ModelFormat


class SWAN(ModelExecuter):
    def _get_input_file_writer(self):
        return inputfile.SWAN

    def _get_model_runner(self):
        return model_runners.SWAN

    def _get_default_format(self):
        return ModelFormat.SWAN
