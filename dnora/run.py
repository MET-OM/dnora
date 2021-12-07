from subprocess import Popen
from abc import ABC, abstractmethod

class ModelExecuter(ABC):
    """Runs the model."""

    @abstractmethod
    def _preferred_format(self) -> str:
        """For the file format using defauts.py, e.g. SWAN or WW3"""
        return

    @abstractmethod
    def __call__(self, input_file: str, model_folder: str) -> None:
        """Runs the model executable"""

        return

class SWAN(ModelExecuter):
    def __init__(self):
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'SWAN'

    def __call__(self, input_file: str, model_folder: str) -> None:
        print('Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        p = Popen(['swanrun', '-input', input_file], cwd=model_folder)
        p.wait()
        return

class SWASH(ModelExecuter):
    def __init__(self):
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'SWASH'

    def __call__(self, input_file: str, model_folder: str) -> None:
        print('Running SWASH----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        p = Popen(['swashrun', '-input', input_file], cwd=model_folder)
        p.wait()
