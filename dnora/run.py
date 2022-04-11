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

class HOS_ocean(ModelExecuter):
    def __init__(self):
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'HOS_ocean'

    def __call__(self, input_file: str, model_folder: str) -> None:
        print('Running HOS_ocean------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        p = Popen(['HOS-ocean'],cwd=model_folder)
        p.wait()


class REEF3D(ModelExecuter):
    def __init__(self):
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'REEF3D'

    #def __call__(self, input_file: str, model_folder: str) -> None:
        #print(input_file, model_folder)
        #model_folder = '/home/konstac/Programs/REEF3D_v22/'
        #print('Copying REEF3D------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        #p = call(['cp ',model_folder+'/REEF3D/bin/REEF3D','REEF3D'],cwd=model_folder)
        #print('Copying DiveMESH------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        #p = call(['cp ',model_folder+'/REEF3D/bin/DiveMESH','DiveMESH'],cwd=model_folder)
        #p.wait()
