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
    def __init__(self, nproc=4):
        self.nproc = nproc
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'SWAN'

    def __call__(self, input_file: str, model_folder: str) -> None:
        print('Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>')
        p = Popen(['swanrun', '-input', input_file,'-omp', str(self.nproc)], cwd=model_folder)
        p.wait()

        return

class SWASH(ModelExecuter):
    def __init__(self):
        pass

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
    def __init__(self, nproc=1):
        self.nproc = nproc
        return

    def _preferred_format(self) -> str:
        """For generation of file name."""
        return 'REEF3D'

    def __call__(self, input_file: str, model_folder: str) -> None:
        p = Popen(['DiveMESH'],cwd=model_folder)
        p.wait()

        if self.nproc == 1:
            p = Popen(['REEF3D'],cwd=model_folder)
        else:
            p = Popen(['/usr/bin/mpirun -n '+str(self.nproc)+' --oversubscribe REEF3D'],cwd=model_folder)
        p.wait()
