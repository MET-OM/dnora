from subprocess import Popen
from abc import ABC, abstractmethod
from dnora.file_module import FileNames
from dnora.type_manager.model_formats import ModelFormat
from .post_processors import PostProcessor, SwashMatToNc, HosOceanToNc
from dnora import msg
import shutil
class ModelRunner(ABC):
    """Runs the model."""

    def __init__(self, model):
        self.model = model

    @abstractmethod
    def preferred_format(self) -> str:
        """For the file format using defauts.py, e.g. ModelFormat.SWAN"""
        return

    def post_processors(self) -> list[PostProcessor]:
        return []

    @abstractmethod
    def __call__(
        self,
        input_file: str = None,
        model_folder: str = None,
    ) -> None:
        """Runs the model executable"""

        return


class SWAN(ModelRunner):
    def __init__(self):
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.SWAN

    def __call__(self, file_object, nproc=4, **kwargs) -> None:

        print("Running SWAN----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>")
        p = Popen(
            ["swanrun", "-input", file_object.get_filename(), "-omp", f"{nproc}"],
            cwd=file_object.get_folder(),
        )
        p.wait()

        return


class SWASH(ModelRunner):
    def __init__(self):
        pass

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.SWASH

    def post_processors(self) -> list[PostProcessor]:
        return [SwashMatToNc()]

    def __call__(self, file_object: FileNames) -> None:
        print("Running SWASH----------------------->>>>>>>>>>>>>>>>>>>>>>>>>>")
        p = Popen(
            ["swashrun", "-input", file_object.get_filename()],
            cwd=file_object.get_folder(),
        )
        p.wait()


class WW3(ModelRunner):
    def __init__(self, program: str):
        """E.g. program = 'grid' to run ww3_grid etc."""
        self.program = program
        if program == 'shel':
            self._post_processors = [WW3('ounf')]
        else:
            self._post_processors = []
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.WW3

    def post_processors(self) -> list[PostProcessor]:
        return self._post_processors

    def __call__(self, file_object, model_folder, nproc=4, **kwargs) -> None:

        if model_folder:
            from_file = f"{model_folder}/ww3_{self.program}"
            to_file = file_object.get_folder()
            msg.copy_file(from_file,to_file)
            shutil.copy(from_file, to_file)

        filename_out = f'{file_object.get_folder()}/ww3_{self.program}.out'
        msg.info(f"Running ww3_{self.program}...")
        msg.to_file(filename_out)
        with open(filename_out, 'w') as outfile:
            p = Popen(
                [f"ww3_{self.program}"],
                cwd=file_object.get_folder(),
                stdout=outfile,
            )
            p.wait()
        
        
        return




class HOS_ocean(ModelRunner):
    def __init__(self):
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.HOS_OCEAN

    def post_processors(self) -> list[PostProcessor]:
        return [HosOceanToNc()]

    def __call__(self, input_file: str, model_folder: str) -> None:
        print("Running HOS_ocean------------------->>>>>>>>>>>>>>>>>>>>>>>>>>")
        p = Popen(["HOS-ocean"], cwd=model_folder)
        p.wait()


class REEF3D(ModelRunner):
    def __init__(self, nproc=1):
        self.nproc = nproc
        return

    def preferred_format(self) -> str:
        """For generation of file name."""
        return ModelFormat.REEF3D

    def __call__(self, input_file: str, model_folder: str) -> None:
        p = Popen(["DiveMESH"], cwd=model_folder)
        p.wait()

        if self.nproc == 1:
            p = Popen(["REEF3D"], cwd=model_folder)
        else:
            p = Popen(
                ["/usr/bin/mpirun -n " + str(self.nproc) + " --oversubscribe REEF3D"],
                cwd=model_folder,
            )
        p.wait()
