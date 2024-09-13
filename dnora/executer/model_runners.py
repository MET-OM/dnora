from subprocess import Popen
from abc import ABC, abstractmethod
from dnora.file_module import FileNames
from dnora.type_manager.model_formats import ModelFormat
from .post_processors import PostProcessor, SwashMatToNc, HosOceanToNc


class ModelRunner(ABC):
    """Runs the model."""

    def __init__(self, model):
        self.model = model

    @abstractmethod
    def _preferred_format(self) -> str:
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

    def _preferred_format(self) -> str:
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

    def _preferred_format(self) -> str:
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


class HOS_ocean(ModelRunner):
    def __init__(self):
        return

    def _preferred_format(self) -> str:
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
