from __future__ import annotations  # For TYPE_CHECKING


# Import abstract classes and needed instances of them
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from dnora.modelrun.modelrun import ModelRun
    from dnora.file_module import FileNames

from .spectra_writers import SpectraWriter
from dnora.type_manager.dnora_types import DnoraDataType


class REEF3D(SpectraWriter):
    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        **kwargs,
    ) -> tuple[str, str]:
        # Take first spectra and first time step for now
        filename = file_object.get_filepath()
        spectra = model.spectra1d()
        with open(filename, "w") as f:
            spec = spectra.isel(time=0, inds=0).spec(angular=True)
            freq = spectra.freq(angular=True)
            for i, w in enumerate(freq):
                f.write(f"{w:.7f} {spec[i]:.7f}\n")

        return filename
