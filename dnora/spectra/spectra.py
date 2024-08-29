from copy import copy

# Import objects
from dnora.type_manager.spectral_conventions import SpectralConvention
from .process import boundary_processor_for_convention_change

# Import abstract classes and needed instances of them
from .process import BoundaryProcessor

# Import default values and aux_funcsiliry functions
from dnora import msg

from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import (
    add_time,
    add_frequency,
    add_direction,
    add_datavar,
    add_magnitude,
)

import geo_parameters as gp


@add_datavar(gp.wave.Efth("spec"), coord_group="all", default_value=0.0)
@add_direction(grid_coord=False)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra(PointSkeleton):
    def process_boundary(
        self, boundary_processors: list[BoundaryProcessor] | None = None
    ):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if boundary_processors is None:
            msg.info("No BoundaryProcessor provided. Doing Nothing.")
            return

        if not isinstance(boundary_processors, list):
            boundary_processors = [boundary_processors]

        for processor in boundary_processors:
            msg.process(f"Processing spectra with {type(processor).__name__}")
            print(processor)
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(
                        f"Boundary convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!"
                    )

            new_spec, new_dirs, new_freq, new_inds = processor(
                self.spec(), self.dirs(), self.freq(), self.inds()
            )
            new_inds = list(new_inds)

            if not new_inds:
                self.ds_manager.set_new_ds(None)
                msg.warning(
                    f"No boundary spectra left after processing. Removing all data."
                )
                return

            del_inds = list(set(self.inds()) - set(new_inds))
            if del_inds:
                lon, lat = self.lonlat(native=True)
                msg.info(f"Removing the following points:")
                for ind in list(set(self.inds()) - set(new_inds)):
                    msg.plain(
                        f"ind = {ind}, {self.x_str} = {lon[ind]}, {self.y_str} = {lat[ind]}"
                    )

            metadata = self.meta.get()

            self._init_structure(
                x=self.x(strict=True, inds=new_inds),
                y=self.y(strict=True, inds=new_inds),
                lon=self.lon(strict=True, inds=new_inds),
                lat=self.lat(strict=True, inds=new_inds),
                time=self.time(),
                freq=new_freq,
                dirs=new_dirs,
            )
            self.set_spec(new_spec)
            self.meta.set(metadata)  # Global attributes

            # Set new convention if the processor changed it
            new_convention = processor._convention_out()
            if new_convention is not None:
                self._mark_convention(new_convention)
                msg.info(
                    f"Changing convention from {old_convention} >>> {new_convention}"
                )

        return

    def set_convention(self, convention: SpectralConvention) -> None:
        """Processes boundary to new directional spectral convention and updates metadata."""
        if self.convention() is None:
            msg.info(
                "set_convention changes the convention AND the data. No convention is currently set. Use _mark_convention to define a convention without touching the data."
            )
            return

        if isinstance(convention, str):
            convention = SpectralConvention[convention.upper()]

        boundary_processor = boundary_processor_for_convention_change(
            current_convention=self.convention(), wanted_convention=convention
        )

        if boundary_processor is None:
            msg.info(
                f"Convention ({self.convention()}) already equals wanted convention ({convention})."
            )
            return

        self.process_boundary(boundary_processor)

    def _mark_convention(self, convention: SpectralConvention) -> None:
        """Marks new convention in metadata etc. but does nothing to the spectra"""
        self._convention = convention
        self.meta.append({"spectral_convention": self.convention().value})
        print(f"Spectral convention is now: {self.convention()}")

    def convention(self) -> SpectralConvention | None:
        """Returns the convention (WW3/OCEAN/MET/MATH/MATHVEC) of the spectra"""
        if not hasattr(self, "_convention"):
            return None
        return copy(self._convention)
