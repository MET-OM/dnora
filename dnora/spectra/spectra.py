from copy import copy

# Import objects
from dnora.type_manager.spectral_conventions import (
    SpectralConvention,
    spectral_convention_from_string,
)
from dnora.process.spectra import spectral_processor_for_convention_change

from dnora.process.spectra import SpectralProcessor

from dnora import msg

from geo_skeletons import PointSkeleton
from geo_skeletons.decorators import (
    add_time,
    add_frequency,
    add_direction,
    add_datavar,
)

import geo_parameters as gp
from typing import Union, Optional

import numpy as np


@add_datavar(gp.wave.Efth("spec"), coord_group="all", default_value=0.0)
@add_direction(grid_coord=False)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Spectra(PointSkeleton):
    def process(self, spectral_processors: Optional[list[SpectralProcessor]] = None):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if spectral_processors is None:
            msg.info("No BoundaryProcessor provided. Doing Nothing.")
            return

        if not isinstance(spectral_processors, list):
            spectral_processors = [spectral_processors]

        for processor in spectral_processors:
            msg.process(f"Processing spectra with {type(processor).__name__}")
            print(processor)
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(
                        f"Spectral convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!"
                    )

            new_spec, new_dirs, new_freq, new_inds, new_times = processor(
                self.spec(squeeze=False),
                self.dirs(),
                self.freq(),
                self.inds(),
                self.time(),
            )
            new_inds = list(new_inds)

            if not new_inds:
                self._ds_manager.set_new_ds(None)
                msg.warning(f"No spectra left after processing. Removing all data.")
                return

            del_inds = list(set(self.inds()) - set(new_inds))
            if del_inds:
                lon, lat = self.lonlat(native=True)
                if len(del_inds) <= 10:
                    msg.info(f"Removing the following points:")
                    for ind in list(set(self.inds()) - set(new_inds)):
                        msg.plain(
                            f"ind = {ind}, {self.core.x_str} = {lon[ind]}, {self.core.y_str} = {lat[ind]}"
                        )
                else:
                    msg.info(f"Removing a total of {len(del_inds)} points...")

            if len(new_times) < len(self.time()):
                msg.info(
                    f"Removing a total of {len(self.time())-len(new_times)} times..."
                )

            metadata = self.meta.get()

            self._init_structure(
                x=self.x(strict=True, inds=new_inds),
                y=self.y(strict=True, inds=new_inds),
                lon=self.lon(strict=True, inds=new_inds),
                lat=self.lat(strict=True, inds=new_inds),
                time=new_times,
                freq=new_freq,
                dirs=new_dirs,
            )
            self.set_spec(new_spec)
            self.meta.set(metadata)  # Global attributes

            # Set new convention if the processor changed it
            new_convention = processor._convention_out()
            if new_convention is not None:
                msg.process(
                    f"Changing convention from {old_convention} >>> {new_convention}"
                )
                self._mark_convention(new_convention)

        return

    def set_convention(self, convention: SpectralConvention) -> None:
        """Processes boundary to new directional spectral convention and updates metadata."""
        if convention is None:
            return

        if self.convention() is None:
            msg.info(
                "set_convention changes the convention AND the data. No convention is currently set. Use _mark_convention to define a convention without touching the data."
            )
            return

        if isinstance(convention, str):
            convention = SpectralConvention[convention.upper()]

        spectral_processor = spectral_processor_for_convention_change(
            current_convention=self.convention(), wanted_convention=convention
        )

        if spectral_processor is None:
            return

        self.process(spectral_processor)

    def _mark_convention(
        self, convention: SpectralConvention, silent: bool = False
    ) -> None:
        """Marks new convention in metadata etc. but does nothing to the spectra"""
        directions_consistent_with_convention(self.dirs(), convention)
        self._convention = convention
        self.meta.append({"dnora_spectral_convention": self.convention().value})
        if not silent:
            msg.plain(f"Spectral convention is now: {self.convention()}")

    def convention(self) -> Union[SpectralConvention, None]:
        """Returns the convention (WW3/OCEAN/MET/MATH/MATHVEC) of the spectra"""
        if not hasattr(self, "_convention"):
            return None
        return copy(self._convention)


def directions_consistent_with_convention(dirs, convention) -> bool:
    convention = spectral_convention_from_string(convention)
    if convention in [SpectralConvention.UNDEFINED]:
        raise Warning(
            f"Marking convention {convention}, so cannot keep track of wave directionality! {dirs}"
        )
    if convention in [SpectralConvention.WW3, SpectralConvention.MATHVEC]:
        if np.all(np.diff(dirs) > 0):
            raise Warning(
                f"Marking convention {convention}, but directional vector is monotonically increasing! ({dirs})"
            )

    else:
        if not np.all(np.diff(dirs) > 0):
            raise Warning(
                f"Marking convention {convention}, but directional vector is NOT monotonically increasing! ({dirs})"
            )
