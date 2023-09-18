from __future__ import annotations # For TYPE_CHECKING
from copy import copy

from typing import List
import numpy as np
import pandas as pd
# Import objects
from .conventions import SpectralConvention
from .process import boundary_processor_for_convention_change
# Import abstract classes and needed instances of them
from .process import BoundaryProcessor
from typing import TYPE_CHECKING

# Import default values and aux_funcsiliry functions
from .. import msg

from skeletons.point_skeleton import PointSkeleton
from skeletons.coordinate_factory import add_time, add_frequency, add_direction
from skeletons.mask_factory import add_mask
from skeletons.datavar_factory import add_datavar

#@add_mask(name='bad', coords='all', default_value=0)
@add_datavar(name='spec', coords='all', default_value=0.)
@add_direction(grid_coord=False)
@add_frequency(grid_coord=False)
@add_time(grid_coord=True)
class Boundary(PointSkeleton):
    def __init__(self, x=None, y=None, lon=None, lat=None, time=pd.date_range('1990-01-01 00:00', '1990-01-01 01:00', freq='H'), freq=np.linspace(0.1,1,10), dirs=np.linspace(0,350,36), name='LonelyBoundary', **kwargs):
        if np.all([a is None for a in [x,y,lon,lat]]):
            x, y = 0, 0
        super().__init__(x=x, y=y, lon=lon, lat=lat, name=name, time=time, freq=freq, dirs=dirs, **kwargs)

    def process_boundary(self, boundary_processors: List[BoundaryProcessor]=None):
        """Process all the individual spectra of the boundary object.

        E.g. change convention form WW3 to Oceanic, interpolate spectra to
        new directional grid, or multiply everything with a constant.
        """

        if boundary_processors is None:
            msg.info("No BoundaryProcessor provided. Doing Nothing.")
            return

        if not isinstance(boundary_processors, list):
            boundary_processors = [boundary_processors]

        convention_warning = False

        for processor in boundary_processors:

            msg.process(f"Processing spectra with {type(processor).__name__}")
            print(processor)
            #self._history.append(copy(processor))
            old_convention = processor._convention_in()
            if old_convention is not None:
                if old_convention != self.convention():
                    msg.warning(f"Boundary convention ({self.convention()}) doesn't match that expected by the processor ({old_convention})!")
                    convention_warning=True


            new_spec, new_dirs, new_freq, new_inds = processor(self.spec(), self.dirs(), self.freq(), self.inds())
            new_inds = list(new_inds)

            if new_inds:
                del_inds = list(set(self.inds())-set(new_inds))
                if del_inds:
                    lon, lat = self.lonlat(native=True)
                    msg.info(f"Removing the following points:")
                    for ind in list(set(self.inds())-set(new_inds)):
                        msg.plain(f"ind = {ind}, {self.x_str} = {lon[ind]}, {self.y_str} = {lat[ind]}")

                metadata = self.metadata()

                self._init_structure(x=self.x(strict=True, inds=new_inds), y=self.y(strict=True, inds=new_inds),
                                lon=self.lon(strict=True, inds=new_inds), lat=self.lat(strict=True, inds=new_inds),
                                time=self.time(), freq=new_freq, dirs=new_dirs)
                self.set_spec(new_spec)
                self.set_metadata(metadata) # Global attributes

                # Set new convention if the processor changed it
                new_convention = processor._convention_out()
                if new_convention is not None:
                    self._set_convention(new_convention, process=False)
                    if convention_warning:
                        msg.warning(f"Convention variable set to {new_convention}, but this might be wrong...")
                    else:
                        msg.info(f"Changing convention from {old_convention} >>> {new_convention}")
            else:
                self.ds_manager.set_new_ds(None)
                msg.warning(f"No boundary spectra left after processing. Removing all data.")

            msg.blank()
        return

    def _set_convention(self, convention: SpectralConvention, process: bool=True) -> None:
        """Sets a new spectral directional convention. To not touch spectra, use process=False."""
        if isinstance(convention, str):
            convention = SpectralConvention[convention.upper()]

        if convention is None:
            msg.info(f"Non new convention given. Keeping convention as {self.convention()}.")
            return

        boundary_processor = boundary_processor_for_convention_change(
                            current_convention = self.convention(),
                            wanted_convention = convention)

        if boundary_processor is None:
            msg.info(f"Convention ({self.convention()}) already equals wanted convention ({convention}).")
            return

        if process:
            self.process_boundary(boundary_processor)
        else:
            self._convention = convention
            self.set_metadata({'spectral_convention': self.convention().value}, append=True)
            print(f'Spectral convention is now: {self.convention()}')

    def convention(self):
        """Returns the convention (WW3/OCEAN/MET/MATH/MATHVEC) of the spectra"""
        if not hasattr(self, '_convention'):
            return None
        return copy(self._convention)
