from __future__ import annotations
from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from .. import msg
from copy import copy
from ..aux import check_if_folder, create_filename_obj, create_filename_time, add_file_extension, add_folder_to_filename

#from .wnd_mod import ForcingWriter # Abstract class
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from .wnd_mod import Forcing # Boundary object
from ..defaults import dflt_frc

class ForcingWriter(ABC):
    """Writes the forcing data to a certain file format.

    This object is provided to the .export_forcing() method.
    """

    @abstractmethod
    def __call__(self, forcing_out: Forcing) -> Tuple[str, str]:
        """Writed the data from the Forcing object and returns the file and
        folder where data were written."""

        return output_file, output_folder

    # def create_filename(self, forcing_out: Forcing, forcing_in_filename: bool=True, grid_in_filename: bool=True, time_in_filename: bool=True) -> str:
    #     """Creates a filename based on the boolean swithes set in __init__ and the meta data in the objects"""
    #
    #     forcing_fn = ''
    #     grid_fn = ''
    #     time_fn = ''
    #
    #     if forcing_in_filename:
    #         forcing_fn = f"_{forcing_out.name()}"
    #
    #     if grid_in_filename:
    #         grid_fn = f"_{forcing_out.grid.name()}"
    #
    #     if time_in_filename:
    #         time_fn = f"_{str(forcing_out.time()[0])[0:10]}_{str(forcing_out.time()[-1])[0:10]}"
    #
    #     filename = forcing_fn + grid_fn + time_fn
    #
    #     return filename

class WW3(ForcingWriter):
    """Writes wind forcing data to WAVEWATH III netcdf format."""

    def __init__(self, folder: str=dflt_frc['fldr']['WW3'], filestring: str=dflt_frc['fs']['WW3'], datestring: str=dflt_frc['ds']['WW3']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> Tuple[str, str]:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = forcing.filename(filestring=self.filestring, datestring=self.datestring, extension='nc')

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)

        msg.to_file(output_path)
        forcing.data.to_netcdf(output_path)

        return output_file, self.folder


class SWAN(ForcingWriter):
    """Writes wind forcing data to SWAN ascii format."""

    def __init__(self, folder: str=dflt_frc['fldr']['SWAN'], filestring: str=dflt_frc['fs']['SWAN'], datestring: str=dflt_frc['ds']['SWAN']) -> None:
        self.folder = copy(folder)
        self.filestring = copy(filestring)
        self.datestring = copy(datestring)

        return

    def __call__(self, forcing: Forcing) -> Tuple[str, str]:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = forcing.filename(filestring=self.filestring, datestring=self.datestring, extension='asc')

        # Add folder
        output_path = add_folder_to_filename(output_file, folder=self.folder)
        msg.to_file(output_path)

        days = forcing.days()
        with open(output_path, 'w') as file_out:
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = forcing.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.u()
                               [n, :, :]*1000, fmt='%i')
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing.v()
                               [n, :, :]*1000, fmt='%i')


        return output_file, self.folder
