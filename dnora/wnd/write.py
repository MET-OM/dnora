import numpy as np
import pandas as pd
from .. import msg

from ..aux import check_if_folder

from .wnd_mod import ForcingWriter # Abstract class
from .wnd_mod import Forcing # Forcing object


class WW3(ForcingWriter):
    def __init__(self, folder: str='', forcing_in_filename: bool=True, time_in_filename: bool=True, grid_in_filename: bool=True) -> None:
        if (not folder == '') and (not folder[-1] == '/'):
            folder = folder + '/'
        self.folder = folder

        self.forcing_in_filename = forcing_in_filename
        self.grid_in_filename = grid_in_filename
        self.time_in_filename = time_in_filename

        return

    def __call__(self, forcing_out: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing_out.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = self.folder + 'wind' + super().create_filename(forcing_out, self.forcing_in_filename, self.grid_in_filename, self.time_in_filename) + '.nc'
        msg.to_file(output_file)

        forcing_out.data.to_netcdf(output_file)

        return


class SWAN(ForcingWriter):
    def __init__(self, folder: str='', forcing_in_filename: bool=True, time_in_filename: bool=True, grid_in_filename: bool=True) -> None:
        if (not folder == '') and (not folder[-1] == '/'):
            folder = folder + '/'
        self.folder = folder

        self.forcing_in_filename = forcing_in_filename
        self.grid_in_filename = grid_in_filename
        self.time_in_filename = time_in_filename

        return

    def __call__(self, forcing_out: Forcing) -> None:
        msg.header(f'{type(self).__name__}: writing wind forcing from {forcing_out.name()}')

        existed = check_if_folder(folder=self.folder, create=True)
        if not existed:
            msg.plain(f"Creating folder {self.folder}")

        output_file = self.folder + 'wind' + super().create_filename(forcing_out, self.forcing_in_filename, self.grid_in_filename, self.time_in_filename) + '.asc'
        msg.to_file(output_file)

        days = forcing_out.days()
        #output_file = f"{forcing_out.grid.name()}_wind{days[0].strftime('%Y%m%d')}_{days[-1].strftime('%Y%m%d')}.asc"

        with open(output_file, 'w') as file_out:
            for day in days:
                msg.plain(day.strftime('%Y-%m-%d'))
                times = forcing_out.times_in_day(day)
                for n in range(len(times)):
                    time_stamp = pd.to_datetime(
                        times[n]).strftime('%Y%m%d.%H%M%S')+'\n'
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing_out.u()
                               [n, :, :]*1000, fmt='%i')
                    file_out.write(time_stamp)
                    np.savetxt(file_out, forcing_out.v()
                               [n, :, :]*1000, fmt='%i')
        return
