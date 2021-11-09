import numpy as np
import pandas as pd
from .. import msg

from .wnd_mod import ForcingWriter # Abstract class
from .wnd_mod import Forcing # Forcing object


class WW3(ForcingWriter):
    def __init__(self):
        pass

    def __call__(self, forcing_out: Forcing):
        msg.header(f"Writing output with {type(self).__name__}")
        output_file = f"wind_{forcing_out.name}_{forcing_out.grid.name()}_{str(forcing_out.time()[0])[0:10]}_{str(forcing_out.time()[-1])[0:10]}.nc"
        msg.to_file(output_file)
        forcing_out.data.to_netcdf(output_file)

        return


class SWAN(ForcingWriter):
    def __init__(self):
        pass

    def __call__(self, forcing_out: Forcing):
        msg.header(
            f'{type(self).__name__}: writing wind forcing from {forcing_out.name}')

        days = forcing_out.days()
        output_file = f"{forcing_out.grid.name()}_wind{days[0].strftime('%Y%m%d')}_{days[-1].strftime('%Y%m%d')}.asc"
        msg.info(f'Writing wind forcing to: {output_file}')
        #print(output_file)
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
