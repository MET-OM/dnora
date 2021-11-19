from abc import ABC
from copy import copy
import os
import numpy as np
from . import msg
import pandas as pd

from .wnd.wnd_mod import Forcing # Forcing object
from .grd.grd_mod import Grid # Grid object
from .bnd.bnd_mod import Boundary # Boundary object

from .defaults import dflt_frc, dflt_bnd, dflt_grd


from .aux import create_filename_obj, create_filename_time, add_folder_to_filename


def set_run_times(start_time, end_time, forcing, boundary):
    if start_time is not None:
        DATE_START = pd.Timestamp(start_time)
    elif forcing is None and boundary is None:
        raise Exception('Provide start_time or at least one object (Forcing/Boundary) that has a start time')
    else:
        if forcing is not None:
            frc_start = forcing.time()[0]
        else:
            frc_start = pd.Timestamp('1970-01-01 00:00:00')

        if boundary is not None:
            bnd_start = boundary.time()[0]
        else:
            bdn_start = pd.Timestamp('1970-01-01 00:00:00')

        DATE_START = max(bnd_start, frc_start)

    if end_time is not None:
        DATE_END = pd.Timestamp(end_time)
    elif forcing is None and boundary is None:
        raise Exception('Provide end_time or at least one object (Forcing/Boundary) that has a end time')
    else:
        if forcing is not None:
            frc_end = forcing.time()[-1]
        else:
            frc_end = pd.Timestamp('2300-01-01 00:00:00')

        if boundary is not None:
            bnd_end = boundary.time()[-1]
        else:
            bdn_end = pd.Timestamp('2300-01-01 00:00:00')

        DATE_END = min(bnd_end, frc_end)

    return DATE_START, DATE_END

class ModelInputFile(ABC):
    def __init__(self):
        pass

    def __call__(self):
        pass

class SWANInputFile(ModelInputFile):
    def __init__(self, grid: Grid, forcing: Forcing=None, boundary: Boundary=None,
                grd_folder: str=None,
                frc_filestring: str=None, frc_datestring: str=None,
                frc_folder: str=None,
                bnd_filestring: str=None, bnd_datestring: str=None,
                bnd_folder: str=None):

        self.grid = grid
        self.forcing = forcing
        self.boundary = boundary

        if Forcing is not None:
            if frc_filestring is not None:
                if frc_datestring is None:
                    frc_datestring = dflt_frc['ds']['SWAN']
                # User provided filename overrides all else
                self.frc_filename = forcing.filename(filestring=frc_filestring, datestring=frc_datestring, extension=dflt_frc['ext']['SWAN'])
            else:
                # Use the file the object has been written to, and fall back on SWAN defaults if this information doesn't exist
                self.frc_filename = forcing.written_as(defaults='SWAN')

            if frc_folder is not None:
                self.frc_folder = frc_folder
            else:
                self.frc_folder = forcing.written_to(folder=frc_folder)

        if Boundary is not None:
            if bnd_filestring is not None:
                if bnd_datestring is None:
                    bnd_datestring = dflt_bnd['ds']['SWAN']
                # User provided filename overrides all else
                self.bnd_filename = boundary.filename(filestring=bnd_filestring, datestring=bnd_datestring, extension=dflt_bnd['ext']['SWAN'])
            else:
                # Use the file the object has been written to, and fall back on SWAN defaults if this information doesn't exist
                self.bnd_filename = boundary.written_as(defaults='SWAN')

            if bnd_folder is not None:
                self.bnd_folder = bnd_folder
            else:
                self.bnd_folder = boundary.written_to(folder=bnd_folder)

        if grd_folder is not None:
            self.grd_folder = grd_folder
        else:
            self.grd_folder = grid.written_to(folder=grd_folder)

        self.grd_filename = grid.written_as(filestring=dflt_grd['fs']['SWAN'], extension=dflt_grd['ext']['SWAN'])


        return

    def __call__(self, start_time=None, end_time=None, folder='', filestring='input_$T0_$Grid.swn', datestring='%Y%m%d', calib_wind=1, calib_wcap=0.5000E-04):
        # path for directory where forcing and boundaries are saved, here it is used the current directory
        #path_forcing = forcing_folder
        #path_forcing = os.getcwd() + '/'

        # Define start and end times of model run
        DATE_START, DATE_END = set_run_times(start_time, end_time, self.forcing, self.boundary)

        STR_START = DATE_START.strftime('%Y%m%d.%M%H%S')
        STR_END = DATE_END.strftime('%Y%m%d.%M%H%S')

        #DATE_START = start_time.replace(
        #    '-', '').replace('T', '.').replace(':', '')+'00'
        #DATE_END = end_time.replace(
        #    '-', '').replace('T', '.').replace(':', '')+'00'
        delta_X = np.round(np.abs(self.grid.lon()[-1] - self.grid.lon()[0]), 5)
        delta_Y = np.round(np.abs(self.grid.lat()[-1] - self.grid.lat()[0]), 5)
        factor_wind = calib_wind*0.001

        # Create input file name
        input_file = create_filename_obj(filestring, objects=[self.grid])
        input_file = create_filename_time(input_file, times=[DATE_START])
        input_file = add_folder_to_filename(input_file, folder)
        #input_file = swan_directory + '/input_' + \
        #    DATE_START.split('.')[0]+'_'+self.grid.name()+'.swn'
        msg.to_file(input_file)
        with open(input_file, 'w') as file_out:
            file_out.write(
                '$************************HEADING************************\n')
            file_out.write('$ \n')
            file_out.write(' PROJ \'' + self.grid.name() + '\' \'T24\' \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************MODEL INPUT*************************\n')
            file_out.write('$ \n')
            file_out.write('SET NAUT \n')
            file_out.write('$ \n')
            file_out.write('MODE NONSTATIONARY TWOD \n')
            file_out.write('COORD SPHE CCM \n')
            file_out.write('CGRID '+str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(delta_X)+' '+str(
                delta_Y)+' '+str(self.grid.nx()-1)+' '+str(self.grid.ny()-1)+' CIRCLE 36 0.04 1.0 31 \n')
            file_out.write('$ \n')

            file_out.write('INPGRID BOTTOM ' + str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(self.grid.nx()-1)+' '+str(
                self.grid.ny()-1)+' ' + str((delta_X/(self.grid.nx()-1)).round(4)) + ' ' + str((delta_Y/(self.grid.ny()-1)).round(4)) + ' EXC 32767\n')
            file_out.write('READINP BOTTOM -1 \''+add_folder_to_filename(self.grd_filename, self.grd_folder) +'\' 3 0 FREE \n')
            file_out.write('$ \n')
            file_out.write('BOU NEST \''+add_folder_to_filename(self.bnd_filename, self.bnd_folder)+'\' OPEN \n')
            file_out.write('$ \n')

            if self.forcing is not None:
                file_out.write('INPGRID WIND '+str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(self.forcing.nx()-1)+' '+str(self.forcing.ny()-1)+' '+str(
                    (delta_X/(self.forcing.nx()-1)).round(4)) + ' '+str((delta_Y/(self.forcing.ny()-1)).round(4)) + ' NONSTATIONARY ' + STR_START + ' 1 HR ' + STR_END + '\n')
                file_out.write('READINP WIND '+str(factor_wind)+'  \''+add_folder_to_filename(self.frc_filename, self.frc_folder)+'\' 3 0 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                file_out.write('OFF QUAD \n')
            file_out.write('GEN3 WESTH cds2='+str(calib_wcap) + '\n')
            file_out.write('FRICTION JON 0.067 \n')
            file_out.write('PROP BSBT \n')
            file_out.write('NUM ACCUR NONST 1 \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************************************************\n')
            file_out.write('$ Generate block-output \n')
            file_out.write('BLOCK \'COMPGRID\' HEAD \''+add_folder_to_filename(self.grid.name()+'_'+STR_START.split('.')[0]+'.nc',self.frc_folder)
                           + '\' & \n')
            file_out.write(
                'LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT ' + STR_START + ' 1 HR \n')
            file_out.write('$ \n')
            file_out.write('COMPUTE '+STR_START+' 10 MIN ' + STR_END + '\n')
            file_out.write('STOP \n')

        return input_file
