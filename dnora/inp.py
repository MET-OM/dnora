from abc import ABC
from copy import copy
import os
import numpy as np
from . import msg
import pandas as pd

from .wnd.wnd_mod import Forcing # Forcing object
from .grd.grd_mod import Grid # Grid object
from .bnd.bnd_mod import Boundary # Boundary object

from .defaults import dflt_frc, dflt_bnd, dflt_grd, dflt_inp, list_of_placeholders


from .aux import add_folder_to_filename, clean_filename

class InputFileWriter(ABC):
    def __init__(self):
        pass

    def _preferred_format(self):
        return 'General'

    def __call__(self):
        pass

class SWAN(InputFileWriter):
    def __init__(self, calib_wind=1, calib_wcap=0.5000E-04, wind=True):

        self.calib_wind = calib_wind
        self.calib_wcap = calib_wcap
        self.wind = wind

        return

    def _preferred_format(self):
        return 'SWAN'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, start_time: str, end_time: str, filename: str, folder: str, grid_path: str, forcing_path: str, boundary_path: str):
        if forcing is None and self.wind == True:
            msg.info('No forcing object provided. Wind information will NOT be written to SWAN input file!')
            self.wind = False



        # Define start and end times of model run
        DATE_START = start_time
        DATE_END = end_time
        STR_START = pd.Timestamp(DATE_START).strftime('%Y%m%d.%H%M%S')
        STR_END = pd.Timestamp(DATE_END).strftime('%Y%m%d.%H%M%S')

        delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 5)
        delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 5)
        factor_wind = self.calib_wind*0.001

        # Create input file name
        output_file = clean_filename(filename, list_of_placeholders)
        output_path = add_folder_to_filename(output_file, folder)

        with open(output_path, 'w') as file_out:
            file_out.write(
                '$************************HEADING************************\n')
            file_out.write('$ \n')
            file_out.write(' PROJ \'' + grid.name() + '\' \'T24\' \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************MODEL INPUT*************************\n')
            file_out.write('$ \n')
            file_out.write('SET NAUT \n')
            file_out.write('$ \n')
            file_out.write('MODE NONSTATIONARY TWOD \n')
            file_out.write('COORD SPHE CCM \n')
            file_out.write('CGRID '+str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(delta_X)+' '+str(
                delta_Y)+' '+str(grid.nx()-1)+' '+str(grid.ny()-1)+' CIRCLE 36 0.04 1.0 31 \n')
            file_out.write('$ \n')

            file_out.write('INPGRID BOTTOM ' + str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(grid.nx()-1)+' '+str(
                grid.ny()-1)+' ' + str((delta_X/(grid.nx()-1)).round(4)) + ' ' + str((delta_Y/(grid.ny()-1)).round(4)) + '\n')
            file_out.write('READINP BOTTOM 1 \''+ grid_path +'\' 3 0 FREE \n')
            file_out.write('$ \n')
            file_out.write('BOU NEST \''+boundary_path+'\' OPEN \n')
            file_out.write('$ \n')

            if self.wind:
                file_out.write('INPGRID WIND '+str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(forcing.nx()-1)+' '+str(forcing.ny()-1)+' '+str(
                    (delta_X/(forcing.nx()-1)).round(4)) + ' '+str((delta_Y/(forcing.ny()-1)).round(4)) + ' NONSTATIONARY ' + STR_START + ' 1 HR ' + STR_END + '\n')
                file_out.write('READINP WIND '+str(factor_wind)+'  \''+forcing_path+'\' 3 0 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                file_out.write('OFF QUAD \n')
            file_out.write('GEN3 WESTH cds2='+str(self.calib_wcap) + '\n')
            file_out.write('FRICTION JON 0.067 \n')
            file_out.write('PROP BSBT \n')
            file_out.write('NUM ACCUR NONST 1 \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************************************************\n')
            file_out.write('$ Generate block-output \n')
            temp_list = forcing_path.split('/')
            forcing_folder = '/'.join(temp_list[0:-1])
            file_out.write('BLOCK \'COMPGRID\' HEAD \''+add_folder_to_filename(grid.name()+'_'+STR_START.split('.')[0]+'.nc',forcing_folder)
                           + '\' & \n')
            file_out.write(
                'LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT ' + STR_START + ' 1 HR \n')
            file_out.write('$ \n')
            file_out.write('COMPUTE '+STR_START+' 10 MIN ' + STR_END + '\n')
            file_out.write('STOP \n')

        return output_file, folder



class SWASH(InputFileWriter):
    def __init__(self, bound_side_command='BOU SIDE W CCW CON REG 0.5 14 270 '):
        self.bound_side_command = bound_side_command

        return

    def _preferred_format(self):
        return 'SWASH'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, start_time: str, end_time: str, filename: str, folder: str, grid_path: str, forcing_path: str, boundary_path: str):

        DATE_START = start_time
        DATE_END = end_time
        STR_START =  pd.Timestamp(DATE_START).strftime('%H%M%S')
        STR_END =  pd.Timestamp(DATE_END).strftime('%H%M%S')

        delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 8)
        delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 8)

        # Create input file name

        output_file = clean_filename(filename, list_of_placeholders)
        output_path = add_folder_to_filename(output_file, folder)

        with open(output_path, 'w') as file_out:
            file_out.write(
                '$************************HEADING************************\n')
            file_out.write('$ \n')
            file_out.write(' PROJ \'' + grid.name() + '\' \'T24\' \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************MODEL INPUT*************************\n')
            file_out.write('$ \n')
            file_out.write('SET NAUT \n')
            file_out.write('$ \n')
            file_out.write('MODE NONSTATIONARY TWOD \n')
            file_out.write('COORD SPHE CCM \n')
            file_out.write('CGRID REG '+str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(delta_X)+' '+str(
                delta_Y)+' '+str(grid.nx()-1)+' '+str(grid.ny()-1)+' \n')
            file_out.write('$ \n')
            file_out.write('VERT 1 \n')
            file_out.write('$ \n')
            file_out.write('INPGRID BOTTOM ' + str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(grid.nx()-1)+' '+str(
                grid.ny()-1)+' ' + str((delta_X/(grid.nx()-1)).round(8)) + ' ' + str((delta_Y/(grid.ny()-1)).round(8)) +  ' EXC -999 \n')
            file_out.write('READINP BOTTOM 1 \''+grid_path +'\' 3 0 FREE \n')
            file_out.write('$ \n')
            file_out.write(self.bound_side_command +' \n')
            #file_out.write('BOU NEST \''+add_folder_to_filename(self.bnd_filename, self.bnd_folder)+'\' OPEN \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************************************************\n')
            file_out.write('$ OUTPUT REQUESTS \n')
            file_out.write('BLOCK \'COMPGRID\' NOHEAD \''+add_folder_to_filename(grid.name()+'.mat',folder)
                           + '\' & \n')
            file_out.write(
                'LAY 3 WATL BOTL OUTPUT ' + STR_START + ' 5 SEC \n')
            file_out.write('$ \n')
            file_out.write('COMPUTE '+STR_START+' 0.001 SEC ' + STR_END + '\n')
            file_out.write('STOP \n')

        return output_file, folder
