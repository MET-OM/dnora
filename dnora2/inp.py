from abc import ABC
from copy import copy
import os
import numpy as np
from dnora2 import msg



class ModelInputFile(ABC):
    def __init__(self):
        pass
    def __call__(self):
        pass
    
    
class SWANInputFile(ModelInputFile):
    def __init__(self, grid):
        self.grid = copy(grid)
        return

    def __call__(self, start_time, end_time, swan_directory = '.', wind =True, calib_wind = 1, calib_wcap = 0.5000E-04):
        path_forcing = os.getcwd() +'/' # path for directory where forcing and boundaries are saved, here it is used the current directory
        DATE_START = start_time.replace('-','').replace('T','.').replace(':','')+'00'
        DATE_END   = end_time.replace('-','').replace('T','.').replace(':','')+'00'
        delta_X = np.round(np.abs(self.grid.lon()[-1] - self.grid.lon()[0]),5)
        delta_Y = np.round(np.abs(self.grid.lat()[-1] - self.grid.lat()[0]),5)
        factor_wind = calib_wind*0.001
        input_file = swan_directory +'/input_'+DATE_START.split('.')[0]+'_'+self.grid.name()+'.swn'
        msg.to_file(input_file)
        with open(input_file, 'w') as file_out:
            file_out.write('$************************HEADING************************\n')
            file_out.write('$ \n')
            file_out.write(' PROJ \'' +self.grid.name()+ '\' \'T24\' \n')
            file_out.write('$ \n')
            file_out.write('$ Topography - Emodnet 2018 \n')
            file_out.write('$ Time of setup: 2021 may \n')
            file_out.write('$ \n')
            file_out.write('$*******************MODEL INPUT*************************\n')
            file_out.write('$ \n')
            file_out.write('SET NAUT \n')
            file_out.write('$ \n')
            file_out.write('MODE NONSTATIONARY TWOD \n')
            file_out.write('COORD SPHE CCM \n')
            file_out.write('CGRID '+str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(delta_X)+' '+str(delta_Y)+' '+str(self.grid.nx()-1)+' '+str(self.grid.ny()-1)+' CIRCLE 36 0.04 1.0 31 \n')
            file_out.write('$ \n')
            file_out.write('INPGRID BOTTOM '  +str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(self.grid.nx()-1)+' '+str(self.grid.ny()-1)+' '+ str((delta_X/(self.grid.nx()-1)).round(4)) +' '+ str((delta_Y/(self.grid.ny()-1)).round(4)) +' EXC 32767\n')
            file_out.write('READINP BOTTOM -1 \''+path_forcing+self.grid.name()+'_SWAN.bot\' 3 0 FREE \n')
            file_out.write('$ \n')
            file_out.write('BOU NEST \''+path_forcing+self.grid.name()+'_spec'+DATE_START.split('.')[0]+'_'+DATE_END.split('.')[0]+'.asc\' OPEN \n')
            file_out.write('$ \n')
            if wind==True:
                file_out.write('INPGRID WIND '+str(self.grid.lon()[0])+' '+str(self.grid.lat()[0])+' 0. '+str(self.grid.nx()-1)+' '+str(self.grid.ny()-1)+' '+str((delta_X/(self.grid.nx()-1)).round(4)) +' '+str((delta_Y/(self.grid.ny()-1)).round(4)) +' NONSTATIONARY '+ DATE_START +' 1 HR ' + DATE_END +'\n')
                file_out.write('READINP WIND '+str(factor_wind)+'  \''+path_forcing+self.grid.name()+'_wind'+DATE_START.split('.')[0]+'_'+DATE_END.split('.')[0]+'.asc\' 3 0 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                file_out.write('OFF QUAD \n')
            file_out.write('GEN3 WESTH cds2='+str(calib_wcap) +'\n')
            file_out.write('FRICTION JON 0.067 \n')
            file_out.write('PROP BSBT \n')
            file_out.write('NUM ACCUR NONST 1 \n')
            file_out.write('$ \n')
            file_out.write('$*******************************************************\n')
            file_out.write('$ Generate block-output \n')
            file_out.write('BLOCK \'COMPGRID\' HEAD \''+path_forcing+self.grid.name()+'_'+DATE_START.split('.')[0]+'.nc\' & \n')
            file_out.write('LAY 1 HSIGN RTP TPS PDIR TM01 DIR DSPR WIND DEP OUTPUT '+ DATE_START +' 1 HR \n')
            file_out.write('$ \n')
            file_out.write('COMPUTE '+DATE_START+' 10 MIN '+ DATE_END + '\n')
            file_out.write('STOP \n')
  
        return input_file
