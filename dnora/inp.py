
from abc import ABC, abstractmethod
from cmath import nan
from copy import copy
import os
import numpy as np
import pandas as pd

# Import objects
from .wnd.wnd_mod import Forcing
from .grd.grd_mod import Grid
from .bnd.bnd_mod import Boundary
from .wlv.wlv_mod import WaterLevel
from .ocr.ocr_mod import OceanCurrent
from .aux_funcs import create_swan_segment_coords
from . import msg
from . import file_module
class InputFileWriter(ABC):
    @abstractmethod
    def _extension(self):
        pass

    def _clean_filename(self):
        """If this is set to False, then the ModelRun object does not clean
        the filename, and possible placeholders (e.g. #T0) can still be
        present.
        """
        return True

    def _im_silent(self) -> bool:
        """Return False if you want to be responsible for printing out the
        file names."""
        return True

    @abstractmethod
    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, oceancurrent: OceanCurrent,
                start_time: str, end_time: str, filename: str,
                grid_path: str, forcing_path: str, boundary_path: str) -> str:
        return output_file

class SWAN(InputFileWriter):
    def __init__(self, calib_wind=1, calib_wcap=0.5000E-04, calib_wlev = 1, calib_ocr = 1, wind=True, waterlevel=True, oceancurrent=True, timestep=10,
                 f_low = 0.04, f_high=1., n_freq=31, n_dir=36, spec_points=None, extension='swn', output_var='HSIGN RTP TPS PDIR TM01 TMM10 DIR DSPR DEP'):

        self.calib_wind = calib_wind
        self.calib_wcap = calib_wcap
        self.calib_wlev = calib_wlev
        self.calib_ocr = calib_ocr
        self.wind = wind
        self.waterlevel = waterlevel
        self.oceancurrent = oceancurrent
        self.spec_points = spec_points # list of (lon, lat) points, e.g.,[(4.4, 60.6),(4.4, 60.8)]
        self._extension_in = extension
        self.swan_timestep = timestep
        self.f_low = f_low
        self.f_high = f_high
        self.n_freq = n_freq
        self.n_dir = n_dir
        self.output_var = output_var
        return

    def _extension(self):
        return self._extension_in

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, waterlevel: WaterLevel, oceancurrent: OceanCurrent,
                start_time: str, end_time: str, filename: str,
                grid_path: str, forcing_path: str, boundary_path: str, waterlevel_path: str, oceancurrent_path: str):

        if forcing is None and self.wind == True:
            msg.info('No forcing object provided. Wind information will NOT be written to SWAN input file!')
            self.wind = False

        if waterlevel is None and self.waterlevel == True:
            msg.info('No waterlevel object provided. Waterlevel information will NOT be written to SWAN input file!')
            self.waterlevel = False

        if oceancurrent is None and self.oceancurrent == True:
            msg.info('No oceancurrent object provided. OceanCurrent information will NOT be written to SWAN input file!')
            self.oceancurrent = False

        # Define start and end times of model run
        DATE_START = start_time
        DATE_END = end_time
        STR_START = pd.Timestamp(DATE_START).strftime('%Y%m%d.%H%M%S')
        STR_END = pd.Timestamp(DATE_END).strftime('%Y%m%d.%H%M%S')
        HOTSTART_FILE = 'hotstart_'+grid.name()+'_'+(pd.Timestamp(DATE_START)-pd.Timedelta(hours=1)).strftime('%Y%m%d%H%M')
        # STR_FORCING_START = STR_START
        # STR_FORCING_END = STR_END

        # # For wind forcing, range of all time steps in .asc file needs to be specified
        # if forcing is not None:
        #     if (forcing.name() == "ERA5"):  # bugfix: for era5 complete 24 hour timesteps are written to asc file.
        #         STR_FORCING_START = pd.Timestamp(DATE_START).strftime('%Y%m%d') + '.000000'
        #         STR_FORCING_END = pd.Timestamp(DATE_END).strftime('%Y%m%d') + '.230000'

        delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 5)
        delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 5)

        factor_wind = self.calib_wind*0.001
        factor_waterlevel = self.calib_wlev*0.001
        factor_oceancurrent = self.calib_ocr*0.001

        with open(filename, 'w') as file_out:
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
                delta_Y)+' '+str(grid.nx()-1)+' '+str(grid.ny()-1)+' CIRCLE %d %f %f %d \n' %(self.n_dir, self.f_low,
                                                                                              self.f_high, self.n_freq))
            file_out.write('$ \n')

            file_out.write('INPGRID BOTTOM ' + str(grid.lon()[0])+' '+str(grid.lat()[0])+' 0. '+str(grid.nx()-1)+' '+str(
                grid.ny()-1)+' ' + str((delta_X/(grid.nx()-1)).round(8)) + ' ' + str((delta_Y/(grid.ny()-1)).round(8)) + '\n')
            file_out.write('READINP BOTTOM 1 \''+ grid_path.split('/')[-1] +'\' 3 0 FREE \n')
            file_out.write('$ \n')
            if boundary is not None:
                lons, lats = create_swan_segment_coords(grid.boundary_mask(), grid.lon_edges(), grid.lat_edges())

            if boundary is not None:
                lons, lats = create_swan_segment_coords(grid.boundary_mask(), grid.lon_edges(), grid.lat_edges())

                bound_string = "BOUNDSPEC SEGMENT XY"
                for lon, lat in zip(lons, lats):
                    bound_string += f" {lon:.4f} {lat:.4f}"
                bound_string += " VARIABLE FILE 0 "
                bound_string += f"'{boundary_path.split('/')[-1]}'\n"
                file_out.write(bound_string)

                #file_out.write('BOU NEST \''+boundary_path.split('/')[-1]+'\' OPEN \n')
                file_out.write('$ \n')

            if self.wind:
                self.output_var = self.output_var + ' WIND'
                delta_Xf = np.round(np.abs(forcing.lon()[-1] - forcing.lon()[0]), 5)
                delta_Yf = np.round(np.abs(forcing.lat()[-1] - forcing.lat()[0]), 5)

                file_out.write('INPGRID WIND ' + str(forcing.lon()[0].round(3)) + ' ' + str(forcing.lat()[0].round(3)) + ' 0. ' + str(
                    forcing.nx() - 1) + ' ' + str(forcing.ny() - 1) + ' ' + str(
                    (delta_Xf / (forcing.nx() - 1)).round(6)) + ' ' + str((delta_Yf / (forcing.ny() - 1)).round(
                    6)) + ' NONSTATIONARY ' + STR_START + f" {forcing.dt():.0f} HR " + STR_END + '\n')

                file_out.write('READINP WIND '+str(factor_wind)+'  \''+forcing_path.split('/')[-1]+'\' 3 0 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                file_out.write('WIND 0 0 \n') # no wind forcing

            if self.waterlevel:
                self.output_var = self.output_var + ' WATLEV'
                delta_Xf = np.round(np.abs(waterlevel.lon()[-1] - waterlevel.lon()[0]), 5)
                delta_Yf = np.round(np.abs(waterlevel.lat()[-1] - waterlevel.lat()[0]), 5)

                file_out.write('INPGRID WLEV ' + str(waterlevel.lon()[0]) + ' ' + str(waterlevel.lat()[0]) + ' 0. ' + str(
                    waterlevel.nx() - 1) + ' ' + str(waterlevel.ny() - 1) + ' ' + str(
                    (delta_Xf / (waterlevel.nx() - 1)).round(6)) + ' ' + str((delta_Yf / (waterlevel.ny() - 1)).round(
                    6)) + ' NONSTATIONARY ' + STR_START + f" {waterlevel.dt():.0f} HR " + STR_END + '\n')

                file_out.write('READINP WLEV '+str(factor_waterlevel)+'  \''+waterlevel_path.split('/')[-1]+'\' 3 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                pass

            if self.oceancurrent:
                self.output_var = self.output_var + ' VEL'
                delta_Xf = np.round(np.abs(oceancurrent.lon()[-1] - oceancurrent.lon()[0]), 5)
                delta_Yf = np.round(np.abs(oceancurrent.lat()[-1] - oceancurrent.lat()[0]), 5)

                file_out.write('INPGRID CUR ' + str(oceancurrent.lon()[0]) + ' ' + str(oceancurrent.lat()[0]) + ' 0. ' + str(
                    oceancurrent.nx() - 1) + ' ' + str(oceancurrent.ny() - 1) + ' ' + str(
                    (delta_Xf / (oceancurrent.nx() - 1)).round(6)) + ' ' + str((delta_Yf / (oceancurrent.ny() - 1)).round(
                    6)) + ' EXC 32767 NONSTATIONARY ' + STR_START + f" {oceancurrent.dt():.0f} HR " + STR_END + ' \n')

                file_out.write('READINP CUR '+str(factor_oceancurrent)+'  \''+oceancurrent_path.split('/')[-1]+'\' 3 0 0 1 FREE \n')
                file_out.write('$ \n')
            else:
                pass

            if os.path.isfile(grid_path.split('/')[0]+'/'+HOTSTART_FILE) is True:
                file_out.write('INITIAL HOTSTART \''+HOTSTART_FILE+'\''  '\n')


            file_out.write('GEN3 WESTH cds2='+str(self.calib_wcap) +' AGROW'+ '\n')
            file_out.write('FRICTION JON 0.067 \n')
            file_out.write('PROP BSBT \n')
            file_out.write('NUM ACCUR NONST 1 \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************************************************\n')

            file_out.write('$ Generate block-output \n')
            temp_list = forcing_path.split('/')
            forcing_folder = '/'.join(temp_list[0:-1])
            file_out.write('BLOCK \'COMPGRID\' HEAD \''+grid.name()+'_'+STR_START.split('.')[0]+'.nc'
                           + '\' & \n')
            file_out.write(
                'LAY 1 '+self.output_var+' OUTPUT ' + STR_START + ' 1 HR \n')
            file_out.write('$ \n')
            if self.spec_points is  not None:
                file_out.write('POINTS \'pkt\' &\n')
                for i in range(len(self.spec_points)):
                    file_out.write(str(self.spec_points[i][0])+' '+str(self.spec_points[i][1])+ ' &\n')
                file_out.write('SPECOUT \'pkt\' SPEC2D ABS \''+grid.name()+'_'+STR_START.split('.')[0]+'_spec.nc'+ '\' & \n')
                file_out.write('OUTPUT ' + STR_START + ' 1 HR \n')
            else:
                pass
            file_out.write('COMPUTE '+STR_START+ ' %d MIN ' % self.swan_timestep + STR_END + '\n')
            file_out.write('HOTFILE \'hotstart_'+grid.name()+'_'+STR_END.replace('.','')[:-2]+'\'' +' FREE \n')
            file_out.write('STOP \n')

        return filename



class SWASH(InputFileWriter):
    def __init__(self,bound_side_command='BOU SIDE W CCW CON REG 0.5 14 270 '):
        self.bound_side_command = bound_side_command



        return

    def _extension(self):
        return 'sws'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, waterlevel: WaterLevel, oceancurrent: OceanCurrent,
                start_time: str, end_time: str, filename: str,
                grid_path: str, forcing_path: str, boundary_path: str, waterlevel_path: str, oceancurrent_path: str):

        DATE_START = start_time
        DATE_END = end_time
        STR_START =  pd.Timestamp(DATE_START).strftime('%H%M%S')
        STR_END =  pd.Timestamp(DATE_END).strftime('%H%M%S')

        delta_X = np.round(np.abs(grid.lon()[-1] - grid.lon()[0]), 8)
        delta_Y = np.round(np.abs(grid.lat()[-1] - grid.lat()[0]), 8)

        with open(filename, 'w') as file_out:
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
            file_out.write('READINP BOTTOM 1 \''+grid_path.split('/')[-1] +'\' 3 0 FREE \n')
            file_out.write('$ \n')
            file_out.write(self.bound_side_command +' \n')
            #file_out.write('BOU NEST \''+add_folder_to_filename(self.bnd_filename, self.bnd_folder)+'\' OPEN \n')
            file_out.write('$ \n')
            file_out.write(
                '$*******************************************************\n')
            file_out.write('$ OUTPUT REQUESTS \n')
            temp_list = grid_path.split('/')
            forcing_folder = '/'.join(temp_list[0:-1])
            file_out.write('BLOCK \'COMPGRID\' NOHEAD \''+grid.name()+'.mat'
                           + '\' & \n')
            file_out.write(
                'LAY 3 WATL BOTL OUTPUT ' + STR_START + ' 5 SEC \n')
            file_out.write('$ \n')
            file_out.write('COMPUTE '+STR_START+' 0.001 SEC ' + STR_END + '\n')
            file_out.write('STOP \n')

        return filename

class REEF3D(InputFileWriter):
    def __init__(self, option = 'REEF3D', edges = ['W'], nproc = 1,rot_angle=0, wave_input='SPEC1D', Hs = nan, Tp = nan , Sprm = nan):
        self.option = option
        self.nproc = nproc # number of processors
        self.edges = edges
        self.rot_angle = rot_angle # anlge for domain rotation
        self.wave_input = wave_input # 'SPEC1D' for NORA3 frequency spectrum or 'JONSWAP' providing Hs and Tp
        self.Hs = Hs # valid only for wave_input = 'JONSWAP'
        self.Tp = Tp # valid only for wave_input = 'JONSWAP'
        self.Sprm = Sprm # valid only for wave_input = 'JONSWAP'
        return

    def _extension(self):
        return 'txt'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary, waterlevel: WaterLevel, oceancurrent: OceanCurrent,
                start_time: str, end_time: str, filename: str, forcing_path: str,waterlevel_path: str, oceancurrent_path: str,
                grid_path: str, boundary_path: str):

        geodat = pd.read_csv(grid_path, sep = ' ') # read geo.dat
        geodat.columns = ['x', 'y', 'z']

        if self.option == 'DiveMESH':
            filename =  '/'.join(filename.split('/')[:-1])+'/control.txt'
            with open(filename, 'w') as file_out:
                if 'W' in self.edges:
                    file_out.write('C 11 6 // West side: wave generation' '\n')
                    file_out.write('C 12 7 // side: numerical beach' '\n')
                    file_out.write('C 13 7 // side: numerical beach' '\n')
                    file_out.write('C 14 7 // side: numerical beach' '\n')
                elif 'N' in self.edges:
                    file_out.write('C 12 6 // North side: wave generation' '\n')
                    file_out.write('C 13 3 // side: numerical beach' '\n')
                    file_out.write('C 14 3 // side: numerical beach' '\n')
                    file_out.write('C 11 7 // side: numerical beach' '\n')
                elif 'E' in self.edges:
                    file_out.write('C 14 6 // East side: wave generation' '\n')
                    file_out.write('C 11 7 // side: numerical beach' '\n')
                    file_out.write('C 12 7 // side: numerical beach' '\n')
                    file_out.write('C 13 7 // side: numerical beach' '\n')
                elif 'S' in self.edges:
                    file_out.write('C 13 6 // South side: wave generation' '\n')
                    file_out.write('C 11 7 // side: numerical beach' '\n')
                    file_out.write('C 12 7 // side: numerical beach' '\n')
                    file_out.write('C 14 7 // side: numerical beach' '\n')

                file_out.write('C 15 21 // bottom: wall boundary' '\n')
                file_out.write('C 16 3 // top: symmetry plane' '\n')
                file_out.write(' \n')

                file_out.write('B 1 '+str(int(grid.dx().round(0)))+' // horizontal mesh size dx' '\n')
                file_out.write('B 2 '+str(int(geodat.x.max()/int(grid.dx().round(0))))+
                ' '+str(int(geodat.y.max()/int(grid.dx().round(0))))
                +' 10 // number of cells in x, y and z directions' '\n')
                file_out.write('B 10 0.0 '+str(int(geodat.x.max()))+' 0.0 '+str(int(geodat.y.max()))+' 0.0 1.0 // rectangular domain size' '\n')
                file_out.write(' \n')

                file_out.write('B 103 5 // vertical grid clustering' '\n')
                file_out.write('B 113 2.5 // the stretching factor for the vertical grid clustering' '\n')
                file_out.write('B 116 1.0 // the focal point for the vertical grid clustering, which is water depth here' '\n')
                file_out.write(' \n')

                file_out.write('G 10 1 // turn geodat on/off' '\n')
                file_out.write('G 13 '+str(self.rot_angle)+' // rotation angle of geo coordinates around vertical axis' '\n')
                file_out.write('G 15 2 // local inverse distance interpolation' '\n')
                file_out.write('G 20 0 // use automatic grid size off' '\n')
                file_out.write('G 31 14 // number of smoothing iterations' '\n')
                #file_out.write('G 41 1' '\n')
                file_out.write(' \n')

                file_out.write('M 10 '+str(self.nproc)+' // number of processors' '\n')
                file_out.write('M 20 2 // decomposition method 2' '\n')
        elif self.option == 'REEF3D':
            with open(filename, 'w') as file_out:
                file_out.write('A 10 3  // choose the model reef::fnpf' '\n')
                file_out.write('A 310 3 // 3rd-order runge-kutta for fsfbc time treatment' '\n')
                file_out.write('A 311 5 // 5th-order weno for fsfbc spatial treatment including wetting-drying' '\n')
                file_out.write('A 320 1 // 2nd-order laplace' '\n')
                file_out.write(' \n')

                file_out.write('A 341 2.0 // size of coastal relaxation zone by a factor of the horizontal cell size' '\n')
                file_out.write('A 343 1   // turn on wetting-drying' '\n')
                file_out.write('A 345 0.001 // wetting-drying water depth threshold' '\n')
                file_out.write('A 346 2.1   // added viscosity within the coastal relaxation zone' '\n')
                file_out.write(' \n')

                file_out.write('A 350 1 // viscosity damping wave breaking algorithm' '\n')
                file_out.write('A 351 3 // breaking wave detection for both deep and shallow water' '\n')
                file_out.write('A 352 3 // additional filtering for viscosity based breaking for both deep and shallow water' '\n')
                file_out.write('A 361 5 // filtering outer iterations' '\n')
                file_out.write('A 362 2 // filtering inner iterations' '\n')
                file_out.write('A 365 1.86 // artificial viscosity for breaking wave energy dissipation' '\n')
                file_out.write(' \n')
                if self.wave_input =='SPEC1D':
                    file_out.write('B 85 10 // spectrum file' '\n')
                    file_out.write('B 90 1 // wave input' '\n')
                    file_out.write('B 92 31 // 1st-order irregular wave' '\n')

                elif self.wave_input =='JONSWAP':
                    file_out.write('B 85 2 // jonswap' '\n')
                    file_out.write('B 90 1 // wave input' '\n')
                    file_out.write('B 92 31 // 1st-order irregular wave' '\n')
                    file_out.write('B 93 '+str(self.Hs)+' '+ str(self.Tp)+' // wave height, wave period' '\n')
                    file_out.write('B 134 '+str(self.Sprm)+' // spreading parameter for the directional spreading functions' '\n')


                file_out.write(' \n')
                file_out.write('B 96 200.0 400.0 // wave generation zone length and numerical beach length' '\n')
                #file_out.write('B 107 0.0 '+str(int(geodat.x.max()))+' 0.0 0.0 200.0 // wave generation zone length and numerical beach length' '\n')
                #file_out.write('B 107 0.0 '+str(int(geodat.x.max()))+' '+str(int(geodat.y.max()))+' '+str(int(geodat.x.max()))+' 200.0 // customised numerical beach at the side walls' '\n')
                #file_out.write('B 107 25000.0 12000.0 0.0 16000.0 200.0 // customised numerical beach at the end of the tank' '\n')
                #file_out.write('B 107 0.0 0.0 2900.0 3500.0 200.0 // customised numerical beach at the side walls' '\n')
                #file_out.write('B 108 0.0 0.0 0.0 '+str(int(geodat.y.max()))+' 200.0 // customised wave generation zone' '\n')
                file_out.write('B 98 2 // relaxation method 2 for wave generation' '\n')
                file_out.write('B 99 2 // relaxation method 2 for numerical beach' '\n')
                file_out.write(' \n')

                file_out.write('F 60 '+str(grid.data.topo.max().round(1).values)+' // still water depth' '\n')
                file_out.write(' \n')

                file_out.write('G 50 1 // read in geo bathymetry' '\n')
                file_out.write(' \n')

                file_out.write('I 30 0 // turn off full tank initialisation, one can turn it on for a quick check of the setup' '\n')
                file_out.write(' \n')

                file_out.write('N 41 1800.0 // simulation time' '\n')
                file_out.write('N 47 1.0 // cfl number' '\n')
                file_out.write(' \n')

                file_out.write('M 10 '+str(self.nproc)+' // number of processors' '\n')
                file_out.write(' \n')

                file_out.write('P 180 1 // turn on .vtp free surface printout' '\n')
                file_out.write('P 185 0.0 1800.0 0.5 // print out .vtp files interval based on simulation time window' '\n')
                file_out.write(' \n')

                file_out.write('W 22 -9.81 // gravity' '\n')
                file_out.write(' \n')

        return filename



class HOS_ocean(InputFileWriter):
    def __init__(self,n1=256, n2=64, xlen=None, ylen=80,T_stop=100,f_out=1,toler=1.0e-7,n=4,Ta=0,
    depth = 100, Tp_real=10,Hs_real=4.5,gamma=3.3,beta=0.78,random_phases=1,
    tecplot=11,i_out_dim=1,i_3d=1,i_a_3d=0,i_2d=0,i_prob=0,i_sw=0):

        self.n1 = n1 # default is 256 at HOS-ocean-1.5/sources/HOS/variables_3D.f90
        self.n2 = n2 # default is 256 at HOS-ocean-1.5/sources/HOS/variables_3D.f90
        self.T_stop = T_stop
        self.f_out = f_out
        self.toler = toler
        self.n = n
        self.Ta = Ta
        self.depth = depth #np.mean(self.topo()[self.land_sea_mask()])
        self.Tp_real = Tp_real
        self.Hs_real = Hs_real
        self.gamma = gamma
        self.beta = beta
        self.random_phases = random_phases
        self.tecplot = tecplot
        self.i_out_dim = i_out_dim
        self.i_3d = i_3d
        self.i_a_3d = i_a_3d
        self.i_2d = i_2d
        self.i_prob = i_prob
        self.i_sw = i_sw
        if xlen is None:
            self.xlen = (self.n1*9.81*self.Tp_real**2)/(8*2*np.pi) # according to n1*2*np.pi/xlen = 5 k_p
            self.ylen = (self.n2*9.81*self.Tp_real**2)/(8*2*np.pi)
        else:
            self.xlen = xlen
            self.ylen = ylen
        return

    def _extension(self):
        return 'dat'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary,
                start_time: str, end_time: str, filename: str,
                grid_path: str, forcing_path: str, boundary_path: str):

        # Create input file name
        __, folder = file_module.split_filepath(filename)
        filename.create_folder(folder=folder+'/Results')

        with open(filename, 'w') as file_out:
            file_out.write(
                'Restart previous computation :: i_restart        :: 0\n')
            file_out.write('Choice of computed case      :: i_case           :: 3\n')

            file_out.write('--- Geometry of the horizontal domain\n')
            file_out.write(
                'Length in x-direction        :: xlen             :: '+format(self.xlen,".1f")+'\n')
            file_out.write(
                'Length in y-direction        :: ylen             :: '+format(self.ylen,".1f")+'\n')

            file_out.write('--- Time stuff \n')
            file_out.write(
                'Duration of the simulation   :: T_stop           :: '+format(self.T_stop,".1f")+'\n')
            file_out.write(
                'Sampling frequency (output)  :: f_out            :: '+format(self.f_out,".1f")+'\n')
            file_out.write(
                'Tolerance of RK scheme       :: toler            :: '+format(self.toler,".2e")+'\n')
            file_out.write(
                'Dommermuth initialisation    :: n                :: '+format(self.n,".0f")+'\n')
            file_out.write(
                'Dommermuth initialisation    :: Ta               :: '+format(self.Ta,".1f")+'\n')

            file_out.write('--- Physical dimensional parameters \n')
            file_out.write(
                'Gravity                      :: grav             :: 9.81\n')
            file_out.write(
                'Water depth                  :: depth            :: '+format(self.depth,".1f")+'\n')

            file_out.write('--- Irregular waves (i_case=3) \n')
            file_out.write(
                'Peak period in s             :: Tp_real          :: '+format(self.Tp_real,".1f")+'\n')
            file_out.write(
                'Significant wave height in m :: Hs_real          :: '+format(self.Hs_real,".1f")+'\n')
            file_out.write(
                'JONSWAP Spectrum             :: gamma            :: '+format(self.gamma,".1f")+'\n')
            file_out.write(
                'Directionality (Dysthe)      :: beta             :: '+format(self.beta,".5f")+'\n')
            file_out.write(
                'Random phases generation     :: random_phases    :: '+format(self.random_phases,".0f")+'\n')

            file_out.write('--- Output files \n')
            file_out.write(
                'Tecplot version              :: tecplot          :: '+format(self.tecplot,".0f")+'\n')
            file_out.write(
                'Output: 1-dim. ; 0-nondim.   :: i_out_dim        :: '+format(self.i_out_dim,".0f")+'\n')
            file_out.write(
                '3d free surface quantities   :: i_3d             :: '+format(self.i_3d,".0f")+'\n')
            file_out.write(
                '3d modes                     :: i_a_3d           :: '+format(self.i_a_3d,".0f")+'\n')
            file_out.write(
                '2d free surface, center line :: i_2d             :: '+format(self.i_2d,".0f")+'\n')
            file_out.write(
                'Wave probes in domain        :: i_prob           :: '+format(self.i_prob,".0f")+'\n')
            file_out.write(
                'Swense output 1="yes",0="no" :: i_sw             :: '+format(self.i_sw,".0f")+'\n')
        return filename

class WW3_grid(InputFileWriter):
    def __init__(self):
        self.scaling = 10**6
        return

    def _extension(self):
        return 'nml'

    def __call__(self, grid: Grid, forcing: Forcing, boundary: Boundary,
                start_time: str, end_time: str, filename: str,
                grid_path: str, forcing_path: str, boundary_path: str):
#         &RECT_NML
#   RECT%NX                =  147
#   RECT%NY                =  126
#   RECT%SX               = 965753.       ! grid increment along x-axis
#   RECT%SY               = 448000.       ! grid increment along y-axis
#   RECT%SF               = 100000000.       ! scaling division factor for x-y axis
#   RECT%X0               = 5.39       ! x-coordinate of lower-left corner (deg)
#   RECT%Y0               = 62.05       ! y-coordinate of lower-left corner (deg)
#   RECT%SF0              = 1.       ! scaling division factor for x0,y0 coord
#
# /
        nx = grid.nx()
        ny = grid.ny()

        sx = round(grid.dlon()*self.scaling)
        sy = round(grid.dlat()*self.scaling)

        sf = self.scaling

        x0 = min(grid.lon())
        y0 = min(grid.lat())
        sf0 = 1.


        with open(filename, 'w') as file_out:
            file_out.write('&RECT_NML\n')
            file_out.write(f'  RECT%NX          = {nx:.0f}\n')
            file_out.write(f'  RECT%NY          = {ny:.0f}\n')
            file_out.write(f'  RECT%SX          = {sx:.0f}.\n')
            file_out.write(f'  RECT%SY          = {sy:.0f}.\n')
            file_out.write(f'  RECT%SF          = {sf:.0f}.\n')
            file_out.write(f'  RECT%X0          = {x0}\n')
            file_out.write(f'  RECT%Y0          = {y0}\n')
            file_out.write(f'  RECT%SF0         = {sf0:.0f}.\n')
            file_out.write('/')

        return filename
