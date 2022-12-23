import os

"""Default filestrings and datestrings for Grid"""
dflt_grd = {'fs': { 'SWAN': '#Grid_SWAN',
                    'WW3': '#Grid',
                    'SWASH': '#Grid_SWASH',
                    'HOS_ocean': '#Grid_HOS_ocean',
                    'REEF3D': 'geo',
                    'General': 'topo_#Grid'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'HOS_ocean': f"{os.getcwd()}/output",
                    'REEF3D': f"{os.getcwd()}/output",
                    'General': ''},
            'info': {   'SWAN': '#Grid_info.txt',
                        'WW3': '#Grid_info.txt',
                        'SWASH': '#Grid_info.txt',
                        'HOS_ocean': '#Grid_info.txt',
                        'REEF3D': '#Grid_info.txt',
                        'General': '#Grid_info.txt'},
            'ext': { 'SWAN': 'bot',
                            'WW3': 'txt',
                            'SWASH': 'bot',
                            'HOS_ocean': 'dat',
                            'REEF3D': 'dat',
                            'General': 'txt'}
        }

"""Default filestrings and datestrings for Forcing """
dflt_frc = {'fs': { 'SWAN': 'wind#Forcing#Grid#T0_#T1',
                    'WW3': 'wind_#Forcing_#Grid_#T0-#T1',
                    'SWASH': 'wind#Forcing#Grid#T0_#T1',
                    'General': 'wind_#Forcing_#Grid_#T0-#T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'General': ''},
            'ext': { 'SWAN': 'asc',
                    'WW3': 'nc',
                    'SWASH': 'asc',
                    'General': 'nc'}
        }

"""Default filestrings and datestrings for WaterLevel """
dflt_wlv = {'fs': { 'SWAN': 'waterlevel#Grid#T0_#T1',
                    'WW3': 'waterlevel_#Grid_#T0-#T1',
                    'SWASH': 'waterlevel#Grid#T0_#T1',
                    'General': 'waterlevel_#Grid_#T0-#T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'General': ''},
            'ext': { 'SWAN': 'asc',
                    'WW3': 'nc',
                    'SWASH': 'asc',
                    'General': 'nc'}
        }

"""Default filestrings and datestrings for OceanCurrent """
dflt_ocr = {'fs': { 'SWAN': 'oceancurrent#Grid#T0_#T1',
                    'WW3': 'oceancurrent_#Grid_#T0-#T1',
                    'SWASH': 'oceancurrent#Grid#T0_#T1',
                    'General': 'oceancurrent_#Grid_#T0-#T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'General': ''},
            'ext': { 'SWAN': 'asc',
                    'WW3': 'nc',
                    'SWASH': 'asc',
                    'General': 'nc'}
        }


"""Default filestrings and datestrings for Boundary"""
dflt_bnd = {'fs': { 'SWAN': 'spec#Boundary#Grid#T0_#T1',
                    'WW3': 'ww3_spec_E#LonN#Lat_#Boundary_#Grid_#T0-#T1',
                    'SWASH': 'spec#Boundary#Grid#T0_#T1',
                    'General': 'spec_#Boundary_#Grid_#T0-#T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'General': ''},
            'ext': { 'SWAN': 'asc',
                    'WW3': 'nc',
                    'SWASH': 'asc',
                    'General': 'nc'}
        }


"""Default filestrings and datestrings for Spectra"""
dflt_spc = {'fs': { 'SWAN': 'spec1d#Spectra#Grid#T0_#T1',
                    'WW3': 'ww3_spec1d_E#LonN#Lat_#Spectra_#Grid_#T0-#T1',
                    'SWASH': 'spec1d#Spectra#Grid#T0_#T1',
                    'REEF3D': 'spectrum-file',
                    'General': 'spec1d_#Spectra_#Grid_#T0-#T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'REEF3D': '%Y%m%dT%H%M',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': f"{os.getcwd()}/output",
                    'WW3': f"{os.getcwd()}/output",
                    'SWASH': f"{os.getcwd()}/output",
                    'REEF3D': f"{os.getcwd()}/output",
                    'General': ''},
            'ext': { 'SWAN': 'asc',
                    'WW3': 'nc',
                    'SWASH': 'asc',
                    'REEF3D': 'dat',
                    'General': 'nc'}
        }

"""Default filestrings and datestrings for inp-module"""
dflt_inp = {'fs': { 'SWAN': 'input_#T0_#Grid',
                    'SWASH': 'input_#T0_#T1_#Grid',
                    'HOS_ocean': 'input_HOS',
                    'WW3': 'ww3_',
                    },
            'ds': { 'SWAN': '%Y%m%d',
                    'SWASH': '%H%M%S',
                    'HOS_ocean': '',
                    'WW3': '',
                    },
            'fldr': {'SWAN': 'MySWANFolder',
                    'SWASH': 'MySWASHFolder',
                    'HOS_ocean': 'MyHOS_oceanFolder',
                    'WW3': 'MyWW3Folder',
                    },
            'ext': { 'SWAN': 'swn',
                    'WW3': 'nml',
                    'SWASH': 'sws',
                    'HOS_ocean': 'dat',
                    'General': 'txt'}
        }


"""Default filestrings and datestrings for dnplot-module"""
dflt_plt = {'fs': { 'Grid': '#Grid.pdf'},
            'ds': { 'Grid': ''},
            'fldr': {'Grid': ''}
        }

"""Default filestrings and datestrings for mdl-module"""
dflt_mdl = {'ds': { 'General': '%Y%m%dT%H%M'},
            'fldr': {'General': ''}
        }



# Used to clean up filenames
list_of_placeholders = ['#Grid', '#Forcing', '#WaterLevel', '#OceanCurrent', '#Boundary', '#Spectra', '#ModelRun', '#T0', '#T1', 'E#Lon', 'N#Lat', '#Lon', '#Lat']
