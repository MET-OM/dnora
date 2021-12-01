"""Default filestrings and datestrings for Grid"""
dflt_grd = {'fs': { 'SWAN': '#Grid_SWAN.bot',
                    'WW3': '#Grid.txt',
                    'SWASH': '#Grid_SWASH.bot',
                    'General': 'topo_#Grid.txt'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'SWASH': '',
                    'General': ''},
            'info': {   'SWAN': '#Grid_info.txt',
                        'WW3': '#Grid_info.txt',
                        'SWASH': '#Grid_info.txt',
                        'General': '#Grid_info.txt'}
        }

"""Default filestrings and datestrings for Forcing """
dflt_frc = {'fs': { 'SWAN': 'wind#Forcing#Grid#T0_#T1.asc',
                    'WW3': 'wind_#Forcing_#Grid_#T0-#T1.nc',
                    'SWASH': 'wind#Forcing#Grid#T0_#T1.asc',
                    'General': 'wind_#Forcing_#Grid_#T0-#T1.nc'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'SWASH': '',
                    'General': ''}
        }

"""Default filestrings and datestrings for Boundary"""
dflt_bnd = {'fs': { 'SWAN': 'spec#Boundary#Grid#T0_#T1.asc',
                    'WW3': 'ww3_spec_E#LonN#Lat_#Boundary_#Grid_#T0-#T1.nc',
                    'SWASH': 'spec#Boundary#Grid#T0_#T1.asc',
                    'General': 'spec_#Boundary_#Grid_#T0-#T1.nc'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'SWASH': '%Y%m%d',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'SWASH': '',
                    'General': ''}
        }

"""Default filestrings and datestrings for inp-module"""
dflt_inp = {'fs': { 'SWAN': 'input_#T0_#Grid.swn',
                    'SWASH': 'input_#T0_#T1_#Grid.sws'},
            'ds': { 'SWAN': '%Y%m%d',
                    'SWASH': '%H%M%S'},
            'fldr': {'SWAN': 'MySWANFolder',
                    'SWASH': 'MySWASHFolder'}
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
list_of_placeholders = ['#Grid', '#Forcing', '#Boundary', '#ModelRun', '#T0', '#T1', 'E#Lon', 'N#Lat', '#Lon', '#Lat']
