"""Default filestrings and datestrings for Grid"""
dflt_grd = {'fs': { 'SWAN': '$Grid_SWAN',
                    'WW3': '$Grid',
                    'General': 'topo_$Grid'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'General': ''},
            'ext': {'SWAN': 'bot',
                    'WW3': 'txt',
                    'General': 'txt'}
        }

"""Default filestrings and datestrings for Forcing """
dflt_frc = {'fs': { 'SWAN': 'wind$Forcing$Grid$T0_$T1',
                    'WW3': 'wind_$Forcing_$Grid_$T0-$T1',
                    'General': 'wind_$Forcing_$Grid_$T0-$T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'General': ''},
            'ext': {'SWAN': 'asc',
                    'WW3': 'nc',
                    'General': 'nc'}
        }

"""Default filestrings and datestrings for Boundary"""
dflt_bnd = {'fs': { 'SWAN': 'spec$Boundary$Grid$T0_$T1',
                    'WW3': 'ww3_spec_E$LonN$Lat_$Boundary_$Grid_$T0-$T1',
                    'General': 'spec_$Boundary_$Grid_$T0-$T1'},
            'ds': { 'SWAN': '%Y%m%d',
                    'WW3': '%Y%m%dT%H%M',
                    'General': '%Y%m%dT%H%M'},
            'fldr': {'SWAN': '',
                    'WW3': '',
                    'General': ''},
            'ext': {'SWAN': 'asc',
                    'WW3': 'nc',
                    'General': 'nc'}
        }
