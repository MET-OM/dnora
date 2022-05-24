import cmocean.cm

default = {'hs': {'name': 'Significant wave height', 'unit':'m', 'cmap': cmocean.cm.amp},
           'ff': {'name': 'Wind', 'unit':'m/s', 'cmap': cmocean.cm.tempo},
           'topo':{'name': 'Topography', 'unit':'m', 'cmap': cmocean.tools.crop_by_percent(cmocean.cm.topo_r,50, which='min')},
           'mask':{'name': ' ', 'unit':' ', 'cmap': 'gray'}
           }
