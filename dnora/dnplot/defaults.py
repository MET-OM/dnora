import cmocean.cm

default = {'hs': {'name': 'Significant wave height', 'unit':'m', 'cmap': cmocean.cm.amp},
           'ff': {'name': 'Wind', 'unit':'m/s', 'cmap': cmocean.cm.tempo},
           'topo':{'name': 'Topography', 'unit':'m', 'cmap': cmocean.cm.topo_r},
           'mask':{'name': ' ', 'unit':' ', 'cmap': 'gray'}
           }
