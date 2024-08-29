"""
This is the Grid-module that is responsible for creating and manipulating the
Grid-object.

The Grid-object and the abstract classes used by the other sub-
modules are defined in grd_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains the classes that are used to read bathymetrical information.

boundary: Contains the classes that define different methods to set boundary
points in the grid.

process: Contain classes that are used to modify the bathymetrical information,
e.g. possible filters.

mesh: Contains classes that are used to mesh the mathymetry to the defined grid,
e.g. using interpolation.

write: Contain methods to write the grid in file formats used by e.g. SWAN or
WAVEWATCH III.

Copyright 2021, Konstantinos Christakos and Jan-Victor Björkqvist, MET Norway
"""

from .grid import Grid
from .trigrid import TriGrid
from . import mask, process, mesh, tri_arangers
from .templates import *
