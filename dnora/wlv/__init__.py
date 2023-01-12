"""
This is the water level-module that is responsible for creating and manipulating
the WaterLevel-object.

The WaterLevel-object and the abstract classes used by the other sub-
modules are defined in wlv_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains classes that are used to read wind data from different sources.

write: Contain methods to write the wind data file formats used by e.g. SWAN

Copyright 2022, Øystein Lande, DNV
"""

from .wlv_mod import WaterLevel
from .read import *
#from .read_metno import *
from .read_ec import *
from .write import *
