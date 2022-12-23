"""
This is the ocean current-module that is responsible for creating and manipulating
the Forcing-object.

The Forcing-object and the abstract classes used by the other sub-
modules are defined in wnd_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains classes that are used to read ocean current data from different sources.

write: Contain methods to write the wind data file formats used by e.g. SWAN or
WAVEWATCH III to be used as a ocean current forcing.

Copyright 2022, Konstantinos Christakos and Jan-Victor Björkqvist, MET Norway
"""

from .ocr_mod import OceanCurrent
from .read import *
from .read_metno import *
#from .read_ec import *
from .write import *
