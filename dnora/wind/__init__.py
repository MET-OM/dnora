"""
This is the Forcing-module that is responsible for creating and manipulating
the Forcing-object.

The Forcing-object and the abstract classes used by the other sub-
modules are defined in wnd_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains classes that are used to read wind data from different sources.

write: Contain methods to write the wind data file formats used by e.g. SWAN or
WAVEWATCH III to be used as a wind forcing.

Copyright 2021, Konstantinos Christakos and Jan-Victor Björkqvist, MET Norway
"""

from .wind import Wind
