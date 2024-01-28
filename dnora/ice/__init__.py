"""
This is the Ice-module that is responsible for creating and manipulating
the Ice-object.

The Ice-object and the abstract classes used by the other sub-
modules are defined in ice_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains classes that are used to read ice data from different sources.

write: Contain methods to write the ice data file formats used by e.g. SWAN or
WAVEWATCH III to be used as a ice forcing.

Copyright 2023, Konstantinos Christakos and Jan-Victor Björkqvist, MET Norway
"""

from .ice_mod import Ice
from . import read_metno