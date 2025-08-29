"""
This is the Boundary-module that is responsible for creating and manipulating
the Boundary-object.

The Boundary-object and the abstract classes used by the other sub-
modules are defined in bnd_mod.py, and this should not be modified.

In addition it contains the following sub-modules, which can be expanded by
the user:

read: Contains classes that are used to read boundary spectra from
different sources.

pick: Contains classes that define which of the spectra are chosen to be used
as boundary forcing, e.g. within a certain area around the grid.

process: Contain classes that are used to modify the spectr,
e.g. change of directional convention or scaling with a constant.

write: Contain methods to write the spectra in file formats used by e.g. SWAN or
WAVEWATCH III to be used as a boundary forcing.

Copyright 2021, Konstantinos Christakos and Jan-Victor Björkqvist, MET Norway
"""

from .spectra import Spectra
