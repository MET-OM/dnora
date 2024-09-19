from dnora.grid import Grid
from dnora.read.grid.grid_readers import EMODNET, GEBCO


class EMODNET(Grid):
    _default_reader = EMODNET()


class GEBCO(Grid):
    _default_reader = GEBCO()

class NCHMF(Grid):
    _default_reader = dnora.read.grid.grid_readers.nchmf.DepFile()
