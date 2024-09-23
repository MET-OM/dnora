from dnora.grid import Grid
from dnora.read.grid.grid_readers import EMODNET, GEBCO, KartverketNo50m
from dnora.read.grid.nchmf import DepFile
from dnora.read.generic.generic_readers import ConstantData


class EMODNET(Grid):
    _default_reader = EMODNET()


class GEBCO(Grid):
    _default_reader = GEBCO()


class NCHMF(Grid):
    _default_reader = DepFile()


class Constant(Grid):
    _default_reader = ConstantData()


class Kartverket(Grid):
    _default_reader = KartverketNo50m()
