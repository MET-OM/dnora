from dnora.grid import Grid
from .read import ConstantTopo, EMODNET


class Constant(Grid):
    _default_reader = ConstantTopo()


class EMODNET(Grid):
    _default_reader = EMODNET()
