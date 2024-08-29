from dnora.grid import Grid
from .read import EMODNET, GEBCO

# from dnora.read.generic_readers import ConstantGriddedData


# class Constant(Grid):
#     _default_reader = ConstantGriddedData(topo=999.0)


class EMODNET(Grid):
    _default_reader = EMODNET()


class GEBCO(Grid):
    _default_reader = GEBCO()
