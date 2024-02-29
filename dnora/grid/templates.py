from dnora.grid import Grid
from .read import EMODNET

# from dnora.readers.generic_readers import ConstantGriddedData


# class Constant(Grid):
#     _default_reader = ConstantGriddedData(topo=999.0)


class EMODNET(Grid):
    _default_reader = EMODNET()
