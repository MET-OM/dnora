from ..run import model_executers
from ..wlv.read import WaterLevelReader
from .mdl_mod import ModelRun

from .. import bnd, wnd, grd, wlv


class SWAN(ModelRun):
    def _get_point_picker(self):
        return bnd.pick.Area()

    def _get_model_executer(self):
        return model_executers.SWAN()


class SWASH(ModelRun):
    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_model_executer(self):
        return model_executers.SWASH()


class WW3(ModelRun):
    def _get_point_picker(self):
        return bnd.pick.Area()


class OnePoint(ModelRun):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid().set_mask(grd.mask.All())

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()


class HOS_ocean(ModelRun):
    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_model_executer(self):
        return model_executers.HOS_ocean()


class REEF3D(ModelRun):
    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_model_executer(self):
        return model_executers.REEF3D()


class NORA3(ModelRun):
    def _get_boundary_reader(self):
        return bnd.read_metno.NORA3()

    def _get_forcing_reader(self):
        return wnd.read_metno.NORA3()


class ERA5(ModelRun):
    def _get_boundary_reader(self):
        return bnd.read_ec.ERA5()

    def _get_forcing_reader(self):
        return wnd.read_ec.ERA5()

    def _get_waterlevel_reader(self) -> WaterLevelReader:
        return wlv.read_ec.GTSM_ERA5()


class WAM4km(ModelRun):
    def _get_boundary_reader(self):
        return bnd.read_metno.WAM4km()

    def _get_forcing_reader(self):
        return wnd.read_metno.MEPS()


class SWAN_NORA3(SWAN, NORA3):
    pass


class SWAN_ERA5(SWAN, ERA5):
    pass


class SWAN_WAM4km(SWAN, WAM4km):
    pass


class SWASH_NORA3(SWASH, NORA3):
    pass


class WW3_NORA3(WW3, NORA3):
    pass


class WW3_WAM4km(WW3, WAM4km):
    pass


class OnePoint_NORA3(OnePoint, NORA3):
    pass


class OnePoint_ERA5(OnePoint, ERA5):
    pass
