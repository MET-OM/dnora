from .mdl_mod import ModelRun

from .. import bnd, wnd

class SWAN(ModelRun):
    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()


class WW3(ModelRun):
    def _get_boundary_writer(self):
        return bnd.write.WW3()

    def _get_forcing_writer(self):
        return wnd.write.WW3()

    def _get_point_picker(self):
        return bnd.pick.Area()


class SWAN_NORA3(SWAN):
    def _get_boundary_reader(self):
        return bnd.read.MetNo_NORA3()

    def _get_forcing_reader(self):
        return wnd.read.MetNo_NORA3()


class WW3_NORA3(WW3):
    def _get_boundary_reader(self):
        return bnd.read.MetNo_NORA3()

    def _get_forcing_reader(self):
        return wnd.read.MetNo_NORA3()
