from .mdl_mod import ModelRun

from .. import bnd, wnd, grd, inp

class SWAN(ModelRun):
    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_grid_writer(self):
        return grd.write.SWAN()

    def _get_input_file_writer(self):
        return inp.SWAN()

class SWASH(ModelRun):
    def _get_boundary_writer(self):
        return bnd.write.SWAN(out_format='SWASH')

    def _get_forcing_writer(self):
        return wnd.write.SWAN(out_format='SWASH')

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_grid_writer(self):
        return grd.write.SWAN(out_format='SWASH')

    def _get_input_file_writer(self):
        return inp.SWASH()

class WW3(ModelRun):
    def _get_boundary_writer(self):
        return bnd.write.WW3()

    def _get_forcing_writer(self):
        return wnd.write.WW3()

    def _get_point_picker(self):
        return bnd.pick.Area()

    def _get_grid_writer(self):
        return grd.write.WW3()

class SWAN_NORA3(SWAN):
    def _get_boundary_reader(self):
        return bnd.read.MetNo_NORA3()

    def _get_forcing_reader(self):
        return wnd.read.MetNo_NORA3()

class SWASH_NORA3(SWASH):
    def _get_boundary_reader(self):
        return bnd.read.MetNo_NORA3()

    def _get_forcing_reader(self):
        return wnd.read.MetNo_NORA3()

class WW3_NORA3(WW3):
    def _get_boundary_reader(self):
        return bnd.read.MetNo_NORA3()

    def _get_forcing_reader(self):
        return wnd.read.MetNo_NORA3()
