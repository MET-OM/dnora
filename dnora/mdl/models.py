from .mdl_mod import ModelRun

from .. import bnd, wnd, grd, inp, run, spc, wlv, ocr

class SWAN(ModelRun):
    def _get_default_format(self):
        return 'SWAN'

    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_spectral_writer(self):
        return spc.write.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_waterlevel_writer(self):
        return wlv.write.SWAN()

    def _get_oceancurrent_writer(self):
        return ocr.write.SWAN()

    def _get_point_picker(self):
        return bnd.pick.Area()

    def _get_grid_writer(self):
        return grd.write.SWAN()

    def _get_input_file_writer(self):
        return inp.SWAN()

    def _get_model_executer(self):
        return run.SWAN()

class SWASH(ModelRun):
    def _get_default_format(self):
        return 'SWASH'

    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_grid_writer(self):
        return grd.write.SWAN()

    def _get_input_file_writer(self):
        return inp.SWASH()

    def _get_model_executer(self):
        return run.SWASH()


class WW3(ModelRun):
    def _get_default_format(self):
        return 'WW3'

    def _get_boundary_writer(self):
        return bnd.write.WW3()

    def _get_forcing_writer(self):
        return wnd.write.WW3()

    def _get_point_picker(self):
        return bnd.pick.Area()

    def _get_grid_writer(self):
        return grd.write.WW3()

class HOS_ocean(ModelRun):
    def _get_default_format(self):
        return 'HOS_ocean'

    def _get_boundary_writer(self):
        return bnd.write.HOS_ocean()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_input_file_writer(self):
        return inp.HOS_ocean()

    def _get_model_executer(self):
        return run.HOS_ocean()

class REEF3D(ModelRun):
    def _get_default_format(self):
        return 'REEF3D'

    def _get_spectral_writer(self):
        return spc.write.REEF3D()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_input_file_writer(self):
        return inp.REEF3D()

    def _get_model_executer(self):
        return run.REEF3D()

class SWAN_NORA3(SWAN):
    def _get_boundary_reader(self):
        return bnd.read_metno.NORA3()

    def _get_forcing_reader(self):
        return wnd.read_metno.NORA3()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint(max_dist=7, remove_duplicate=True)

class SWAN_ERA5(SWAN):
    def _get_boundary_reader(self):
        return bnd.read_ec.ERA5()

    def _get_forcing_reader(self):
        return wnd.read_ec.ERA5()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint(max_dist=15, remove_duplicate=True)

    def _get_waterlevel_reader(self):
        return wlv.read_ec.GTSM_ERA5()

class SWAN_WAM4km(SWAN):
    def _get_boundary_reader(self):
        return bnd.read_metno.WAM4km()

    def _get_forcing_reader(self):
        return wnd.read_metno.MEPS()

class SWAN_WW3_4km(SWAN):
    def _get_boundary_reader(self):
        return bnd.read_metno.WW3_4km()

    def _get_forcing_reader(self):
        return wnd.read_metno.MEPS()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint(max_dist=7, remove_duplicate=True)

class SWASH_NORA3(SWASH):
    def _get_boundary_reader(self):
        return bnd.read_metno.NORA3()

    def _get_forcing_reader(self):
        return wnd.read_metno.NORA3()

class WW3_NORA3(WW3):
    def _get_boundary_reader(self):
        return bnd.read_metno.NORA3()

    def _get_forcing_reader(self):
        return wnd.read_metno.NORA3()
