from .mdl_mod import ModelRun

from .. import bnd, wnd, grd, inp, run, spc, wsr

class SWAN(ModelRun):
    def _get_default_format(self):
        return 'SWAN'

    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_spectral_writer(self):
        return spc.write.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

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

    def _get_spectral_writer(self):
        return spc.write.DumpToNc()

    def _get_waveseries_writer(self):
        return wsr.write.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.write.WW3()

    def _get_point_picker(self):
        return bnd.pick.Area()

    def _get_grid_writer(self):
        return grd.write.WW3()

class OnePoint(ModelRun):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid().set_mask(grd.mask.All())

    def _get_default_format(self):
        return 'WW3'

    def _get_boundary_writer(self):
        return bnd.write.DumpToNc()

    def _get_spectral_writer(self):
        return spc.write.DumpToNc()

    def _get_waveseries_writer(self):
        return wsr.write.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.write.DumpToNc()

    def _get_point_picker(self):
        return bnd.pick.NearestGridPoint()

    def _get_grid_writer(self):
        return grd.write.DumpToNc()


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
