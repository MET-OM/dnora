from .exp_mod import DataExporter
from .. import grd, bnd, wnd, spc, wsr, wlv, inp

class SWAN(DataExporter):
    def _get_default_format(self):
        return 'SWAN'

    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_spectral_writer(self):
        return spc.write.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_grid_writer(self):
        return grd.write.SWAN()

    def _get_input_file_writer(self):
        return inp.SWAN()

class SWASH(DataExporter):
    def _get_default_format(self):
        return 'SWASH'

    def _get_boundary_writer(self):
        return bnd.write.SWAN()

    def _get_forcing_writer(self):
        return wnd.write.SWAN()

    def _get_grid_writer(self):
        return grd.write.SWAN()

    def _get_input_file_writer(self):
        return inp.SWASH()

class WW3(DataExporter):
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

    def _get_grid_writer(self):
        return grd.write.WW3()

class HOS_ocean(DataExporter):
    def _get_default_format(self):
        return 'HOS_ocean'

    def _get_boundary_writer(self):
        return bnd.write.HOS_ocean()

    def _get_input_file_writer(self):
        return inp.HOS_ocean()

class REEF3D(DataExporter):
    def _get_default_format(self):
        return 'REEF3D'

    def _get_spectral_writer(self):
        return spc.write.REEF3D()

    def _get_input_file_writer(self):
        return inp.REEF3D()
