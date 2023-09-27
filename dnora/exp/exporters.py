from .exp_mod import DataExporter
from . import grd, bnd, wnd, spc, wsr, wlv, inp


class Null(DataExporter):
    def _get_boundary_writer(self):
        return bnd.Null()

    def _get_spectra_writer(self):
        return spc.Null()

    def _get_forcing_writer(self):
        return wnd.Null()

    def _get_grid_writer(self):
        return grd.Null()

    def _get_waterlevel_writer(self):
        return wlv.Null()

    def _get_input_file_writer(self):
        return inp.Null()


class SWAN(DataExporter):
    def _get_default_format(self):
        return "SWAN"

    def _get_boundary_writer(self):
        return bnd.SWAN()

    def _get_spectra_writer(self):
        return spc.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.SWAN()

    def _get_grid_writer(self):
        return grd.SWAN()

    def _get_input_file_writer(self):
        return inp.SWAN()


class SWASH(DataExporter):
    def _get_default_format(self):
        return "SWASH"

    def _get_boundary_writer(self):
        return bnd.SWAN()

    def _get_forcing_writer(self):
        return wnd.SWAN()

    def _get_grid_writer(self):
        return grd.SWAN()

    def _get_input_file_writer(self):
        return inp.SWASH()


class WW3(DataExporter):
    def _get_default_format(self):
        return "WW3"

    def _get_boundary_writer(self):
        return bnd.WW3()

    def _get_spectra_writer(self):
        return spc.DumpToNc()

    def _get_waveseries_writer(self):
        return wsr.DumpToNc()

    def _get_forcing_writer(self):
        return wnd.WW3()

    def _get_grid_writer(self):
        return grd.WW3()

    def _get_trigrid_writer(self):
        return grd.WW3Triangular()


class HOS_ocean(DataExporter):
    def _get_default_format(self):
        return "HOS_ocean"

    def _get_boundary_writer(self):
        return bnd.HOS_ocean()

    def _get_input_file_writer(self):
        return inp.HOS_ocean()


class REEF3D(DataExporter):
    def _get_default_format(self):
        return "REEF3D"

    def _get_spectral_writer(self):
        return spc.REEF3D()

    def _get_input_file_writer(self):
        return inp.REEF3D()
