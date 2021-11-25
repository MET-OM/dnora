from copy import copy
from ..grd.grd_mod import Grid
from ..bnd.read import BoundaryReader
from ..bnd.write import BoundaryWriter
from ..wnd.read import ForcingReader
from ..wnd.write import ForcingWriter

from ..bnd.pick import PointPicker

from ..bnd.bnd_mod import Boundary
from ..wnd.wnd_mod import Forcing

class ModelRun:
    def __init__(self, grid: Grid, start_time: str, end_time: str,
                name='AnonymousModel'):
        self._name = copy(name)
        # boundary_reader: BoundaryReader=None, boundary_writer: BoundaryWriter=None,
        # forcing_reader: ForcingReader=None, forcing_writer: ForcingWriter=None,
        # point_picker: PointPicker=None,

        # self._boundary_reader = copy(boundary_reader)
        # self._boundary_writer = copy(boundary_writer)
        # self._forcing_reader = copy(forcing_reader)
        # self._forcing_writer = copy(forcing_writer)
        # self._point_picker = copy(point_picker)
        self.grid = copy(grid)
        self.start_time = start_time
        self.end_time = end_time
        return

    # def set_run_time(self, start_time: str, end_time: str):
    #     self.start_time = start_time
    #     self.end_time = end_time
    #     return
    #
    # def set_grid(self, grid: Grid):
    #     self.grid = copy(grid)
    #     return

    def import_boundary(self, boundary_reader: BoundaryReader=None, point_picker: PointPicker=None, name: str=''):
        if boundary_reader is None:
            self._boundary_reader = self._get_boundary_reader()
        else:
            self._boundary_reader = boundary_reader

        if point_picker is None:
            self._point_picker = self._get_point_picker()
        else:
            self._point_picker = point_picker

        if self._boundary_reader is None :
            raise Exception('Define a BoundaryReader!')
        elif self._point_picker is None:
            raise Exception('Define a PointPicker!')
        else:
            if self.grid is None:
                self.grid = self._get_grid()

            if self.grid is None:
                raise Exception('Define a Grid!')
            else:
                if name:
                    self.boundary = Boundary(grid=self.grid, name=name)
                else:
                    self.boundary = Boundary(grid=self.grid, name=type(self._boundary_reader).__name__)
                self.boundary.import_boundary(start_time=self.start_time, end_time=self.end_time, boundary_reader=self._boundary_reader, point_picker=self._point_picker)

        return None

    def import_forcing(self, forcing_reader: ForcingReader=None, name: str=''):
        if forcing_reader is None:
            self._forcing_reader = self._get_forcing_reader()
        else:
            self._forcing_reader = forcing_reader

        if self._forcing_reader is None:
            raise Exception('Define a ForcingReader!')
        else:
            if self.grid is None:
                self.grid = self._get_grid()

            if self.grid is None:
                raise Exception('Define a Grid!')
            else:
                if name:
                    self.forcing = Forcing(grid=self.grid, name=name)
                else:
                    self.forcing = Forcing(grid=self.grid, name=type(self._forcing_reader).__name__)
                self.forcing.import_forcing(start_time=self.start_time, end_time=self.end_time, forcing_reader=self._forcing_reader)

        return None

    def export_boundary(self, boundary_writer: BoundaryWriter=None):
        if boundary_writer is None:
             self._boundary_writer = self._get_boundary_writer()
        else:
            self._boundary_writer = boundary_writer

        if self._boundary_writer is None:
            raise Exception('Define a BoundaryWriter!')
        else:
            self.boundary.export_boundary(boundary_writer=self._boundary_writer)

        return None

    def export_forcing(self, forcing_writer: ForcingWriter=None):
        if forcing_writer is None:
             self._forcing_writer = self._get_forcing_writer()
        else:
            self._forcing_writer = forcing_writer

        if self._forcing_writer is None:
            raise Exception('Define a ForcingWriter!')
        else:
            self.forcing.export_forcing(forcing_writer=self._forcing_writer)

        return None

    def _get_boundary_reader(self):
        return None

    def _get_boundary_writer(self):
        return None

    def _get_forcing_reader(self):
        return None

    def _get_forcing_writer(self):
        return None

    def _get_point_picker(self):
        return None

    def _get_grid(self):
        return None
