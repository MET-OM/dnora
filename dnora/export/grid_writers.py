from __future__ import annotations  # For TYPE_CHECKING

import numpy as np
from copy import copy
from .data_writers import DataWriter

# Import objects
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from dnora.modelrun.modelrun import ModelRun
    from dnora.file_module import FileNames

from dnora import msg
from dnora import file_module
from dnora.type_manager.dnora_types import DnoraDataType


class SWAN(DataWriter):
    """Writes the grid to SWAN format."""

    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type, **kwargs
    ) -> str:
        filename = file_object.get_filepath()
        grid = model.grid()

        if grid.boundary_mask().size > 0:
            msg.info(
                f"Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.sea_mask()))):d} boundary points in grid..."
            )
        np.savetxt(filename, grid.topo(), delimiter="\t", fmt="%1.2f")

        return filename


class WW3(DataWriter):
    """Writes the grid to WAVEWATCH III format."""

    def __init__(self, matrix=False) -> None:
        self.matrix = matrix
        return

    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: DnoraDataType, **kwargs
    ) -> list[str]:
        filename = file_object.get_filepath()
        grid = model.grid()

        mask_out = np.zeros(grid.topo().shape)
        mask_out[grid.sea_mask()] = 1
        if grid.boundary_mask().size > 0:
            msg.info(
                f"Setting {sum(sum(np.logical_and(grid.boundary_mask(), grid.sea_mask()))):d} boundary points in grid..."
            )
            mask_out[np.logical_and(grid.boundary_mask(), grid.sea_mask())] = 2

        output_files = []
        if self.matrix:
            output_file = file_module.add_prefix(filename, "mat")
            output_file = file_module.add_suffix(output_file, "bathy")
            output_files.append(output_file)

            # msg.to_file(output_path)
            np.savetxt(output_file, grid.topo(), delimiter=",", fmt="%1.6f")

            output_file = file_module.add_prefix(filename, "mat")
            output_file = file_module.add_suffix(output_file, "mapsta")
            output_files.append(output_file)

            # msg.to_file(output_path)
            np.savetxt(output_file, mask_out, delimiter=",", fmt="%1.0f")

        else:
            output_file = file_module.add_suffix(filename, "bathy")
            output_files.append(output_file)

            # msg.to_file(output_path)
            np.savetxt(output_file, grid.topo().ravel(), delimiter=",", fmt="%1.6f")

            output_file = file_module.add_suffix(filename, "mapsta")
            output_files.append(output_file)

            # msg.to_file(output_path)
            np.savetxt(output_file, mask_out.ravel(), delimiter=",", fmt="%1.0f")

        return output_files


class WW3Triangular(DataWriter):
    """Writes the grid to WAVEWATCH III unstructured format."""

    def __call__(
        self, model: ModelRun, file_object: FileNames, obj_type: DnoraDataType, **kwargs
    ) -> str:
        filename = file_object.get_filepath(extension="msh")
        grid = model.grid()
        output_file = file_module.add_suffix(filename, "bathy")

        with open(output_file, "w") as f:
            # Write header
            f.write("$MeshFormat\n")
            f.write("2 0 8\n")
            f.write("$EndMeshFormat\n")

            # Write Nodes
            arr_node = np.stack(
                [
                    grid.inds() + 1,
                    grid.lon(),
                    grid.lat(),
                    grid.topo(),
                ]
            ).transpose()
            fmt_node = "%10.0f%22.8f%22.8f%22.8f"

            f.write("$Nodes\n")
            f.write(f"{len(grid.lon()):12.0f}\n")
            np.savetxt(f, arr_node, fmt=fmt_node)
            f.write("$EndNodes\n")

            # Write Elements
            N_bound = len(grid.boundary_points()[0])
            arr_bound = np.stack(
                [
                    np.arange(1, N_bound + 1),
                    N_bound * [15],
                    N_bound * [2],
                    N_bound * [1],
                    N_bound * [0],
                    grid.inds()[grid.boundary_mask()] + 1,
                ]
            ).transpose()
            fmt_bound = "%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f"

            N_ele = len(grid.tri())
            arr_ele = np.stack(
                [
                    np.arange(1, N_ele + 1) + N_bound,
                    N_ele * [2],
                    N_ele * [3],
                    N_ele * [0],
                    np.arange(1, N_ele + 1),
                    N_ele * [0],
                    grid.tri()[:, 0] + 1,
                    grid.tri()[:, 1] + 1,
                    grid.tri()[:, 2] + 1,
                ]
            ).transpose()
            fmt_ele = "%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f"

            f.write("$Elements\n")
            f.write(f"{len(grid.tri())+len(grid.boundary_points()[0]):12.0f}\n")
            np.savetxt(f, arr_bound, fmt=fmt_bound)
            np.savetxt(f, arr_ele, fmt=fmt_ele)
            f.write("$EndElements\n")

        return output_file


class Xyz(DataWriter):
    """Writes the grid to Xyz-format."""

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        use_raw: bool = False,
        cartesian: bool = False,
        **kwargs,
    ) -> str:
        filename = file_object.get_filepath(extension="xyz")
        grid = model.grid()

        if self._use_raw:
            grid = model.grid().raw()
        else:
            grid = model.grid()

        if cartesian:
            x = grid.x()
            y = grid.y()
            fmt = ".4f"
        else:
            x = grid.lon()
            y = grid.lat()
            fmt = ".9f"

        z = grid.topo()

        with open(filename, "w") as f:
            f.write(f"{filename}\n")
            for nx, lon in enumerate(x):
                for ny, lat in enumerate(y):
                    lon_out = lon
                    lat_out = lat

                    if ~np.isnan(z[ny, nx]):
                        f.write(f"{lon_out:{fmt}},{lat_out:{fmt}},{z[ny,nx]:.1f}\n")

        return filename


class REEF3D(DataWriter):
    """Writes the grid to Xyz-format in relative Cartesian grid."""

    def __call__(
        self,
        model: ModelRun,
        file_object: FileNames,
        obj_type: DnoraDataType,
        use_raw: bool = False,
        **kwargs,
    ) -> str:
        filename = file_object.get_filepath()

        if self._use_raw:
            grid = model.grid().raw()
        else:
            grid = model.grid()

        z = grid.topo().ravel()
        x, y = grid.xy(normalize=True)

        print("Max depth(m):", np.nanmax(z))
        z = -1 * z + np.nanmax(z)  # set at zero the maximum depth.
        z[np.isnan(z)] = np.nanmax(z) + 14  # + 14 for land points

        fmt = ".5f"
        with open(filename, "w") as f:
            for i, __ in enumerate(x):
                f.write(f"{x[i]:{fmt}} {y[i]:{fmt}} {z[i]:.1f}\n")

        return filename
