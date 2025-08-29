import numpy as np
from dnora.utils.grid import create_swan_segment_coords, identify_boundary_edges

from dnora import msg


def swan_header(file_out, grid_name: str) -> None:
    """Writes header information to SWAN input file"""
    file_out.write("$************************HEADING************************\n")
    file_out.write("$ \n")
    file_out.write(" PROJ '" + grid_name + "' 'T24' \n")
    file_out.write("$ \n")
    file_out.write("$*******************MODEL INPUT*************************\n")
    file_out.write("$ \n")
    file_out.write("SET NAUT \n")
    file_out.write("$ \n")


def swan_grid(
    file_out, grid, grid_path: str, n_dir: int, f_low: float, f_high: float, n_freq: int
) -> None:
    """Writes grid specifications to SWAN input file"""
    delta_X = np.round(np.diff(grid.edges("lon")), 5)[0]
    delta_Y = np.round(np.diff(grid.edges("lat")), 5)[0]
    file_out.write("COORD SPHE CCM \n")
    file_out.write(
        "CGRID "
        + str(grid.lon()[0])
        + " "
        + str(grid.lat()[0])
        + " 0. "
        + str(delta_X)
        + " "
        + str(delta_Y)
        + " "
        + str(grid.nx() - 1)
        + " "
        + str(grid.ny() - 1)
        + " CIRCLE %d %f %f %d \n" % (n_dir, f_low, f_high, n_freq)
    )
    file_out.write("$ \n")

    file_out.write(
        "INPGRID BOTTOM "
        + str(grid.lon()[0])
        + " "
        + str(grid.lat()[0])
        + " 0. "
        + str(grid.nx() - 1)
        + " "
        + str(grid.ny() - 1)
        + " "
        + str((delta_X / (grid.nx() - 1)).round(8))
        + " "
        + str((delta_Y / (grid.ny() - 1)).round(8))
        + "\n"
    )
    file_out.write("READINP BOTTOM 1 '" + grid_path.split("/")[-1] + "' 3 0 FREE \n")
    file_out.write("$ \n")


def swan_spectra(file_out, grid, spectra, boundary_path: str) -> None:
    """Writes information about boundary spectra to SWAN input file"""

    if spectra is None:
        msg.info(
            "No spectra object provided. Spectra information will NOT be written to SWAN input file!"
        )
        return

    lons, lats = create_swan_segment_coords(
        grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
    )

    bound_string = "BOUNDSPEC SEGMENT XY"

    for lon, lat in zip(lons, lats):
        bound_string += f" {lon:.4f} {lat:.4f}"
    bound_string += " VARIABLE FILE 0 "
    bound_string += f"'{boundary_path.split('/')[-1]}'\n"
    file_out.write(bound_string)

    file_out.write("$ \n")


def swan_homog_spectra(file_out, grid, homog: dict) -> None:
    """Writes stationary boundary conditions to SWAN input file"""
    gamma = homog.get("gamma", 3.3)
    file_out.write(f"BOUND SHAP JON {gamma} PEAK DSPR POWER\n")
    for side in identify_boundary_edges(grid.boundary_mask()):
        boundary = homog.get(side)
        if boundary is None:
            boundary = (
                homog  ## If only one set of parameters have been giiven for all sides
            )

        hs = boundary.get("hs")
        if hs is None:
            raise ValueError(
                f"No stationary boundary condition 'hs' given for boundary '{side}'"
            )
        tp = boundary.get("tp")
        if tp is None:
            raise ValueError(
                f"No stationary boundary condition 'tp' given for boundary '{side}'"
            )
        dirp = boundary.get("dirp")
        if dirp is None:
            raise ValueError(
                f"No stationary boundary condition 'dirp' given for boundary '{side}'"
            )
        file_out.write(f"BOUND SIDE {side} CONST PAR {hs:.1f} {tp:.1f} {dirp:.0f}\n")


def swan_wind(
    file_out,
    forcing,
    STR_START,
    STR_END,
    factor: float,
    forcing_path: str,
) -> None:
    """Writes wind information to SWAN input file"""

    if forcing is None:
        msg.info(
            "No wind object provided. Wind information will NOT be written to SWAN input file!"
        )
        return

    delta_Xf = np.round(np.diff(forcing.edges("lon")), 5)[0]
    delta_Yf = np.round(np.diff(forcing.edges("lat")), 5)[0]

    file_out.write(
        "INPGRID WIND "
        + str(forcing.lon()[0].round(3))
        + " "
        + str(forcing.lat()[0].round(3))
        + " 0. "
        + str(forcing.nx() - 1)
        + " "
        + str(forcing.ny() - 1)
        + " "
        + str((delta_Xf / (forcing.nx() - 1)).round(6))
        + " "
        + str((delta_Yf / (forcing.ny() - 1)).round(6))
        + " NONSTATIONARY "
        + STR_START
        + f" {forcing.dt():.0f} HR "
        + STR_END
        + "\n"
    )

    file_out.write(
        "READINP WIND "
        + str(factor)
        + "  '"
        + forcing_path.split("/")[-1]
        + "' 3 0 0 1 FREE \n"
    )
    file_out.write("$ \n")


def swan_waterlevel(
    file_out, waterlevel, STR_START, STR_END, factor: float, waterlevel_path: str
) -> None:
    """Writes waterlevel information to SWAN input file"""
    if waterlevel is None:
        msg.info(
            "No waterlevel object provided. Waterlevel information will NOT be written to SWAN input file!"
        )
        return

    delta_Xf = np.round(np.diff(waterlevel.edges("lon")), 5)[0]
    delta_Yf = np.round(np.diff(waterlevel.edges("lat")), 5)[0]

    file_out.write(
        "INPGRID WLEV "
        + str(waterlevel.lon()[0])
        + " "
        + str(waterlevel.lat()[0])
        + " 0. "
        + str(waterlevel.nx() - 1)
        + " "
        + str(waterlevel.ny() - 1)
        + " "
        + str((delta_Xf / (waterlevel.nx() - 1)).round(6))
        + " "
        + str((delta_Yf / (waterlevel.ny() - 1)).round(6))
        + " NONSTATIONARY "
        + STR_START
        + f" {waterlevel.dt():.0f} HR "
        + STR_END
        + "\n"
    )

    file_out.write(
        "READINP WLEV "
        + str(factor)
        + "  '"
        + waterlevel_path.split("/")[-1]
        + "' 3 0 1 FREE \n"
    )
    file_out.write("$ \n")


def swan_current(
    file_out, current, STR_START, STR_END, factor: float, current_path: str
) -> None:
    """Writes current information to SWAN input file"""

    if current is None:
        msg.info(
            "No current object provided. OceanCurrent information will NOT be written to SWAN input file!"
        )
        return

    delta_Xf = np.round(np.diff(current.edges("lon")), 5)[0]
    delta_Yf = np.round(np.diff(current.edges("lat")), 5)[0]

    file_out.write(
        "INPGRID CUR "
        + str(current.lon()[0].round(3))
        + " "
        + str(current.lat()[0].round(3))
        + " 0. "
        + str(current.nx() - 1)
        + " "
        + str(current.ny() - 1)
        + " "
        + str((delta_Xf / (current.nx() - 1)).round(6))
        + " "
        + str((delta_Yf / (current.ny() - 1)).round(6))
        + " EXC 32767 NONSTATIONARY "
        + STR_START
        + f" {current.dt():.0f} HR "
        + STR_END
        + " \n"
    )

    file_out.write(
        "READINP CUR "
        + str(factor)
        + "  '"
        + current_path.split("/")[-1]
        + "' 3 0 0 1 FREE \n"
    )
    file_out.write("$ \n")


def swan_ice(file_out, ice, STR_START, STR_END, factor: float, ice_path: str) -> None:
    """Writes ice information to SWAN input file"""

    if ice is None:
        msg.info(
            "No ice object provided. Ice information will NOT be written to SWAN input file!"
        )
        return

    delta_Xf = np.round(np.diff(ice.edges("lon")), 5)[0]
    delta_Yf = np.round(np.diff(ice.edges("lat")), 5)[0]
    for i in range(len(ice_path)):
        if ice_path[i].split("/")[-1].startswith("sic") or ice_path[i].split("/")[
            -1
        ].startswith("ice"):
            ICE_NAME = "AICE"
        elif ice_path[i].split("/")[-1].startswith("sit"):
            ICE_NAME = "HICE"
        # self.output_var = self.output_var + " " + ICE_NAME # comment due to AICE/HICE not available as output in nc-format in SWAN
        file_out.write(
            "INPGRID "
            + ICE_NAME
            + " "
            + str(ice.lon()[0].round(3))
            + " "
            + str(ice.lat()[0].round(3))
            + " 0. "
            + str(ice.nx() - 1)
            + " "
            + str(ice.ny() - 1)
            + " "
            + str((delta_Xf / (ice.nx() - 1)).round(6))
            + " "
            + str((delta_Yf / (ice.ny() - 1)).round(6))
            + " NONSTATIONARY "
            + STR_START
            + f" {ice.dt():.0f} HR "
            + STR_END
            + " \n"
        )
        file_out.write(
            "READINP "
            + ICE_NAME
            + " "
            + str(factor)
            + "  '"
            + ice_path[i].split("/")[-1]
            + "' \n"
        )
        file_out.write("$ \n")


def swan_structures(file_out, structures: list[dict]) -> None:
    trans = None
    refl = 0.0
    closed = False
    if structures:
        file_out.write("$ Define obstacles (structures) \n")
    for structure in structures:
        trans = structure.get("trans") or trans
        if trans is None:
            raise ValueError(
                "Provide a value for the transparency of the structure with the dict key 'trans' in every structure. To use a constant transparency, only provide it in the first structure."
            )
        refl = structure.get("refl") or refl

        if structure.get("name") is not None:
            file_out.write(f"$ --- {structure.get('name')}")
            if structure.get("closed"):
                file_out.write(f" [closed]")
            file_out.write(f" ---\n")

        file_out.write(f"OBSTACLE TRANS {trans:.2f} REFL {refl:.2f} LINE")
        for lon, lat in zip(structure.get("lon"), structure.get("lat")):
            file_out.write(f" {lon:.6f} {lat:.6f}")

        if structure.get("closed"):
            file_out.write(
                f" {structure.get('lon')[0]:.6f} {structure.get('lat')[0]:.6f}"
            )
        file_out.write("\n")


def swan_spectral_output_points(file_out, grid, STR_START: str, homog: bool) -> None:
    """Writes a list of spectral output points set in grid"""
    spec_lon, spec_lat = grid.output_points()
    if len(spec_lon) == 0:
        return
    file_out.write("$ Generate spectral output \n")
    file_out.write("POINTS 'pkt' &\n")
    for slon, slat in zip(spec_lon, spec_lat):
        file_out.write(str(slon) + " " + str(slat) + " &\n")
    file_out.write(
        "SPECOUT 'pkt' SPEC2D ABS '"
        + grid.name
        + "_"
        + STR_START.split(".")[0]
        + "_spec"
        + ".nc'"
    )
    if not homog:
        file_out.write(" & \n OUTPUT " + STR_START + " 1 HR")

    file_out.write(" \n")
