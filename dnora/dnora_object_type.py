from enum import Enum


class DnoraObjectType(Enum):
    ModelRun = "modelrun"
    Forcing = "forcing"
    Boundary = "boundary"
    Grid = "grid"
    TriGrid = "trigrid"
    WaveSeries = "waveseries"
    OceanCurrent = "oceancurrent"
    WaterLevel = "waterlevel"
    Spectra = "spectra"
    Ice = "ice"
    SpectralGrid = "spectral_grid"
    InputFile = "input_file"


def object_type_from_string(obj_str: str) -> DnoraObjectType:
    if isinstance(obj_str, DnoraObjectType):
        return obj_str

    try:
        obj_type = DnoraObjectType[obj_str]
        return obj_type
    except KeyError:
        types = [o for o in DnoraObjectType if o.value == obj_str]
        return types[0]
