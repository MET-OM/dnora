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
    IceForcing = "iceforcing"
    SpectralGrid = "spectral_grid"
    InputFile = "input_file"
