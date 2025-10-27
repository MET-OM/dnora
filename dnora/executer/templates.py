from .executer import ModelExecuter
from . import model_runners
from . import inputfile
from dnora.type_manager.model_formats import ModelFormat
from dnora.type_manager.dnora_types import DnoraFileType
import geo_parameters as gp


class SWAN(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWAN()}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWAN()}
    _output_var_aliases = {
        gp.wave.Hs: "HSIGN",
        gp.wave.Tp: "TPS",
        gp.wave.Tm01: "TM01",
        gp.wave.Tm_10: "TMM10",
        gp.wave.Dirm: "DIR",
        gp.wave.Dirp: "PDIR",
        gp.ocean.WaterDepth: "DEP",
    }
    output_vars = [
        gp.wave.Hs,
        gp.wave.Tp,
        "RTP",
        gp.wave.Tm01,
        gp.wave.Tm_10,
        gp.wave.Dirm,
        gp.wave.Dirp,
        gp.ocean.WaterDepth,
    ]

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWASH()}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWASH()}

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(ModelExecuter):
    output_vars = [
        "HS",
        "LM",
        "TP",
        "DIR",
        "SPR",
        "DP",
        "T02",
        "T0M1",
        "T01",
        "UST",
        "CHA",
        "DPT",
        "WND",
        "USS",
        "TUS",
        "TAW",
        "TWO",
        "TOC",
        "FAW",
        "FOC",
        "PHS",
        "PTP",
        "PTM10",
        "PT01",
        "PT02",
        "PDIR",
        "PDP",
        "MXE",
        "MXH",
        "MXHC",
        "SDMH",
        "SDMHC",
        "ABR",
        "UBR",
        "FBB",
        "TBB",
        "CGE",
        "WCC",
        "WBT",
    ]
    _input_file_writers = {
        DnoraFileType.INPUT: inputfile.WW3(),
        DnoraFileType.GRID: inputfile.WW3Grid(),
        DnoraFileType.WIND: inputfile.WW3Forcing(DnoraFileType.WIND),
        DnoraFileType.CURRENT: inputfile.WW3Forcing(DnoraFileType.CURRENT),
        DnoraFileType.WATERLEVEL: inputfile.WW3Forcing(DnoraFileType.WATERLEVEL),
        DnoraFileType.ICE: inputfile.WW3Forcing(DnoraFileType.ICE),
        DnoraFileType.SPECTRA: inputfile.WW3Spectra(),
    }
    _model_runners = {
        DnoraFileType.GRID: model_runners.WW3("grid"),
        DnoraFileType.WIND: model_runners.WW3("prnc"),
        DnoraFileType.SPECTRA: model_runners.WW3("bounc"),
        DnoraFileType.INPUT: model_runners.WW3("shel"),
        DnoraFileType.CURRENT: model_runners.WW3("prnc"),
        DnoraFileType.WATERLEVEL: model_runners.WW3("prnc"),
        DnoraFileType.ICE: model_runners.WW3("prnc"),
    }

    def _get_default_format(self):
        return ModelFormat.WW3
