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
        gp.wave.Spr: "DSPR",
        gp.ocean.WaterDepth: "DEPTH",
        gp.wind.Wind: "WIND",
        gp.wave.Tm02: "TM02",
        gp.ocean.Current: "VEL",
        gp.ocean.WaterLevel: "WATLEV",
    }
    _output_vars = [
        gp.wave.Hs,
        gp.wave.Tp,
        "RTP",
        gp.wave.Tm01,
        gp.wave.Tm_10,
        gp.wave.Tm02,
        gp.wave.Dirm,
        gp.wave.Dirp,
        gp.wave.Spr,
        gp.ocean.WaterDepth,
        gp.wind.Wind,
    ]

    def _get_default_format(self):
        return ModelFormat.SWAN


class SWASH(ModelExecuter):

    _input_file_writers = {DnoraFileType.INPUT: inputfile.SWASH()}
    _model_runners = {DnoraFileType.INPUT: model_runners.SWASH()}

    def _get_default_format(self):
        return ModelFormat.SWASH


class WW3(ModelExecuter):
    _output_var_aliases = {
        gp.wave.Hs: "HS",
        gp.wave.Tp: "TP",
        gp.wave.Lm_10: "LM",
        gp.wave.Dirm: "DIR",
        gp.wave.Spr: "SPR",
        gp.wave.Dirp: "DP",
        gp.wave.Tm02: "T02",
        gp.wave.Tm_10: "T0M1",
        gp.wave.Tm01: "T01",
        gp.wind.FrictionVelocity: "UST",
        gp.ocean.WaterDepth: "DPT",
        gp.wind.Wind: "WND",
        gp.wave.Stokes: "USS",
        gp.wave.Fp: "FP",
        gp.wave.Km: "WNM",
        gp.wave.HsIG: "HIG",
        gp.wave.Current: "CUR",
        gp.ocean.WaterLevel: "WLV",
        gp.ocean.IceFraction: "ICE",
        gp.ocean.IceThickness: "IC1",
    }
    _output_vars = [
        gp.wave.Hs,
        gp.wave.Tp,
        gp.wave.Tm01,
        gp.wave.Tm02,
        gp.wave.Tm_10,
        gp.wave.Dirm,
        gp.wave.Dirp,
        gp.wave.Spr,
        gp.wind.Wind,
        gp.wind.FrictionVelocity,
        gp.ocean.WaterDepth,
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
