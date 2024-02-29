from .metaparameter import MetaParameter
import pint

ureg = pint.UnitRegistry()


class WaterDepth(MetaParameter):
    _short_name = "depth"
    _long_name = "water_depth"
    _standard_name = "sea_floor_depth_below_sea_surface"
    _unit = ureg.m


class Hs(MetaParameter):
    _short_name = "hs"
    _long_name = "significant_wave_height"
    _standard_name = [
        "sea_surface_wave_significant_height",
        "significant_height_of_wind_and_swell_waves",
    ]
    _unit = ureg.m


class HsSwell(MetaParameter):
    _short_name = "hs_swell"
    _long_name = "significant_swell_height"
    _standard_name = [
        "sea_surface_swell_wave_significant_height",
        "significant_height_of_swell_waves",
    ]
    _unit = ureg.m


class HsSwell1(MetaParameter):
    _short_name = "hs_swell1"
    _long_name = "significant_primary_swell_height"
    _standard_name = "sea_surface_primary_swell_wave_significant_height"
    _unit = ureg.m


class Tm01(MetaParameter):
    _short_name = "tm01"
    _long_name = "first_moment_mean_wave_period"
    _standard_name = "sea_surface_wave_mean_period_from_variance_spectral_density_first_frequency_moment"
    _unit = ureg.s


class Tm_10(MetaParameter):
    _short_name = "tm_10"
    _long_name = "inverse_moment_mean_wave_period"
    _standard_name = "sea_surface_wave_mean_period_from_variance_spectral_density_inverse_frequency_moment"
    _unit = ureg.s


class DirmSwell(MetaParameter):
    _short_name = "dirm_swell"
    _long_name = "mean_swell_direction"
    _standard_name = "sea_surface_swell_wave_from_direction"
    _unit = ureg.deg


class Dirm(MetaParameter):
    _short_name = "dirm"
    _long_name = "mean_wave_direction"
    _standard_name = "sea_surface_wave_from_direction"
    _unit = ureg.deg


class Dirp(MetaParameter):
    _short_name = "dirp"
    _long_name = "peak_wave_direction"
    _standard_name = (
        "sea_surface_wave_from_direction_at_variance_spectral_density_maximum"
    )
    _unit = ureg.deg


class Spr(MetaParameter):
    _short_name = "spr"
    _long_name = "wave_directional_spread"
    _standard_name = "sea_surface_wave_directional_spread"
    _unit = ureg.deg


class SprP(MetaParameter):
    _short_name = "sprp"
    _long_name = "peak_wave_directional_spread"
    _standard_name = (
        "sea_surface_wave_directional_spread_at_variance_spectral_density_maximum"
    )
    _unit = ureg.deg


class M0(MetaParameter):
    _short_name = "m0"
    _long_name = "zeroth_wave_moment"
    _standard_name = "sea_surface_zeroth_wave_moment"
    _unit = ureg.m**2
    _cf = False


class XWind(MetaParameter):
    _short_name = "x_wind"
    _long_name = "eastward_wind_component"
    _standard_name = [
        "x_wind",
        "grid_eastward_wind",
    ]
    _unit = ureg.m / ureg.s


class YWind(MetaParameter):
    _short_name = "y_wind"
    _long_name = "northward_wind_component"
    _standard_name = [
        "y_wind",
        "grid_northward_wind",
    ]
    _unit = ureg.m / ureg.s


class SeaLevel(MetaParameter):
    _short_name = "eta"
    _long_name = "sea_surface_height"
    _standard_name = [
        "sea_surface_elevation",
        "sea_surface_elevation_anomaly",
        "sea_surface_height_above_geoid",
    ]
    _unit = ureg.m


class XCurrent(MetaParameter):
    _short_name = "x_current"
    _long_name = "eastward_current_component"
    _standard_name = [
        "sea_water_x_velocity",
        "x_sea_water_velocity",
    ]
    _unit = ureg.m / ureg.s


class YCurrent(MetaParameter):
    _short_name = "y_current"
    _long_name = "northward_current_component"
    _standard_name = [
        "sea_water_y_velocity",
        "y_sea_water_velocity",
    ]
    _unit = ureg.m / ureg.s


class IceFraction(MetaParameter):
    _short_name = "ice_fraction"
    _long_name = "sea_ice_fraction"
    _standard_name = "sea_ice_area_fraction"
    _unit = ureg.percent


class IceThickness(MetaParameter):
    _short_name = "ice_thickness"
    _long_name = "sea_ice_thickness"
    _standard_name = "sea_ice_thickness"
    _unit = ureg.m


class Ef(MetaParameter):
    _short_name = "ef"
    _long_name = "spectral_density"
    _standard_name = "sea_surface_wave_variance_spectral_density"
    _unit = ureg.m * ureg.m * ureg.s


class Efth(MetaParameter):
    _short_name = "efth"
    _long_name = "directional_spectral_density"
    _standard_name = "sea_surface_wave_directional_variance_spectral_density"
    _unit = ureg.m * ureg.m * ureg.s / ureg.rad
