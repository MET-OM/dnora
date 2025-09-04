import xarray as xr
import geo_parameters as gp
from dnora import msg


def compile_data_vars(data_vars, aliases: dict):
    """Compiles data_vars accounting for possible aliases etc."""
    new_vars = []
    for var in data_vars:
        name, param = gp.decode(var)
        if param is not None:
            name = param.standard_name()
        if aliases.get(name) is not None:
            if gp.is_gp_class(param):
                param = param(aliases.get(name))  # Create e.g. Dirp('dp')
            else:
                param.name = aliases.get(name)
            msg.plain(
                f"Identified that {type(param).__name__} is called {aliases.get(name)} in the source!"
            )
            new_vars.append(param)
        else:
            new_vars.append(name)
    return new_vars


def read_data_vars(
    data_vars: list[str, gp.metaparameter.MetaParameter],
    ds: xr.Dataset,
    keep_gp_names: bool,
    keep_source_names: bool,
    decode_cf: bool,
) -> dict:
    """Read the given variables from the xarray dataset

    Uninitilized geo-parameter, e.g. gp.wave.Hs
        - Decode using standard-name. If not found, try default name of parameter (e.g. 'hs').
    Uninitilized geo-parameter, e.g. gp.wave.Hs('hsig')
        - Decode using standard-name. If not found, try given name of parameter ('hsig').
    String, e.g. 'hs'
        - Try to read that variable from dataset

    keep_source_names = True
        - Keep the variable name that is used by the data source even if a geo-parameter is given
        - E.g. gp.wave.Hs can match to 'hsig' by standard_name and 'hsig' will be used ('hs' if false)

    keep_gp_names = True
        - Create the variable names using the default class name even if it is initialized
        - E.g. gp.wave.Hs('swh') matches 'swh' variable with no metadata, but the variable created is still called 'hs'

    decode_cf = True
        - Try to decode standard_name from source and create geo-parameters if only string given
        - E.g. given 'hs' will try to find standard name in metadata of 'hs' and create gp.wave.Hs
    """

    if keep_gp_names and keep_source_names:
        raise ValueError("keep_source_names and keep_gp_names cannot both be True!")
    data_dict = {}
    for var in data_vars:
        var, param = gp.decode(var, init=True)
        if param is not None:
            var_in_ds = param.find_me_in_ds(ds)
            if var_in_ds is None:
                var_in_ds = var
            if keep_gp_names:
                param.name = gp.get(param.standard_name()).name
            elif keep_source_names:
                param.name = var_in_ds

            if hasattr(ds[var_in_ds], "standard_name"):
                std_name = ds[var_in_ds].standard_name
            else:
                std_name = ""

            msg.plain(
                f"[{type(param).__name__}('{param.name}')] << '{var_in_ds}' {std_name}"
            )
            data_dict[param] = ds.get(var_in_ds).data
        elif ds.get(var) is not None:
            if hasattr(ds[var], "standard_name"):
                std_name = ds[var].standard_name
            else:
                std_name = ""
            if std_name and decode_cf:
                param = gp.get(std_name)
            else:
                param = None
            if param is not None:
                if keep_gp_names:
                    param = param()
                else:
                    param = param(var)
                msg.plain(
                    f"[{type(param).__name__}('{param.name}')] << '{var}' {std_name}"
                )
                data_dict[param] = ds.get(var).data
            else:
                msg.plain(f"'{var}' << {var}")
                data_dict[var] = ds.get(var).data
        else:
            msg.plain(f"'{var}' unknown. Skipping.")
    return data_dict
