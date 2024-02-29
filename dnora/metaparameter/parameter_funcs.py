from . import parameters
import inspect
from .metaparameter import MetaParameter


def list_of_parameters() -> list:
    parameter_list = []
    for __, obj in inspect.getmembers(parameters):
        if inspect.isclass(obj):
            if issubclass(obj, MetaParameter) and not obj == MetaParameter:
                parameter_list.append(obj)
    return parameter_list


def dict_of_parameters(short: bool = False, alias: bool = False) -> dict:
    if short:
        return {c.short_name(): c for c in list_of_parameters()}
    return {c.standard_name(alias=alias): c for c in list_of_parameters()}


def get(key: str):
    return (
        dict_of_parameters().get(key)
        or dict_of_parameters(alias=True).get(key)
        or dict_of_parameters(short=True).get(key)
    )


def create_metaparameter_dict(parameter_strings: list[str]):
    metaparameter_dict = {}
    for param in parameter_strings:
        val = get(param)
        if val is not None:
            metaparameter_dict[param] = val
    return metaparameter_dict


def set_metaparameters_in_object(obj, metaparameter_dict, data_dict):
    for key, value in data_dict.items():
        metaparameter = metaparameter_dict.get(
            key
        )  # Check if metaparameter provided by reader

        if metaparameter is None:
            # DNORA object usually has specified the metaparameters
            if hasattr(obj, "meta_dict"):
                metaparameter = obj.meta_dict.get(key)

        if metaparameter is not None:
            obj.set_metadata(metaparameter.meta_dict(), data_array_name=key)

    return obj
