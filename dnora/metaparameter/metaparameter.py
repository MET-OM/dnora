from abc import ABC, abstractstaticmethod
import numpy as np

from . import parameters
import inspect
from pint import Unit


def list_of_parameters() -> list:
    parameter_list = []
    for name, obj in inspect.getmembers(parameters):
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


class MetaParameter(ABC):
    _cf = True

    def __init__(self, name: str = ""):
        self._name = name or self.short_name()

    @property
    @abstractstaticmethod
    def _short_name() -> str:
        pass

    @property
    @abstractstaticmethod
    def _long_name() -> str:
        pass

    @property
    @abstractstaticmethod
    def _standard_name() -> str:
        pass

    @property
    @abstractstaticmethod
    def _unit() -> str:
        pass

    @classmethod
    def cf(cls) -> bool:
        return cls._cf

    @classmethod
    def short_name(cls) -> str:
        return cls._short_name

    @classmethod
    def long_name(cls) -> str:
        return cls._long_name

    @classmethod
    def standard_name(cls, strict: bool = False, alias: bool = False) -> str:
        names = np.atleast_1d(cls._standard_name)
        if cls.cf() or not strict:
            if alias:
                return names[-1]
            return names[0]

        return ""

    @classmethod
    def standard_aliases(cls, strict=False) -> list[str]:
        if cls.cf() or not strict:
            return cls._standard_name

        return [""]

    @classmethod
    def unit(cls) -> Unit:
        return cls._unit

    def quantify(self):
        return self * self.unit()

    def set_name(self, name: str = None) -> None:
        if name is None:
            self._name = self.short_name()
        self._name = name

    def name(self) -> str:
        return self._name or self.short_name()

    @classmethod
    def cf(cls) -> str:
        return cls._cf

    @classmethod
    def meta_dict(cls, alias: bool = False) -> dict:
        return {
            "short_name": cls.short_name(),
            "long_name": cls.long_name(),
            "standard_name": cls.standard_name(alias=alias),
            "unit": str(cls.unit()),
        }
