from abc import ABC, abstractstaticmethod
import numpy as np

from pint import Unit


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

    @classmethod
    def find_me_in_ds(cls, ds) -> str | None:
        """Takes an Xarray Dataset and returns name of variable that matches the parameter based on standard_name"""
        data_vars = list(ds.data_vars)

        for var in data_vars:
            if hasattr(ds.get(var), "standard_name"):
                if ds.get(var).standard_name in cls.standard_aliases():
                    return var

        return None
