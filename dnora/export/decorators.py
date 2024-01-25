from .generic.generic_writers import DnoraNc

from ..dnora_types import DnoraDataType
from ..model_formats import ModelFormat


def add_export_method(obj_type: DnoraDataType):
    def wrapper(c):
        def export(self, **kwargs) -> None:
            exec(f"self.export(obj_type={obj_type}, **kwargs)")

        exec(f"c.export_{obj_type.name.lower()} = export")
        return c

    return wrapper
