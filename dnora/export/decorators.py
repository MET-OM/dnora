from .generic.generic_writers import DnoraNc

from ..dnora_object_type import DnoraObjectType
from ..model_formats import ModelFormat


def add_export_method(obj_type: DnoraObjectType):
    def wrapper(c):
        def export(self, **kwargs) -> None:
            exec(f"self.export(obj_type={obj_type}, **kwargs)")

        exec(f"c.export_{obj_type.value} = export")
        return c

    return wrapper
