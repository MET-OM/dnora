from dnora.dnora_type_manager.dnora_types import DnoraDataType


def add_export_method(obj_type: DnoraDataType):
    def wrapper(c):
        def export(self, **kwargs) -> None:
            self.export(obj_type=obj_type, **kwargs)

        exec(f"c.export_{obj_type.name.lower()} = export")
        return c

    return wrapper
