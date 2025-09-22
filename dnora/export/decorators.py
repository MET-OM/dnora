from dnora.type_manager.dnora_types import DnoraDataType


def add_export_method(obj_type: DnoraDataType):
    def wrapper(c):
        def export(self, *args, **kwargs) -> None:
            self.export(obj_type, *args, **kwargs)

            for __, single_nest in self._nest.items():
                exec(f"single_nest.export_{obj_type.name.lower()}(*args, **kwargs)")

        exec(f"c.export_{obj_type.name.lower()} = export")
        return c

    return wrapper
