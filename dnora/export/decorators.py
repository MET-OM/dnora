from dnora.type_manager.dnora_types import DnoraDataType


def add_export_method(obj_type: DnoraDataType):
    def wrapper(c):
        def export(self, *args, **kwargs) -> None:
            self.export(obj_type, *args, **kwargs)

            nest = self._nest

            while nest is not None:
                nest.export(obj_type, *args, **kwargs)
                nest = nest._nest

        exec(f"c.export_{obj_type.name.lower()} = export")
        return c

    return wrapper
