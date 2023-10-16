from .general.general_writing_functions import DnoraNc

from ..dnora_object_type import DnoraObjectType


def add_export_method(obj_type: DnoraObjectType):
    def wrapper(c):
        def export(
            self,
            writer: str = None,
            filename: str = None,
            folder: str = None,
            dateformat: str = None,
            format: str = None,
            dry_run=False,
            **kwargs,
        ) -> None:
            # dnora_obj_type = type(
            #     self.model[obj_type]
            # ).__name__  # Can be Grid or TriGrid fro grid export
            # if dnora_obj_type in [
            #     "DummyDnoraObject",
            #     "NoneType",
            # ]:  # Dry running can create Dummy objects, gives None is no data exists.
            #     dnora_obj_type = obj_type
            writer_function = self._setup_export(obj_type, writer, dry_run)

            if not self.dry_run():
                try:  # GeneralWritingFunction might not have this method defined
                    wanted_convention = writer_function.convention()
                except AttributeError:
                    wanted_convention = self._get_spectral_convention()

                try:
                    self.model[obj_type]._set_convention(wanted_convention)
                except AttributeError:  # Can only be done for spectra
                    pass

            self._export_object(
                obj_type,
                filename,
                folder,
                dateformat,
                writer_function=writer_function,
                format=format,
                **kwargs,
            )

        def cache(self) -> None:
            exec(f'self.export_{obj_type.value}(writer=DnoraNc(), format="Cache")')

        exec(f"c.export_{obj_type.value} = export")
        exec(f"c.cache_{obj_type.value} = cache")

        return c

    return wrapper
