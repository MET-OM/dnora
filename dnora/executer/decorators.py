from dnora.type_manager.dnora_types import DnoraFileType


def add_write_method(file_type: DnoraFileType):
    def wrapper(c):
        def write(self, **kwargs) -> None:
            self._write(file_type=file_type, **kwargs)
            for __, single_nest in self._nest.items():
                exec(f"single_nest.write_{file_type.name.lower()}_file(**kwargs)")

        exec(f"c.write_{file_type.name.lower()}_file = write")
        return c

    return wrapper


def add_run_method(file_type: DnoraFileType):
    def wrapper(c):
        def run(self, **kwargs) -> None:
            self._run(file_type=file_type, **kwargs)
            for __, single_nest in self._nest.items():
                exec(f"single_nest.run_{file_type.name.lower()}(**kwargs)")

        exec(f"c.run_{file_type.name.lower()} = run")
        return c

    return wrapper
