def print_line(marker="-", length=75):
    """Use to print a line '-----...-----'"""

    print(marker * length)

def copy_file(file_from, file_to):
    plain(f"Copying {file_from} >>> {file_to}")

def to_file(filename):
    plain(f"Writing to file >>> {filename}")


def to_multifile(filenames: list[str], max_length: int = 5):
    if len(filenames) <= max_length:
        for file in filenames:
            to_file(file)
        return
    if max_length % 2 == 0:
        max_length -= 1
    inds = list(range(len(filenames)))

    inds[max_length // 2 : -max_length // 2 + 1] = []

    for ind in inds[0 : len(inds) // 2]:
        to_file(filenames[ind])
    plain(f"... (total of {len(filenames)} files) ...")
    for ind in inds[len(inds) // 2 :]:
        to_file(filenames[ind])


def from_file(filename):
    plain(f"Reading from file <<< {filename}")


def from_multifile(filenames: list[str], max_length: int = 5):
    if len(filenames) <= max_length:
        for file in filenames:
            from_file(file)
        return
    if max_length % 2 == 0:
        max_length -= 1
    inds = list(range(len(filenames)))

    inds[max_length // 2 : -max_length // 2 + 1] = []

    for ind in inds[0 : len(inds) // 2]:
        from_file(filenames[ind])
    plain(f"... (total of {len(filenames)} files) ...")
    for ind in inds[len(inds) // 2 :]:
        from_file(filenames[ind])


def plain(msg):
    print(msg)


def info(msg):
    print(f"*** {msg} ***")


def advice(msg):
    print(f"!!! {msg} !!!")


def warning(msg, length=90):
    blank()
    print_line("!", length)
    advice(msg)
    print_line("!", length)
    blank()


def header(obj, msg):
    """Use with objects that are an instance of an abstract class."""

    blank()
    print_line()
    if isinstance(obj, str):
        plain(f"{obj}: {msg}")
    else:
        plain(f"{type(obj).__name__} ({obj.__class__.__bases__[0].__name__}): {msg}")
    print_line()


def blank():
    print("")


def process(msg):
    print(f">>> {msg} <<<")
