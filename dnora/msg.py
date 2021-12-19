def print_line(marker='-', length=75):
    """Use to print a line '-----...-----'"""

    print(marker * length)

def to_file(filename):
    plain(f"Writing to file >>> {filename}")

def from_file(filename):
    plain(f"Reading from file <<< {filename}")

def plain(msg):
    print(msg)

def info(msg):
    print(f"*** {msg} ***")

def advice(msg):
    print(f"!!! {msg} !!!")

def warning(msg, length=90):
    blank()
    print_line('!', length)
    advice(msg)
    print_line('!', length)
    blank()

def header(obj, msg):
    """Use with objects that are an instance of an abstract class."""

    blank()
    print_line()
    plain(f"{type(obj).__name__} ({obj.__class__.__bases__[0].__name__}): {msg}")
    print_line()

def blank():
    print('')

def process(msg):
    print(f">>> {msg} <<<")

def templates(code):
    if code == 'no_spacing':
        info("No information about grid spacing exists")
        advice("First set spacing with .set_spacing(dlon, dlat) or .set_spacing(dm) (in metres)")
    elif code == 'no_topo':
        info("No topography exists")
        advice("First load topography with .import_topo(topography_reader)")
    elif code == 'no_mask':
        info("No land-sea mask exists")
        advice("First load topography with .import_topo(topography_reader)")
