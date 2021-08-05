def print_line(length = 75):
    print("-" * length)

def to_file(filename):
    info(f"Writing to file {filename}")
    
def plain(msg):
    print(msg)

def info(msg):
    print(f"*** {msg} ***")

def advice(msg):
    print(f"!!! {msg} !!!")

def header(msg):
    print_line()
    plain(msg)
    print_line()
    
def templates(code):
    if code == 'no_spacing':
        info("No information about grid spacing exists")
        advice("First set spacing with .set_spacing(dlon, dlat) or .set_spacing(dm) (in metres)")
    elif code == 'no_topo':
        info("No topography exists")
        advice("First load topography with .import_topo(topography_fetcher)")
    elif code == 'no_mask':
        info("No land-sea mask exists")
        advice("First load topography with .import_topo(topography_fetcher)")

