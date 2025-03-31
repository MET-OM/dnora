def read_ww3_nml(filename: str) -> dict:
    """Reads a WAVEWATCH III namelist file (FORTRAN style namelist)"""
    with open(filename, "r") as file:
        nml_dict = {}
        line = "!"
        while line:
            line = file.readline()
            if len(line) < 1:
                continue
            if line[0] == "&":
                main_key = line.replace("\n", "").replace("&", "")
                nml_dict[main_key] = {}
                line = file.readline()
                while line[0] != "/":
                    key, value = line.replace(" ", "").split("=")
                    key, subkey = key.split("%")
                    if nml_dict[main_key].get(key) is None:
                        nml_dict[main_key][key] = {}
                    nml_dict[main_key][key][subkey] = value.replace("\n", "")
                    line = file.readline()
    return nml_dict
