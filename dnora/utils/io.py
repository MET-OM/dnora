def set_nested_value(dictionary, keys, value):
    for key in keys[:-1]:  # Iterate through all keys except the last one
        dictionary = dictionary[key]
    dictionary[keys[-1]] = value  # Set the value for the last key


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
                    keys, value = line.replace(" ", "").split("=")
                    keys = keys.split("%")
                    dictionary = nml_dict[main_key]
                    for key in keys[:-1]:
                        if dictionary.get(key) is None:
                            dictionary[key] = {}
                        dictionary = dictionary[key]

                    value = value.replace("\n", "")
                    set_nested_value(nml_dict[main_key], keys, value)
                    line = file.readline()

    return nml_dict
