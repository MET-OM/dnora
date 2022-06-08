from pathlib import Path
import yaml
from dataclasses import dataclass
import pandas as pd
import re
from . import msg

def add_prefix(filename: str, prefix: str) -> str:
    """Adds a prefix to a filename, e.g. FileName.txt -> new_FileName.txt"""
    if prefix == '':
        return filename

    if prefix[-1] == '_':
        prefix = prefix[:-1]

    if filename == '':
            return prefix

    if filename[0] == '_':
        filename = filename[1:]

    return f'{prefix}_{filename}'

def add_suffix(filename: str, suffix: str) -> str:
    """Adds a suffix to a filename, e.g. FileName.txt -> FileName_new.txt"""
    if suffix == '':
        return filename

    if suffix[0] == '_':
        suffix = suffix[1:]

    if filename == '':
        return suffix

    filename_list = filename.split('.')

    if len(filename_list) == 1:
        return f'{filename}_{suffix}'
    else:
        filename = '.'.join(filename_list[0:-1])
        extension = f'{filename_list[-1]}'
        return f'{filename}_{suffix}.{extension}'

def replace_times(filename: str, dateformat: str, times: list) -> str:
    """Substitutes the strings #T0, #T1, #T2... etc. in filname with time
    stamps from a list of times, using the format in dateformat.

    e.g. #T0_file.txt, ['2020-05-04 18:00'], %Y%m%d%H%M -> 202005041800_file.txt
    """

    for ct, t in enumerate(times):
        filename = re.sub(f"#T{ct}", pd.Timestamp(t).strftime(dateformat), filename)

    return filename

def replace_lonlat(filename: str, lon: float, lat: float) -> str:
    """Substitutes the strings #Lon, #Lat in filename with values of lon and
    lat.

    e.g. #Lon_#Lat_file.txt, 8.0, 60.05 -> 08.0000000_60.05000000_file.txt
    """

    filename = re.sub("#Lon", f"{lon:010.7f}", filename)
    filename = re.sub("#Lat", f"{lat:010.7f}", filename)

    return filename

def replace_objects(filename: str, dict_of_object_names: dict[str: str]) -> str:
    """Substitutes the strings #{Object} in filename with the name given to
    the object.

    e.g. #Grid_#Forcing_file.txt, [Grid(..., name="Sula"), Forcing(..., name='NORA3')]
        -> Sula_NORA3_file.txt
    """

    for obj_type, obj_name in dict_of_object_names.items():
        if obj_name is not None:
            filename = re.sub(f"#{obj_type}", obj_name, filename)

    return filename

def get_list_of_placeholders():
    defaults_file = Path(__file__).parent.joinpath(Path('defaults.yml'))
    with open(defaults_file, 'r') as file:
      defaults = yaml.safe_load(file)
    return defaults['list_of_placeholders']

def clean(filename: str, list_of_placeholders: list[str]=None) -> str:
    """ Cleans out the file name from possible used placeholders, e.g. #Grid
    as given in the list.

    Also removes multiple underscores '___' etc.
    """
    if list_of_placeholders is None:
        list_of_placeholders = get_list_of_placeholders()

    for s in list_of_placeholders:
            filename = re.sub(s, '', filename)

    filename = re.sub("_{2,10}", '_', filename)
    filename = re.sub("_-_", '', filename)
    filename = re.sub("_-", '_', filename)
    if filename and filename[-1] == '_':
        filename = filename[:-1]

    return filename

def get_default_value(key: str, dnora_obj: str, primary: dict, fallback: dict):
    """Get a key (e.g. folder) from the defaults list.

    1) Tries Model+dnora_obj specific value (e.g. SWAN-wnd-folder)

    2) Tries Model specifc values (e.g. SWAN-folder)

    3) Returns ModelRun defaults (e.g. ModulRun-wnd-folder)
    """

    dnora_obj = dnora_obj.lower()
    if dnora_obj not in fallback.keys():
        raise ValueError(f'Default values not defined for {dnora_obj}!')
    fallback_filename = fallback[dnora_obj].get(key)

    # Try dnora_obj specific filename is dnora_obj settings defined
    dnora_obj_filename = None
    if primary.get(dnora_obj) is not None:
        dnora_obj_filename = primary[dnora_obj].get(key) or dnora_obj_filename

    # If filename not defined for specific dnora_obj, try Model specific name
    dnora_obj_filename = dnora_obj_filename or primary.get(key)

    return dnora_obj_filename or fallback_filename

def add_folder_to_filename(filename: str, folder: str) -> str:
    return str(Path(folder).joinpath(filename))

def split_filepath(filepath: str) -> tuple[str, str]:
    folder = str(Path(filepath).parent)
    filename = Path(filepath).name
    return filename, folder

@dataclass
class FileNames:
    format: str
    dnora_obj: str
    clean_names: bool
    dict_of_object_names: list
    start_time: str
    end_time: str
    _filename: str
    _folder: str
    _dateformat: str
    extension: str

    def __post_init__(self):
        defaults_file = Path(__file__).parent.joinpath(Path('defaults.yml'))
        with open(defaults_file, 'r') as file:
          self._defaults = yaml.safe_load(file)
        self.fallback = self._defaults['ModelRun']
        self.primary = self._defaults[self.format]
        #self.placeholders = self._defaults['list_of_placeholders']

    def dateformat(self) -> str:
        return self._dateformat or get_default_value('dateformat', self.dnora_obj, self.primary, self.fallback)

    def filename(self, extension: str=None) -> str:
        filename = self._filename or get_default_value('filename', self.dnora_obj, self.primary, self.fallback)
        filename = self.replace_placeholders(filename, self.dateformat())
        extension = extension or self.extension
        return Path(filename).with_suffix(f'.{extension}')

    def folder(self) -> str:
        folder = self._folder or get_default_value('folder', self.dnora_obj, self.primary, self.fallback)
        return Path(self.replace_placeholders(folder, self.dateformat()))

    def filepath(self) -> str:
        return add_folder_to_filename(self.filename(), self.folder())

    def create_folder(self) -> None:
        folder = Path(self.folder())
        if not folder.is_dir():
            msg.plain(f"Creating folder {str(folder)}")
            folder.mkdir()

    def replace_placeholders(self, unclean_string: str, dateformat: str) -> str:
        unclean_string = replace_objects(unclean_string, self.dict_of_object_names)

        clean_string = replace_times(unclean_string, dateformat, [self.start_time, self.end_time])

        if self.clean_names:
            clean_string = clean(clean_string)

        return clean_string
