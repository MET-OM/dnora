from dnora.utils.time import create_time_stamps
from dataclasses import dataclass
from dnora.utils.time import create_monthly_stamps


@dataclass
class FileStructure:
    """stride:
    int in hours (e.g. 24 for daily files)
    'month' for monthly files
    None if all data is in one file, but points are spread across different files
    """

    stride: int | str | None
    hours_per_file: int = None
    lead_time: int = 0
    last_file: str = ""
    offset: int = 0
    tile: str = None

    def create_time_stamps(self, start_time: str, end_time: str):
        if self.stride == "month":
            start_times, end_times = create_monthly_stamps(start_time, end_time)
            file_times = start_times
        else:
            start_times, end_times, file_times = create_time_stamps(
                start_time,
                end_time,
                stride=self.stride,
                hours_per_file=self.hours_per_file,
                last_file=self.last_file,
                lead_time=self.lead_time,
                offset=self.offset,
            )

        return start_times, end_times, file_times
