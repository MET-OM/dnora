from dnora.utils.time import create_time_stamps
from dataclasses import dataclass


@dataclass
class FileStructure:
    stride: int
    hours_per_file: int
    lead_time: int = 0
    last_file: str = ""
    offset: int = 0

    def create_time_stamps(self, start_time: str, end_time: str):
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
