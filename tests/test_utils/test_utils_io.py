import unittest
import sys
import numpy as np
import pandas as pd

from dnora.utils.time import create_time_stamps
from dnora.utils.io import get_url


def test_get_url():
    time = pd.Timestamp("2020-01-01 07:00")
    url = "https://run[6]/%Y-%m-%d %H:%M"
    folder = "const"
    new_url = get_url(folder, url, time)
    assert new_url == "const/https://run06/2020-01-01 07:00"

    time = pd.Timestamp("2020-01-01 07:00")
    url = "https://run[3]/%Y-%m-%d %H:%M"
    folder = "const"
    new_url = get_url(folder, url, time)
    assert new_url == "const/https://run06/2020-01-01 07:00"

    time = pd.Timestamp("2020-01-01 07:00")
    url = "https://run[12]/%Y-%m-%d %H:%M"
    folder = "const"
    new_url = get_url(folder, url, time)
    assert new_url == "const/https://run00/2020-01-01 07:00"

    time = pd.Timestamp("2020-01-01 14:00")
    url = "https://run[12]/%Y-%m-%d %H:%M"
    folder = "const"
    new_url = get_url(folder, url, time)
    assert new_url == "const/https://run12/2020-01-01 14:00"


class TimeStamps(unittest.TestCase):
    def test_stride(self):
        """Test that stride works"""
        start_time = "2020-01-01 00:00"
        end_time = "2020-01-02 06:00"

        st, et, ft = create_time_stamps(start_time, end_time, stride=12)

        self.assertEqual(len(st), 3)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")
        self.assertEqual(st[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 23:00")
        self.assertEqual(et[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 06:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")
        self.assertEqual(ft[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

    def test_lead_time(self):
        """Test that leadtime works"""
        start_time = "2020-01-01 00:00"
        end_time = "2020-01-02 00:00"

        st, et, ft = create_time_stamps(start_time, end_time, stride=12, lead_time=24)

        self.assertEqual(len(st), 3)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")
        self.assertEqual(st[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 23:00")
        self.assertEqual(et[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2019-12-31 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2019-12-31 12:00")
        self.assertEqual(ft[2].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")

    def test_negative_lead_time(self):
        """Test that negative lead_time works"""
        start_time = "2020-01-01 00:00"
        end_time = "2020-01-02 00:00"

        st, et, ft = create_time_stamps(start_time, end_time, stride=12, lead_time=-24)

        self.assertEqual(len(st), 3)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")
        self.assertEqual(st[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 23:00")
        self.assertEqual(et[2].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-02 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2020-01-02 12:00")
        self.assertEqual(ft[2].strftime("%Y-%m-%d %H:%M"), "2020-01-03 00:00")

    def test_last_file(self):
        """Test that last_file works"""
        start_time = "2020-01-01 00:00"
        end_time = "2020-01-02 00:00"

        st, et, ft = create_time_stamps(
            start_time, end_time, stride=12, last_file="2020-01-01 12:00"
        )

        self.assertEqual(len(st), 2)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 23:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

    def test_hours_per_file(self):
        """Test that last_file works"""
        start_time = "2020-01-01 00:00"
        end_time = "2020-01-10 00:00"

        st, et, ft = create_time_stamps(
            start_time,
            end_time,
            stride=12,
            last_file="2020-01-01 12:00",
            hours_per_file=48,
        )

        self.assertEqual(len(st), 2)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-03 11:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

    def test_between_files(self):
        """Test that it works to start and stop between files"""
        start_time = "2020-01-01 02:00"
        end_time = "2020-01-01 17:00"

        st, et, ft = create_time_stamps(start_time, end_time, stride=12)

        self.assertEqual(len(st), 2)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 02:00")
        self.assertEqual(st[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 17:00")

        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")
        self.assertEqual(ft[1].strftime("%Y-%m-%d %H:%M"), "2020-01-01 12:00")

    def test_start_after_last_file(self):
        """Test that it works to start and stop between files"""
        start_time = "2020-01-01 11:00"
        end_time = "2020-01-05 17:00"

        st, et, ft = create_time_stamps(
            start_time,
            end_time,
            stride=12,
            last_file="2020-01-01 00:00",
            hours_per_file=72,
        )

        self.assertEqual(len(st), 1)

        self.assertEqual(st[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 11:00")
        self.assertEqual(et[0].strftime("%Y-%m-%d %H:%M"), "2020-01-03 23:00")
        self.assertEqual(ft[0].strftime("%Y-%m-%d %H:%M"), "2020-01-01 00:00")


if __name__ == "__main__":
    unittest.main()
