"""
Authors:
    Mikhail Razumovskii and Denis Astanin, 2025

Description:
    This module is a part of the MARFA-webapp project.
"""
import os
import sys
from pathlib import Path

import numpy as np

from handlers.constants import INFO_FILENAME, POINTS_PER_RECORD, RECORD_SIZE, RECORD_WV_SPAN, get_pt_name


def base_parser(pttable_file: Path, v1: float, v2: float) -> tuple[list[np.float32], list[np.float32]]:
    """
    Parse data from a pt-table binary file based on left and right spectral boundaries.

    This function reads a binary pt-table file, extracts absorption data within the
    spectral range defined by v1 and v2, and calculates corresponding wavenumbers.
    It processes records sequentially within the specified boundaries and returns
    two lists: one for the calculated wavenumbers and one for the absorption data.

    Args:
        pttable_file (Path): Path to the binary pt-table file containing absorption data.
        v1 (float): Left spectral boundary (starting value).
        v2 (float): Right spectral boundary (ending value).

    Returns:
        tuple[list[np.float32], list[np.float32]]:
            - First element: List of wavenumber values (vw_data) as np.float32.
            - Second element: List of corresponding absorption values as np.float32.
    """
    start_record_number = int(v1 / 10.0)
    end_record_number = int((v2 - 1) / 10.0)
    num_records = end_record_number - start_record_number + 1
    step = RECORD_WV_SPAN / (POINTS_PER_RECORD - 1)

    vw_data = []
    absorption_data = []
    record_number = start_record_number
    try:
        with open(pttable_file, 'rb') as f:
            for j in range(1, num_records + 1):
                seek_position = (record_number - 1) * RECORD_SIZE
                if seek_position >= os.path.getsize(pttable_file):
                    raise IndexError(f"Record number {record_number} exceeds a file size.")
                f.seek(seek_position)
                binary_abs_data = f.read(RECORD_SIZE)  # reads absorption data from one record

                abs_data = np.frombuffer(binary_abs_data, dtype=np.float32)

                in_record_start_wv = RECORD_WV_SPAN * record_number
                for i in range(POINTS_PER_RECORD):
                    vw = in_record_start_wv + i * step
                    vw_data.append(np.float32(vw))
                    absorption_data.append(np.float32(abs_data[i]))

                record_number = start_record_number + j
    except FileNotFoundError:
        print(f"Error: invalid number of atmospheric level")
        sys.exit()
    return vw_data, absorption_data


def plot_parser(directory: Path, level: int, vl: float, vr: float) -> tuple[
        list[np.float32], list[np.float32], tuple[str, int, str, str, int, str]]:
    """
    Validates input boundaries for plots requests

    Args:
        directory (Path): Path to the directory containing pt-table files.
        level (int): The atmospheric level at which absorption data should be extracted.
        vl (float): left boundary in cm-1
        vr (float): right boundary in cm-1

    Returns:
        tuple:
            list[np.float32]: x-axis data
            list[np.float32]: y-axis data
            tuple[str, int, str, str, int, str]: Plotting parameters including molecule name,
              cutoff, y-title, atmospheric profile filename, level, and plot filename.
    """

    pttable_file = get_pt_name(directory, level)
    info_file = directory / INFO_FILENAME

    with open(info_file) as info:
        for line in info.readlines():
            if line.startswith("Start Wavenumber"):
                v1 = float(line.split(':')[1])
            if line.startswith("End Wavenumber"):
                v2 = float(line.split(':')[1])
            if line.startswith("Input Molecule"):
                molecule = line.split(':')[1].strip()
            if line.startswith("Cut Off"):
                cutoff = int(line.split(':')[1].strip())
            if line.startswith("Target Value"):
                y_title = line.split(':')[1].strip()
            if line.startswith("Atmospheric Profile File"):
                atm_filename = line.split(':')[1].strip()

    if not v1 <= vl < vr <= v2:
        raise ValueError(f"Left and right boundaries requested for plotting: {vl}, {vr} cm-1 are "
                         f"outside the spectre range: {v1}, {v2} cm-1.")
    x_data, y_data = base_parser(pttable_file, vl, vr)
    plot_filename = f"{molecule}_{int(vl)}-{int(vr)}_{y_title}.svg"
    parameters = (molecule, cutoff, y_title, atm_filename, level, plot_filename)
    return x_data, y_data, parameters
