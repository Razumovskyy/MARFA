"""
Module for keeping constants and function for getting a PT-file based on a level number and directory
"""
from pathlib import Path

POINTS_PER_RECORD = 20481
RECORD_WV_SPAN = 10  # one record spans 10 cm-1 of absorption data
RECORD_SIZE = POINTS_PER_RECORD * 4  # 4 bytes per each value
PT_EXTENSION = 'ptbin'
INFO_FILENAME = 'info.txt'
HUMAN_READABLE_DIRECTORY = 'output/processedData'
PLOTS_DIRECTORY = 'output/plots'


def get_pt_name(directory: Path, level: int,) -> Path:
    """
    Generates pt-table file name based on a level number and directory name

    Args:
        directory(Path): Directory containing pt-table files
        level(int): The atmospheric level at which absorption data should be extracted.

    Returns:
        Path: PT-table file name
    """
    name = '___.'
    if level < 10:
        name = f"{level}{name[1:]}"  # e.g., '1___.'
    elif level < 100:
        name = f"{level}{name[2:]}"  # e.g., '12_.'
    else:
        name = f"{level}{name[3:]}"  # e.g., '123.'
    pttable_file = directory / f"{name}{PT_EXTENSION}"
    return pttable_file
