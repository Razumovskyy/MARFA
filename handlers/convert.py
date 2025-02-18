"""
Module for converting binary PT-table files into human-readable absorption data.
"""

import shutil
from pathlib import Path

from handlers.constants import HUMAN_READABLE_DIRECTORY, INFO_FILENAME, get_pt_name
from handlers.parsers import base_parser


def convert_pttable(directory: Path, level: int) -> None:
    """
    Converts a PT-table binary file inside a given directory into a human-readable format.

    It reads the start and end wavenumbers from an `info.txt` file in the same directory, extracts 
    the relevant absorption data, and formats it into a readable table.

    Args:
        directory (Path): The directory containing the PT-table binary file 
                                and the associated file.
        level (int): The atmospheric level at which absorption data should be extracted.
    """

    v1, v2 = None, None

    pttable_file = get_pt_name(directory, level)
    info_file = directory / INFO_FILENAME

    with open(info_file) as info:
        for line in info.readlines():
            if line.startswith("Start Wavenumber"):
                v1 = float(line.split(':')[1])
            if line.startswith("End Wavenumber"):
                v2 = float(line.split(':')[1])
            if line.startswith("Input Molecule"):
                molecule = line.split(':')[1]

    if v1 is None or v2 is None:
        raise ValueError(f"Corrupted info file: start or end wavenumbers are not defined")

    output_filename = f"{molecule.strip()}_{int(v1)}-{int(v2)}_{level}level.dat"
    output_file = shutil.copyfile(info_file, Path(HUMAN_READABLE_DIRECTORY)/output_filename)

    vw_data, absorption_data = base_parser(pttable_file, v1, v2)

    with open(output_file, 'a') as output:
        for vw, abs_data in zip(vw_data, absorption_data):
            output.write(f"{vw:15.5f} {abs_data:17.7e}\n")

    print("Success!")
    print(f"File with absorption data in human-readable format for level {level} "
          f"is saved as {HUMAN_READABLE_DIRECTORY}/{output_filename}")
