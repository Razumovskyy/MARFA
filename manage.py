"""
Base module for CLI handling users requests for plotting and files conversion
"""
from pathlib import Path
import sys

BASE_DIR = 'output/ptTables'


def main():
    """
    Handles commands running with  `python manage.py`:

    `python manage.py convert`: CLI-interactive command for converting ptbin PT-files to human-readable formats
    `python manage.py plot`: CLI-interactive command for plotting spectra
    """
    if len(sys.argv) < 2:
        print("Usage: manage.py <command> [options]")
        sys.exit(1)

    command = sys.argv[1]

    latest_run_file = Path(BASE_DIR) / 'latest_run.txt'

    with open(latest_run_file, 'r') as f:
        subdir = f.read().strip()
        if not subdir:
            sys.exit("Error: latest_run.txt is empty.")

    if command == 'convert':
        from handlers.convert import convert_pttable
        level = int(input("Enter the level number: "))
        convert_pttable(Path(BASE_DIR) / subdir, level)
    elif command == 'plot':
        from handlers.plot import generate_plot
        from handlers.parsers import plot_parser
        level = int(input("Enter the level number: "))
        vl = float(input("Enter the left boundary for plot (cm-1): "))
        vr = float(input("Enter the right boundary for plot (cm-1): "))
        x_data, y_data, parameters = plot_parser(Path(BASE_DIR) / subdir, level, vl, vr)
        generate_plot(x_data, y_data, parameters)

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == '__main__':
    main()
