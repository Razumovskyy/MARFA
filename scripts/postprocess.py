import numpy as np
import sys
import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def process_data(subdir_path, V1, V2, level, resolution):
    print('Processing started')
    NT = 20481
    extention = 'ptbin'
    full_subdir_path = os.path.join('output', 'ptTables', subdir_path)

    # Verify the directory exists
    if not os.path.isdir(full_subdir_path):
        sys.exit(f"Error: Directory {full_subdir_path} does not exist.")

    # Path to the info.txt file
    info_file_path = os.path.join(full_subdir_path, 'info.txt')
    
    # Read info.txt content
    if os.path.isfile(info_file_path):
        with open(info_file_path, 'r') as info_file:
            header_content = info_file.read()
    else:
        sys.exit(f"Error: info.txt file not found in {full_subdir_path}")

    # Extract required information from info.txt
    molecule = None
    target_value = None
    atmospheric_file = None
    for line in header_content.splitlines():
        if "Input Molecule:" in line:
            molecule = line.split("Input Molecule:")[1].strip()
        elif "Target Value:" in line:
            target_value = line.split("Target Value:")[1].strip()
        elif "Atmospheric Profile File:" in line:
            atmospheric_file = line.split("Atmospheric Profile File:")[1].strip()
    
    if not molecule or not target_value:
        sys.exit("Error: Unable to extract 'Input Molecule' or 'Target Value' from 'info.txt'.")

    deltaWV = 10 # deltaWV parameter from marfa code

    # Determine deltaWV based on resolution
    if resolution == 'high':
        granularity = 1  # Fine granularity
    elif resolution == 'medium':
        granularity = 10  # Moderate granularity
    elif resolution == 'coarse':
        granularity = 100 # Broad granularity
    else:
        sys.exit(f"Invalid resolution: {resolution}")

    step = granularity * deltaWV / (NT-1)

    # Calculate NZ1 based on V1
    NZ1 = int(V1 / 10.0)

    # Calculate NZ2 based on V2
    NZ2 = int((V2-1) / 10.0)
    count = NZ2 - NZ1 + 1  # Number of records to process

        # Construct DIR_NAME based on level
    NAME = '___.'  # Initialize with underscores and a dot
    if level < 10:
        NAME = f"{level}{NAME[1:]}"  # e.g., '1___.'
    elif level < 100:
        NAME = f"{level}{NAME[2:]}"  # e.g., '12_.'
    else:
        NAME = f"{level}{NAME[3:]}"  # e.g., '123.'

    # Open the binary file
    filename = os.path.join(full_subdir_path, f"{NAME}{extention}")
    if not os.path.exists(filename):
        sys.exit(f"File {filename} does not exist.")

    # Read the binary file
    record_length = NT * 4  # Each real number is 4 bytes (float32)
    data = []

    record_number = NZ1
    with open(filename, 'rb') as f:
        for II in range(1, count + 1):  # II from 1 to count inclusive
            seek_position = (record_number-1) * record_length
            print(f"Processing record {record_number} ...")

            if seek_position >= os.path.getsize(filename):
                sys.exit(f"Record number {record_number} exceeds file size.")

            f.seek(seek_position)
            record_bytes = f.read(record_length)
            if len(record_bytes) < record_length:
                sys.exit(f"Unexpected end of file at record {record_number}.")

            # Unpack NT float32 numbers using little-endian
            RK = np.frombuffer(record_bytes, dtype='<f4')
            RK = RK.astype(np.float32)

            # Calculate base wavenumber for this record
            V11 = V1 + deltaWV * (II - 1)

            for I in range(1, int((NT+1) / granularity)):
                VV = V11 + (I - 1) * step
                RK_index = (I - 1) * granularity
                if RK_index >= NT:
                    break  # Prevent index out of range
                RK_value = RK[RK_index]
                if RK_value > 0:
                    log_RK = np.log10(RK_value).astype(np.float32)
                else:
                    log_RK = float('-inf')  # Handle log10(0) case
                data.append((VV, log_RK))
            
            record_number = NZ1 + II
    
    # Create DataFrame
    df = pd.DataFrame(data, columns=['Wavenumber', 'Log Absorption Coefficient'])

    # Write to file with formatted output
    formatted_filename = os.path.join('output', 'processedData', f'{molecule}_{level}_{target_value}_{int(V1)}-{int(V2)}.dat')
    with open(formatted_filename, 'w') as f_out:
        # Write header lines
        # Filter out unwanted lines from info.txt
        for line in header_content.splitlines():
            if not any(keyword in line for keyword in ["UUID", "Start Wavenumber", "End Wavenumber", "Command-Line Arguments"]):
                f_out.write(f"# {line}\n")
        
        # Add new header lines
        additional_headers = [
            f"V1: {V1:.1f} cm-1",
            f"V2: {V2:.1f} cm-1",
            f"Resolution: {resolution}",
            f"Level Number: {level}"
        ]
        for header in additional_headers:
            f_out.write(f"# {header}\n")
        M = 1
        total_points = len(data)
        # Match Fortran formatting in the header
        f_out.write(f"{M:>16d}{total_points:>12d}\n")
        for VV, log_RK in data:
            # Write with controlled precision: 5 decimal places for VV and 7 for log_RK
            f_out.write(f"{VV:15.5f} {log_RK:17.7f}\n")

    print("Data processing complete. Output written to " + formatted_filename)
    return df, formatted_filename, atmospheric_file

def plot_spectra(df, file_name, atmospheric_file, Vleft, Vright):
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['text.usetex'] = False
    sns.set(style="whitegrid", context='talk')

    level = int(file_name.split('_')[1])

    # read data from the atmospheric file to match pressure, temperature and level
    atmospheric_file_fullname = os.path.join('data', 'atmospheres', atmospheric_file)

    # Read info.txt content
    if os.path.isfile(atmospheric_file_fullname):
        with open(atmospheric_file_fullname, 'r') as atmosphere:
            lines = atmosphere.readlines()
    else:
        sys.exit(f"Error: Atmospheric file {atmospheric_file_fullname} not found")

    # Ensure the file has at least two lines (title and N)
    if len(lines) < 2:
        sys.exit("Error: Atmospheric file does not contain enough header information.")

    # Extract the title and number of levels
    title = lines[0].strip()
    try:
        N = int(lines[1].strip())
    except ValueError:
        sys.exit("Error: The second line of the atmospheric file should be an integer representing the number of levels.")

    # Validate the level
    if level < 1 or level > N:
        sys.exit(f"Error: Level {level} is out of range. The file contains {N} levels.")

    # Ensure the file has enough data lines
    expected_total_lines = 2 + N  # 1 title + 1 N + N data lines
    if len(lines) < expected_total_lines:
        sys.exit(f"Error: Expected {N} data lines, but found {len(lines) - 2}.")

    # Extract the data line for the specified level
    # Assuming levels are 1-based
    data_line_index = 2 + (level - 1)
    data_line = lines[data_line_index].strip()

    # Split the data line into columns
    columns = data_line.split()

    # Ensure there are enough columns (at least 3 for pressure and temperature)
    if len(columns) < 3:
        sys.exit(f"Error: Data line for level {level} does not contain enough columns.")

    # Extract pressure and temperature
    try:
        height = float(columns[0])
        pressure = float(columns[1])      # Second column: Pressure
        temperature = float(columns[2])   # Third column: Temperature
    except ValueError:
        sys.exit(f"Error: Non-numeric data found in pressure or temperature columns for level {level}.")


    # Create plots directory if it doesn't exist
    plots_dir = os.path.dirname(file_name).replace('processedData', 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    # Define the plot filename based on the data filename
    base_filename = os.path.basename(file_name)
    molecule = base_filename[:3]
    plot_image_name = base_filename.replace('.dat', '.png')
    plot_image_path = os.path.join(plots_dir, plot_image_name)

    # Plot using Seaborn
    plt.figure(figsize=(12, 6))
    ax = sns.lineplot(x='Wavenumber', y='Log Absorption Coefficient', data=df, color='b')

    ax.set_xlabel('Wavenumber [cm$^{-1}$]')
    ax.set_ylabel('Log$_{10}$(Absorption Cross-Section) [cm$^2$ mol$^{-1}$]')
    ax.set_title(f'{molecule} Absorption Spectrum for the atmosphere: {atmospheric_file}\n'
                 f'Level {level}: height {height}[km], pressure {pressure:.3E} [atm], temperature {temperature} [K]')
    
    ax.set_xlim(Vleft, Vright)
    ax.grid(True)

    # Customize tick parameters
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    plt.tight_layout()
    try:
        plt.savefig(plot_image_path)
        print(f"Plot saved to '{plot_image_path}'.")
    except Exception as e:
        sys.exit(f"Error saving plot to '{plot_image_path}': {e}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Process and plot PT-table data.')
    parser.add_argument('subdir', type=str, nargs='?', help="Subdirectory name (e.g., CO2_7500-7600_20240922_015244). If omitted, the latest run will be used.")
    parser.add_argument('--v1', type=float, required=True, help='Starting wavenumber V1')
    parser.add_argument('--v2', type=float, required=True, help='Ending wavenumber V2')
    parser.add_argument('--level', type=int, required=True, help='Atmospheric level')
    parser.add_argument(
        '--resolution',
        type=str,
        required=True,
        choices=['high', 'medium', 'coarse'],
        help="Set resolution for the output data. Choose from 'high', 'medium', or 'coarse'."
    )
    parser.add_argument(
        '--plot',
        action='store_true',
        help="Include this flag to generate and display the plot."
    )
    args = parser.parse_args()

    # Determine subdirectory path
    if args.subdir:
        subdir = args.subdir
    else:
        # Read the latest_run.txt to get the latest subdirectory
        latest_run_file = os.path.join('output', 'ptTables', 'latest_run.txt')
        if not os.path.isfile(latest_run_file):
            sys.exit("Error: No subdirectory provided and latest_run.txt not found.")
        with open(latest_run_file, 'r') as f:
            subdir = f.read().strip()
        if not subdir:
            sys.exit("Error: latest_run.txt is empty.")

    V1 = args.v1
    V2 = args.v2
    level = args.level
    resolution = args.resolution
    plot_flag = args.plot

    # Input validation
    if V1 >= V2:
        sys.exit("Error: V1 must be less than V2.")
    if V1 < 0 or V2 < 0:
        sys.exit("Error: Wavenumbers must be non-negative.")
    if level <= 0:
        sys.exit("Error: Atmospheric level must be a positive integer.")

    # Process data
    df, formatted_filename, atmospheric_file = process_data(subdir, V1, V2, level, resolution)

    # plt.show()
    
    # Conditionally plot
    if plot_flag:
        plot_spectra(df, formatted_filename, atmospheric_file, V1, V2)
    else:
        print("Plotting skipped as per the '--plot' flag.")

if __name__ == '__main__':
    main()
