import numpy as np
import sys
import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def process_data(V1, V2, level, resolution):
    NT = 20481
    MET = 'CO2'
    OUTPUT_DIR = './output/PT_CALC/'

    # Construct DIR_NAME based on level
    DIR_NAME = '___.'  # Initialize with underscores and a dot
    if level < 10:
        DIR_NAME = f"{level}{DIR_NAME[1:]}"  # e.g., '1___.'
    elif level < 100:
        DIR_NAME = f"{level}{DIR_NAME[2:]}"  # e.g., '12_.'
    else:
        DIR_NAME = f"{level}{DIR_NAME[3:]}"  # e.g., '123.'

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

    # Open the binary file
    filename = os.path.join(OUTPUT_DIR, f"{DIR_NAME}{MET}")
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

    # Write to 'SPECTR_PYTHON.dat' with formatted output
    with open('SPECTR_PYTHON.dat', 'w') as f_out:
        M = 1
        total_points = len(data)
        # Match Fortran formatting in the header
        f_out.write(f"{M:>16d}{total_points:>12d}\n")
        for VV, log_RK in data:
            # Write with controlled precision: 5 decimal places for VV and 7 for log_RK
            f_out.write(f"{VV:15.5f} {log_RK:17.7f}\n")

    print("Data processing complete. Output written to 'SPECTR_PYTHON.dat'.")
    return df

def plot_spectra(df):
    sns.set(style="whitegrid", context='talk')

    # Plot using Seaborn
    plt.figure(figsize=(12, 6))
    ax = sns.lineplot(x='Wavenumber', y='Log Absorption Coefficient', data=df, color='b')

    ax.set_xlabel('Wavenumber [cm$^{-1}$]')
    ax.set_ylabel('Log$_{10}$(Absorption Coefficient) [cm$^2$ mol$^{-1}$]')
    ax.set_title('Absorption Spectrum')
    ax.grid(True)

    # Customize tick parameters
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)

    plt.tight_layout()
    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Process and plot PT-table data.')
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
    df = process_data(V1, V2, level, resolution)

    # Conditionally plot
    if plot_flag:
        plot_spectra(df)
    else:
        print("Plotting skipped as per the '--plot' flag.")

if __name__ == '__main__':
    main()