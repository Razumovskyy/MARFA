"""
Module for generating plots
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib
import gc
import numpy as np

from handlers.constants import PLOTS_DIRECTORY


def generate_plot(x_data: list[np.float32], y_data: list[np.float32], parameters: tuple) -> None:
    """
    Generates a plot using matplotlib and saves it in SVG format.
    Adds a representative title to the figure.

    Args:
        x_data(list[np.float32]): X-axis data
        y_data(list[np.float32]): Y-axis data
        parameters(tuple): parameters of calculation used for dynamic render of plot title, axis title, and
                            generating a filename

    Returns:
        None
    """

    matplotlib.use("Agg")
    plt.figure(figsize=(12, 6))
    fig, ax = plt.subplots()
            
    molecule = parameters[0]
    cutoff = parameters[1]
    y_title = parameters[2]
    atm_filename = parameters[3]
    level = parameters[4]
    plot_filename = parameters[5]
    
    if y_title == 'ACS':
        yaxis_title = f'Absorption Cross-Section [cm$^{{2}}$ mol$^{{-1}}$]'
    elif y_title == 'VAC':
        yaxis_title = f'Volume Absorption Coefficient [km$^{{-1}}$]'
    else:
        raise ValueError(f"Unknown target value {y_title}")
    y_data = np.log10(y_data)

    ax.plot(x_data, y_data, color='g', linestyle='-')

    plt.title(f"Absorption spectrum of {molecule} at {level} level of atmosphere {atm_filename}\n "
              f"Cutoff is {cutoff} cm-1")
    ax.set_xlabel(r'Wavenumber [$\mathregular{cm^{-1}}$]')
    ax.set_ylabel(yaxis_title)
    ax.grid(which='major', axis='both', color='gray', alpha=0.5)

    ax.minorticks_on()
    ax.tick_params(axis='y', which='both', right=True)

    plt.savefig(Path(PLOTS_DIRECTORY)/plot_filename, format='svg', bbox_inches='tight')
    plt.close(fig)

    gc.collect()
    
    print("Success!")
    print(f"Plot is saved as {PLOTS_DIRECTORY}/{plot_filename}")
