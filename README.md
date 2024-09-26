# MARFA
## Overview
MARFA (Molecular Atmospheric Absorption with Rapid and Flexible Analysis) is a versatile tool designed to calculate volume absorption coefficients or monochromatic absorption cross-sections using initial spectroscopic data from sources such as the HITRAN and HITEMP databases. With MARFA, users can generate absorption look-up tables for each atmospheric level simultaneously within seconds. These look-up tables are produced in a binary format, making them easily integrable into radiative transfer codes. Originally developed to facilitate the modeling of radiative transfer in Venus's atmosphere, MARFA's flexible design allows it to be adapted to a wide range of spectroscopic and atmospheric scenarios.

In addition to using and contributing to the source code, it is recommended to interact with the web interface of the tool to better understand its capabilities. The web interface can be accessed at: <URL>

## Features

- **Efficient Line-by-Line Technique**: Utilizes an effective interpolation method [Fomin, 1995] with nine grids to accelerate the summation of contributions from a large number of spectral lines.

- **Atmospheric Profile Handling**: Computes output absorption spectra for all atmospheric levels in a single runtime.

- **Line Shapes Support**: Supports standard line shapes, including Doppler, Lorentz, and Voigt (default). Additional line shapes can be manually added.

- **Line Wings Corrections**: Implements various chi-factors to accommodate sub-Lorentzian behavior of spectral line wings. Currently chi-factors are implemented for CO₂ used for Venus atmosphere. Custom chi-factors can be manually added.

- **Line Cut-Off Criterion**: Allows users to set input parameters to control accuracy and align with continuum parameters.

- **Line Databases Support**: Includes HITRAN2016 databases for CO₂ and H₂O within the source code. Other spectral databases can be incorporated by preprocessing them into the required format.

- **PT-Tables Generation**: Produces resulting spectra as direct-access files in PT-format (each output PT-file corresponds to one atmospheric level), which can be directly integrated into radiative transfer schemes.

- **Additional Tools**: Provides various scripts for plotting and data processing, facilitating validation and the integration of new data.

## Prerequisites
To build and run the source code on your machine, you need to have GFortran (GNU Fortran compiler) and the Fortran Package Manager (fpm) installed. 
The tool is compatible with macOS, Linux and Windows operating systems. However, installation might be challenging for Windows unexperienced users.

### Fortran compiler: Gfortran
For installing the `gfortran` you can use [GNU Fortran website](https://gcc.gnu.org/fortran/) or use your system's package manager.
#### Other compilers:
Other fortran compilers were not tested and checked.
### fpm
Installation instructions ara available on the [official website](https://fpm.fortran-lang.org/install/index.html) or on the [fpm github page](https://github.com/fortran-lang/fpm).
### Python3
Python3 is mainly needed for running the plotting scripts and converting binary files to a human readable format.

## Quick start instructions
### Clone the repository and navigate to the project directory:

```
git clone https://github.com/Razumovskyy/MARFA
cd MARFA
```
### Build the project
```
fpm build
```
### Choose the atmospheric file
For a quick start you can choose one of the default atmospheric profiles located in the `data/Atmospheres` folder. For example `VenusCO2.dat` which reflects the carbon dioxide profile in the Venus nightside atmosphere.
### Run the project with some command line parameters
For example:
```
fpm run VenusPT-Tables -- CO2 4000 4100 125 tonkov VAC VenusCO2.dat
```

Here is a breakdown of the command-line arguments:
- **`CO2`**: The input molecule
- **`4000`** and **`4100`**: The boundaries of the desired spectral interval (in cm<sup>-1</sup>)
- **`125`**: The line cut-off condition (in cm<sup>-1</sup>)
- **`tonkov`**: The name of the chi-factor correction used for CO₂.
- **`VAC`** Specifies the target calculated value as volume absorption coefficient.
- **`VenusCO2.dat`** The atmospheric profile file to read pressure, tempareature and molecular density data from.

After running this command, the PT-tables for each level from the `VenusCO2.dat` file would be genereated in the `output/ptTables` folder. The output files are genereated in the binary format (direct access files), in order to speed-up the integration to radiative transfer models.

### Converting to a human-readable output and plotting
To convert the specific PT-table file to a human-readable format and to plot the spectra, you can use the python script located in the `scripts` directory. 

## Command line parameters: overview
## Atmospheric profile file structure
## Output PT-table file structure
## Chi-factors
## Spectral databases
## Introducing custom features

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## References
- Fomin, B.A. (1995). _Effective interpolation technique for line-by-line calculations of radiation absorption in gases_. Journal of Quantitative Spectroscopy and Radiative Transfer, 53(6), 663–669.
