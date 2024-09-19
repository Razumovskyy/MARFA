# MARFA (Molecular Absorption with Rapid and Flexible Analysis)
## Overview
MARFA is a Fortran-based tool designed for calculating volume molecular absorption coefficients or monochromatic absorption cross-sections based on initial spectroscopic data and a given atmospheric profile. Originally developed to facilitate modeling of radiative transfer in Venus's atmosphere, its flexible design permits adaptation to various spectroscopic and atmospheric scenarios.

In addition to using and contributing to the source code, it is recommended to interact with the web interface of the tool to better understand its capabilities. The web interface can be accessed at: <URL>

## Features

- **Efficient Line-by-Line Technique**: Employs an effective interpolation method [Fomin, 1995] featuring nine grids to speed up line summation in the core of the tool.
- **Atmospheric Profile Handling**: Computes absorption coefficients or cross-sections based on atmospheric profiles provided via input files.
- **Line Shapes Support**: Supports standard line shapes with the Voigt profile set as the default. Additional line shapes can be manually added by implementing a specific abstract interface.
- **Line Wings Corrections**: Supports various chi-factors for implementing sub-Lorentzian behavior of wings of spectral lines. Custom chi-factors can be added manually.
- **Line Cut-Off Criterion**: Can be set as an input to control accuracy and match continuum parameters.
- **Line Databases Support**: Includes HITRAN2016 databases for CO₂ and H₂O within the source code. Custom databases can be introduced by preprocessing them into the required format.
- **PT-Tables Generation**: Generates resulting spectra as direct-access files in PT-format (each output file corresponds to one atmospheric level), which can be immediately introduced into radiative transfer schemes.
- **Additional Tools**: Provides various scripts for plotting and data pre-processing and post-processing, facilitating validation and the integration of new data.

## Prerequisites
To build and run this tool on your machine, you need to have GFortran (GNU Fortran compiler) and the Fortran Package Manager (fpm) installed. The tool is compatible with Windows, macOS, and Linux operating systems.

### Gfortran
For installing the `gfortran` you can use [GNU Fortran website](https://gcc.gnu.org/fortran/) or use your system's package manager.
### fpm
Installation instructions ara available on the [official website](https://fpm.fortran-lang.org/install/index.html) or on the [fpm github page](https://github.com/fortran-lang/fpm).

## Quick start instructions
### Clone the repository:

```
git clone https://github.com/Razumovskyy/MultiGrid-AtmoSpectra
```
### Navigate to the project directory:
for example:
```
cd MultiGrid-AtmoSpectra
```
### Build the project
```
fpm build
```
### Set the atmospheric file
For a quick start you can choose one of the default atmospheric profiles located in the `data/Atmospheres` folder. For example `CO2_gas_profile.dat` which reflects the carbon dioxide profile in the Venus nightside atmosphere.
### Run the project with some command line parameters
For example:
```
fpm run VenusPT-Tables -- CO2 4000 4100 125 pollack VAC CO2_gas_profile.dat default
```

Here is a breakdown of the command-line arguments:
- **`CO2`**: The input molecule
- **`4000`** and **`4100`**: The boundaries of the desired spectral interval (in cm<sup>-1</sup>)
- **`125`**: The line cut-off condition (in cm<sup>-1</sup>)
- **`pollack`**: The name of the chi-factor correction used for CO₂.
- **`VAC`** Specifies the target value as volume absorption coefficient.
- **`CO2_gas_profile.dat`** The atmospheric profile file to use.
- **`default`** - A technical flag indicating that this run is made locally rather than via the web version.

After running this command, the PT-tables for each level from the `CO2_gas_profile.dat` file would be genereated in the `output` folder. 
## Command line parameters
## Modules overview
### `MolarMasses.f90`
Contains molar masses of 124 isotoplogues of 42 gaseous species common for spectroscopic and atmospheric studies.
### `Constants.f90`
This module contains the main physical constants and other unchangeable figures used throughout the project. The chosen unit system is cgs (centimeter-gram-second). To ensure high numerical accuracy in line-by-line calculations, double-precision floating-point numbers are used for representing wavenumber related values. The precision and range of real numbers can be specified using the kind parameter `DP` which could be found in the `Constants.f90` module.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## References
- Fomin, B.A. (1995). _Effective interpolation technique for line-by-line calculations of radiation absorption in gases_. Journal of Quantitative Spectroscopy and Radiative Transfer, 53(6), 663–669.
