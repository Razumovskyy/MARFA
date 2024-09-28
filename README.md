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
### Note: 
Continuum absorption is not accounted in this project. This functionality might be later added through the effort from new contributors.

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
### Install python packages
In python virtual environment (recommended) run:
```
pip install -r requirements.txt
```
### Choose the atmospheric file
For a quick start you can choose one of the default atmospheric profiles located in the `data/Atmospheres` folder. For example `data/Atmospheres/VenusCO2.dat` which reflects the carbon dioxide profile in the Venus nightside atmosphere.
### Run the project with some command line parameters
For example:
```
fpm run marfa -- CO2 4000 4100 125 tonkov VAC VenusCO2.dat
```

Here is a breakdown of the command-line arguments:
- **`CO2`**: The input molecule
- **`4000`** and **`4100`**: The boundaries of the spectral interval of interest (in cm<sup>-1</sup>)
- **`125`**: The line cut-off condition (in cm<sup>-1</sup>)
- **`tonkov`**: The name of the chi-factor correction used for CO₂.
- **`VAC`** Specifies the target calculated value as volume absorption coefficient.
- **`VenusCO2.dat`** The atmospheric profile file to read pressure, tempareature and molecular density data from.

After running this command, the PT-tables for each level from the `VenusCO2.dat` file will be generated in the `output/ptTables` folder. The output files are created in binary format (direct access files) to facilitate faster integration with radiative transfer models.

### Converting to a human-readable output and plotting
To convert a specific PT-table file to a human-readable format and plot the spectra, use the Python script located in the scripts directory. Execute the following command:
```
python scripts/postprocess.py --v1 4000 --v2 4100 --level 40 --resolution medium --plot
```
You can find the file containing human-readable data in the `processedData` folder, named according to the corresponding format, e.g., `output/processedData/CO2_40_VAC_4000-4100.dat`. Below is an example of the file’s content, which includes log information and data: the first column represents wavenumbers [cm<sup>-1</sup>], and the second column shows the log<sub>10</sub> of the volume absorption coefficient (or cross-section if `targetValue` is set to VAC).
```ini
# Input Molecule: CO2
# Cut Off: 125
# Chi Factor Function Name: tonkov
# Target Value: VAC
# Atmospheric Profile File: VenusCO2.dat
# V1: 4000.0 cm-1
# V2: 4100.0 cm-1
# Resolution: medium
# Level Number: 40
               1       20470
     4000.00000        -1.4572563
     4000.00488        -1.4536923
     4000.00977        -1.4693524
     4000.01465        -1.4902035
      ...
     4099.98535        -5.6071868
     4099.99023        -5.5950541
```
If the `--plot` flag is enabled, a plot of the data set is generated and saved to the `plots` directory with the same name, for example: `output/plots/CO2_40_VAC_4000-4100.png`. Here’s an example:

![](images/CO2_40_VAC_4000-4100.png)

The `V1` and `V2` values do not necessarily need to match the initial boundaries of the spectral interval used to calculate the PT-table. Instead, you can examine a narrower interval with higher resolution to gain more detailed insights:
```
python scripts/postprocess.py --v1 4020 --v2 4022 --level 40 --resolution high --plot
```

![](images/CO2_40_VAC_4020-4022.png)

## Command line parameters: overview
## Atmospheric profile file structure
To correctly run the MARFA code, the atmospheric file must adhere to a specific format and be placed in the `data/Atmospheres/` directory. Below is an example of the required format:

```ini
# Atmospheric file example (Haus2015) CO2                                                                     
81
     0.000 0.90918E+02  733.00  0.8694E+26
     2.000 0.80059E+02  717.00  0.7851E+26
     4.000 0.70286E+02  701.00  0.7075E+26
     6.000 0.61599E+02  685.00  0.6366E+26
```
### File Format Breakdown

1. **Header (First Line):**
   - Describes the atmospheric characteristics, such as planet, authors, and gas species.
   - Example: `# Atmospheric file example (Haus2015) CO2`

2. **Number of Levels (Second Line):**
   - Contains a single number `N` representing the number of atmospheric levels.
   - Example: `81` (indicating 81 levels).

3. **Atmospheric Data (Next N Lines):**
   - Each line contains data in four columns:
     - **Column 1:** Height [km]
     - **Column 2:** Total pressure [atm]
     - **Column 3:** Temperature [K]
     - **Column 4:** Number density [mol/(cm<sup>2</sup>*km)]

   - Example (Venus lower atmosphere):
     ``` 
      0.000 0.90918E+02  733.00  0.8694E+26
      2.000 0.80059E+02  717.00  0.7851E+26
      4.000 0.70286E+02  701.00  0.7075E+26
      6.000 0.61599E+02  685.00  0.6366E+26
     ```

### Important Notes

- **Units:**
  - The data uses [km] for height to ensure the total absorption coefficient is in [km<sup>-1</sup>].

- **Density Calculation:**
  - It is recommended to provide the number density in [mol/(cm<sup>2</sup>*km)] rather than [mol/cm<sup>3</sup>]. This is because the density is also used in calculating the partial pressure of species (for Lorentz HWHM), where the units are fixed. Refer to the `pSelf` calculation in the `app/main.f90` file for more details.

- **Directory Placement:**
  - Ensure the atmospheric file is stored in the `data/Atmospheres/` directory to be recognized by the MARFA code.

## Output PT-table file structure
## Spectral databases
## χ-factors
There is an indication that the far wings of spectral lines tend to diverge from expected Lorentzian or Voigt behavior. To address that χ&#8204;-factor could be applied. 

In MARFA code currently there are several χ&#8204;-factors implemented, which describe sub-Lorentzian behavior of CO<sub>2</sub> rotational lines. They could be found in the `ChiFactors` module. These corrections are primarilly used for **Venus atmosphere conditions**. Below is a reference to the scientific papers from which the analytical expressions were derived (see References section for more details):
| χ-factor | Reference                  |
|----------|----------------------------|
| `tonkov` | Tonkov et al. (1996)       |
| `pollack`| Pollack et al. (1993)      |
| `perrin` | Perrin and Hartmann (1989) |

The χ-factors dataset is intended to be expanded through the effort from other contributors.
## Performance overview
## Introducing custom features
## Troubleshooting

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## References
- Fomin, B. A. _Effective interpolation technique for line-by-line calculations of radiation absorption in gases._ Journal of Quantitative Spectroscopy and Radiative Transfer 53.6 (1995): 663-669.
- Tonkov, M. V., et al. _Measurements and empirical modeling of pure CO<sub>2</sub> absorption in the 2.3-μm region at room temperature: far wings, allowed and collision-induced bands._ Applied optics 35.24 (1996): 4863-4870.
- Pollack, James B., et al. _Near-infrared light from Venus' nightside: A spectroscopic analysis._ Icarus 103.1 (1993): 1-42.
- Perrin, M. Y., and J. M. Hartmann. _Temperature-dependent measurements and modeling of absorption by CO<sub>2</sub>-N<sub>2</sub> mixtures in the far line-wings of the 4.3 μm CO<sub>2</sub> band._ Journal of Quantitative Spectroscopy and Radiative Transfer 42.4 (1989): 311-317.
