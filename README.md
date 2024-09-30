# MARFA
## Table of Contents
- [Introduction](#overview)
- [Features](#features)
- [Usage-scenarious-and-limitations](#usage-scenarios-and-limitations)
- [Prerequisites](#prerequisites)
- [Quick Start Instructions](#quick-start-instructions)
- [Command Line Parameters](#command-line-parameters)
- [Atmospheric Profile File Structure](#atmospheric-profile-file-structure)
- [Output PT-table structure](#output-pt-table-file)
- [Spectral Databases](#spectral-databases)
- [χ-factors](#χ-factors)
- [Other spectroscopic data](#other-spectral-or-molecular-data)
  - [TIPS](#tips)
  - [Molar masses](#molar-masses)
- [Performance Estimations](#performance-estimations)
- [Introducing Custom Features](#introducing-custom-features)
    - [Custom χ-factors](#custom-χ-factors)
    - [Custom line shapes](#custom-line-shapes-for-advanced-users)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [References](#references)

## Overview
MARFA (Molecular Atmospheric Absorption with Rapid and Flexible Analysis) is a versatile tool designed to calculate volume absorption coefficients or monochromatic absorption cross-sections using initial spectroscopic data from sources such as the HITRAN and HITEMP databases. With MARFA, users can generate absorption PT-tables (look-up tables) for each atmospheric level simultaneously within seconds. These PT-tables are produced in a binary format, making them easily integrable into radiative transfer codes. Originally developed to facilitate the modeling of radiative transfer in Venus's atmosphere, MARFA's flexible design allows it to be adapted to a wide range of spectroscopic and atmospheric scenarios.

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

## Usage Scenarios and Limitations

The codes are well-suited for calculating absorption features in scenarios where spectral and atmospheric data are uncertain, such as for terrestrial planets or exoplanets. They can efficiently handle long cut-off conditions without a significant increase in computational time. Additionally, the codes are designed to allow users to introduce their own functionality. Therefore, if you have your own functions but lack a computational core, you can use these codes to calculate absorption features over a wide range of input values.

These codes can be applied to calculate absorption features across the far, mid, and near infrared (IR) regions, as well as the visible spectrum, roughly covering the range from 10 to 20,000 cm<sup>-1</sup>. However, you must verify which spectral range is covered by the database you are using. 

The microwave region can also be considered, but certain adjustments must be made, including applying the Van-Vleck-Weisskopf correction and addressing the issue of "negative wavenumbers" in the line-by-line algorithm.

The recommended number of atmospheric levels is around 100. It is not advised to use atmospheres with more than 200 levels, as high-resolution atmospheres in terms of altitude will result in extended computational times, reducing the efficiency of the codes.

### Note 0:
The current resolution at which PT tables are calculated is fixed and determined by the Doppler half-width in the far-IR. The resolution is `deltaWV/NT = 10/20480 ≈ 5 * 10<sup>-4</sup>`. Such resolution is excessive for calculations in the near-IR. A dynamic resolution based on the spectral interval considered will soon be implemented, allowing it to be set as user input.

### Note 1:
Continuum absorption is not accounted for in this project. This functionality may be added later with contributions from new collaborators.

### Note 2:
Calculations are performed for one molecular species across all atmospheric levels in a single runtime.

### Note 3:
Field of View (FoV) functions are not implemented (see HAPI for reference).


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

## Command line parameters
### Parameters for `fpm run marfa` command
Required syntax: `fpm run marfa -- arg1 arg2 ... arg7 <arg8>`
| № | Argument      | Description                          | Required | Allowed values |
|-|---------------|--------------------------------------|----------|---------------|
|1| Molecule   | Species to calculate absorption features of. More molecules will be added   | Yes      | `CO2`, `H2O`           | 
|2| `Vstart`       | Left boundary of the spectral interval | Yes   | 10-10000 cm<sup>-1</sup> | 
|3| `Vend`          | Right boundary of the spectral interval | Yes      | 10-10000 cm<sup>-1</sup>| 
|4| cut-off condition | Distance from the center of the line from which absorption from this line is neglected | Yes | 10-500 cm<sup>-1</sup>  |
|5| χ-factor | χ-factor function name. Currently only CO2 corrections are employed. For more details see [χ-factors section](χ-factors)| Yes | `none`, `tonkov`, `pollack`, `perrin` |
|6| `targetValue` | Absorption feature to be written in the PT-table: volume absorption coefficient (km<sup>-1</sup>) or absorption cross-section (cm<sup>2</sup>/mol)| Yes       | `VAC` (volume absorption coefficient), `ACS` (absorption cross-section) |
|7| Atmosphere file name | Atmospheric file name, located in the `data/Atmospheres` directory. For the format of the file see: [Atmospheric Profile File Structure](#atmospheric-profile-file-structure)  | Yes | file names from the `data/Atmospheres` directory |
|8| `uuid` | ID for user request | No | Optional and is needed for the web-version of the program |

Examples:
```
fpm run marfa -- CO2 660 670 25 none ACS Venus1CO2.dat
fpm run marfa -- H2O 10 3000 250 none VAC VenusH2O.dat
```

### Parameters for `python scripts/postprocess.py` command
|Argument|Description|Required|Allowed values|
|--------|-----------|--------|--------------|
|subdirectory| Name of the PT-table directory| No, default value is where PT-tables from the latest run are stored. See the `output/ptTables/latest_run.txt` file | Any |
|`v1` | Start wavenumber from which you want to get processed data | Yes | `Vstart` < v1 < v2 < `Vend` |
|`v2` | End wavenumber to which you want to get processed data | Yes | `Vstart` < v1 < v2 < `Vend` |
|level| Atmospheric level at which you want to access data. Essentially means, that you access to the file `<level>.ptbin` | Yes | Normally from 1 to 100 (but see your atmospheric file)|
|resolution| Resolution at which you want to obtain the data. If you consider large intervals, it is not recommeded to use `high` resolution | Yes | `high` (4.8E-4cm<sup>-1</sup>), `medium` (4.8E-3cm<sup>-1</sup>), `coarse`(4.8E-2cm<sup>-1</sup>) |
| plot | Plot the data you postprocessed | No | Provide just a flag |

Examples:
```
python scripts/postprocess.py --v1 4032 --v2 4038 --level 40 --resolution high --plot
python scripts/posprocess.py --v1 10 --v2 3000 --level 50 --resolution coarse
python scripts/postprocess.py directory_name --v1 2500 --v2 --2550 --level 30 --resolution medium
```

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

## Output PT-table file
### General information
- PT-table files are generated and placed somwhere inside `output/ptTables` directory. 
- One PT-table file corresponds to one atmospheric level.
- Name of the PT-table file contains only level number, e.g. `1__.ptbin, 65_.ptbin`.
- PT-table file has an extention `.ptbin`.

### PT-file structure
The file consists of records with data, which could be directly accessed.
- Each record contains an array `RK` of high-resolution absorption data: either cross-section or absorption coefficient 
- Each record contains data about 10cm<sup>-1</sup> interval (`deltaWV` parameter in the code).
- Resolution of the data is defined by: `deltaWV/NT` = 10/20480 ≈ 4.8828×10<sup>-4</sup>cm<sup>-1</sup>.
- Relation between wavenumber of interest and record number: `record_number = int(WV/10)`. For example, if you want to know the value of absorption at 7560 cm<sup>-1</sup>, than you need to access the 756 record.
- Each record is of length `NT * 4` = 20481 * 4 bytes = 81924 bytes = 81,924 kilobytes

Schematic python code snippet for accessing data from this file:
```python
            import numpy as np

            # assuming you need to know absorption in 756 cm^-1
            V1 = 756

            # determining record number where this data stored
            recrod_number = int(V1 / 10.0)

            # number of values in one record
            NT = 20481

            # record length in bytes
            record_length = NT * 4

            with open('1__.ptbin', 'rb') as f:
              # start reading from record with record number "record_number"
              seek_position = (record_number-1) * record_length
              f.seek(seek_position)
  
              # read one record
              record_bytes = f.read(record_length)
  
              # Unpack an array of data using little-endian
              RK = np.frombuffer(record_bytes, dtype='<f4')
  
              # Converting data to float32 format (optional)
              RK = RK.astype(np.float32)
```
### Note

There is current minor vulnerabilty in the code regarding record numbers, which result in excess of PT-files size with increasing wavenumber. For example, you calculated PT-table for 6500-7000 cm<sup>-1</sup> spectral interval. First non-zero data will populate 650th record, leaving 1-649 records unpopulated. This excess of records may require additional and unnecessary storage. This issue would be resolved soon. But because of it, PT-table caculated for 100-110 cm<sup>-1</sup> will be much smaller in size than one calculated for 6500-6510 cm<sup>-1</sup>.
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
## Other spectral or molecular data
#### TIPS 
Total internal partition sums (TIPS) are needed for obtaining temperature-dependent spectral intensities. How the TIPS are implemented in MARFA:
- TIPS data are taken from Gamache work (Gamache 2017, see [References](#references))
- Data are available for first 74 isotopologues
- Covered temperature range is 20 - 1000 K
- For chosen isotope number TIPS as function of temperature can be accessed through the function: `TIPSofT` located in the `Spectroscopy.f90` module.

##### Outlook
- I plan to add recent TIPS (Gamache et al. 2021, see [References](#references))
- It might be better to organize input TIPS as en external subroutine based on Gamache's code Fortran or Python: TIPS_2021_v1p0.for, BD_TIPS_2021_v1p0.for or TIPS_2021_v1p0.py. I will soon return to it.
#### Molar masses
Molar masses are available in the `MolarMasses.f90` module. `WISO` array contains molar masses for 124 isotopolouges of first 42 molecules according to HITRAN numbering system.
## Performance estimations
Execution time at one atmospheric level largerly depends on number of spectral lines and line cut-off condition. Here are some benchmarks for Apple M1 chip:
|species|spectral interval (cm<sup>-1</sup>)|number of lines|cut off condition (cm<sup>-1</sup>)|execution time (s)|
|---------|-----------------|---------------|-----------------|----------------|
| CO<sub>2</sub>| 4000 - 4100 |  | 25 | 0.06 |
| CO<sub>2</sub>|   4000 - 4100|   | 250 | 0.24 |
| CO<sub>2</sub> | 10 - 3000 |  | 25 | 4.08 |
| CO<sub>2</sub> | 10 - 3000 |  | 250 | 25.2 |

Loop over atmospheric levels is currently not parallelized but I plan to do it with OpenMP in near time. So, presently time of processing atmospheric profile is linear to number of atmospheric levels.

Additional room for optimization might be organized with increasing effectiveness of Voigt function estimation algorithm and lastly negliecting weak lines.
## Introducing custom features
### Custom χ-factors
To add custom χ-factor function follow the steps:
1. Write a fortran function with χ-factor. Optionally it might be pressure or tempreature dependent. You can use `pressure` and `temperature` parameters inside a function.
2. Put this function at the end of the `chiFactors` module
3. Add new case to the `select-case` clause in the `fetchChiFactorFunction()` in the `ChiFactors` module.
4. Check that in the `LBL` module `chiCorrectedLorentz` line shape is used for the line wing description. There must be a line, like: `shapeFuncPtr => chiCorrectedLorentz`.

**Note:** Your χ-factor function must match the abstract interface for the χ-factor function: `chiFactorFuncPtr` (see `ShapeFuncInterface` module). Normally, it means that the function takes only one argument - distance from the line center in wavenumber in double precision and return only one value of `real` type. Check how the inputs and outputs are organized in predefined χ-factors functions and adjust accordingly.

#### Example: 
```fortran
    real function myChiFactor(X)
    real(kind=DP), intent(in) :: X
        
        myChiFactor = 1. ! default value of the chi-factor
        if (abs(X) > 4.) then
            if (abs(X) <= 125.) then
                myChiFactor = 1.24 * exp(-0.017*abs(X))
            else
                myChiFactor = 0.31 * exp(-0.036*abs(X))
            end if
        end if
        return
    end function myChiFactor
```
```fortran
    subroutine fetchChiFactorFunction()
        select case(trim(adjustl(chiFactorFuncName)))
        case ('none')
            chiFactorFuncPtr => noneChi
        case ('tonkov')
            chiFactorFuncPtr => tonkov
        case ('pollack')
            chiFactorFuncPtr => pollack
        case ('perrin')
            chiFactorFuncPtr => perrin
        ! ADD YOUR CHI-FACTOR HERE:
        case ('myChiFactor')
            chiFactorFuncPtr => myChiFactor
        end select
    end subroutine fetchChiFactorFunction
```
### Custom line shapes (for advanced users)
Currently in the `Shapes` module only standard line shapes are introduced: `lorentz`, `doppler`, `voigt`, `chiCorrectedLorentz` (Lorentz shape with wings corrected with χ-factor) and `correctedDoppler` (auxillary function, used only in Voigt shape calculation). To add your own line shape, you need to:
1. Provide a Fortran function of a custom line shape to the `Shapes` module. 
2. Go to the `LBL` module and assign `shapeFuncPtr` pointer to your line shape instead of the predefined shape.
3. Adjust the logic of the choice of line shape accordingly. This is required only for introducing "non-standard" (not Voigt, Lorentz or Doppler) line shapes.

**Note 1:** Essentialy, line shape function in `Shapes` module returns individual cross-section, because multiplication of line profile function to line intensity is included inside the function. Thus, multiplication on the temperature-dependent intensity must be provided. You can use predefined line intensity function `intensityOfT` in `Spectroscopy` module. Example:
```fortran
    real function myLineShape(X)
        real(kind=DP), intent(in) :: X

        ! ... YOUR LOGIC HERE

        myLineShape = myLineShape * intensityOfT(temperature) ! THIS LINE MUST BE PROVIDED
    end function myLineShape
```
If you want to use your own intensities, you must add new intensity function to the `Spectroscopy` module and refer to it.

**Note 2:** Your function must match the abstract interface for the line shape (see `ShapeFuncInterface` module). Normally, it means that the function takes only one argument - distance from the line center in wavenumber in double precision and return only one value of `real` type. Check how the inputs and outputs are organized in predefined line shape functions and adjust accordingly.

#### Example 1 
If you want to provide your own Voigt line shape, then it is the most straightforward. 
1. Add your function to the `Shapes` module.
2. Go to `src/LBL.f90` module and substitute code lines with `shapeFuncPtr => voigt` with e.g. `shapeFuncPtr => myVoigt`.

#### Example 2
If you want to introduce non-standard (not Voigt, Loretntz, Doppler) line shapes, then you additionally need to adjust the logic in `LBL` module yourself. In the current state, line shape function for calculation is a specific spectral point depends on how far this point from the center of the spectral line. If it is far enough, the `lorentz` or `chiCorrectedLorentz`, in the middle region `voigt` is applied and in the close vicinity of the line center `doppler` is used. If you want to introduce your own line shape you need to specify where it will be estimated, so adjust this pice of code in `LBL` module:
```fortran
if (shapePrevailFactor > BOUNDL) then
                ! Region utilizing Lorentz or chiCorrectedLorentz shape
                if (shiftedLineWV < startDeltaWV) then
                    shapeFuncPtr => Lorentz
                    call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr) 
                ! rest of the logic ...
else
                if (shapePrevailFactor > BOUNDD) then
                ! Region utilizing Voigt shape
                    if (shiftedLineWV < startDeltaWV) then
                        shapeFuncPtr => voigt
                        call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                    ! rest of the logic ... 
                else
                ! Region utilizing Doppler shape
                    if (shiftedLineWV < startDeltaWV ) then
                        shapeFuncPtr => doppler
                        call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                        ! rest of the logic ...
                end if
end if
```
## Troubleshooting
Feedback is awaited to populate this section. Most potential issues, such as invalid user inputs, in the Fortran source code and Python scripts are handled with clear error messages to facilitate troubleshooting.

If you refer to the χ-factor function which doesn't exist, then the `Segmentation Fault` error will be raised during runtime. This would be fixed later.

## How to cite this work

## License
This project is licensed under the MIT License. See the LICENSE file for more details.

## References
- Fomin, B. A. _Effective interpolation technique for line-by-line calculations of radiation absorption in gases._ Journal of Quantitative Spectroscopy and Radiative Transfer 53.6 (1995): 663-669.
- Tonkov, M. V., et al. _Measurements and empirical modeling of pure CO<sub>2</sub> absorption in the 2.3-μm region at room temperature: far wings, allowed and collision-induced bands._ Applied optics 35.24 (1996): 4863-4870.
- Pollack, James B., et al. _Near-infrared light from Venus' nightside: A spectroscopic analysis._ Icarus 103.1 (1993): 1-42.
- Perrin, M. Y., and J. M. Hartmann. _Temperature-dependent measurements and modeling of absorption by CO<sub>2</sub>-N<sub>2</sub> mixtures in the far line-wings of the 4.3 μm CO<sub>2</sub> band._ Journal of Quantitative Spectroscopy and Radiative Transfer 42.4 (1989): 311-317.
- Gamache, Robert R., et al. _Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres: Application to HITRAN2016 and beyond._ Journal of Quantitative Spectroscopy and Radiative Transfer 203 (2017): 70-87.
- Gordon, Iouli E., et al. _The HITRAN2020 molecular spectroscopic database._ Journal of quantitative spectroscopy and radiative transfer 277 (2022): 107949.
- Gamache, Robert R., et al. "Total internal partition sums for the HITRAN2020 database." Journal of Quantitative Spectroscopy and Radiative Transfer 271 (2021): 107713.
