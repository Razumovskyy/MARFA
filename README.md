## Atmospheric Absorption Coefficient Calculator
### Overview
This repository hosts a Fortran-based tool designed for calculating atmospheric absorption coefficients in the longwave region. Leveraging HITRAN and HITEMP databases, it caters to the specific atmospheric profiles of Venus and Mars, though its flexible design permits adaptation to various celestial bodies.

### Features
- Atmospheric Profile Handling: Computes absorption coefficients based on atmospheric profiles, including height, pressure, temperature, and gas species density.
- Spectroscopic Line Shape Parameters: Supports various line shapes, including Lorentz, Doppler, Voigt, and Tonkov profiles, with the capability to integrate custom line shapes.
- Planetary Flexibility: Primarily focused on Venus and Mars, but adaptable for use in other planetary studies.
- van Vleck-Weisskopf-Huber Correction: Optional implementation for enhanced accuracy.
- PT-Tables Generation: Generates Pressure-Temperature tables usable in a broad range of atmospheric radiative transfer calculations.

## Usage
The program is designed to be straightforward for users familiar with Fortran and atmospheric sciences. Detailed instructions on compiling and running the application, along with the necessary input formats, are provided in the docs folder.

## License
This project is licensed under the MIT License. See the LICENSE file for more details.
