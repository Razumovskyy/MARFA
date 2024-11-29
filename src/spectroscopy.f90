module Spectroscopy
    ! This module contains calculation of spectroscopic features:
    ! pressure-induced shifted line positions, half-widths for various broadening mechanisms,
    ! establishing TIPS as a function of temperature and calculation of temprerature-dependent
    ! line intensities

    ! Add custom spectroscopic functions here in this module after the `contains` statement !
    
    use Constants
    use Atmosphere
    use MolarMasses
    implicit none
    
    ! Line spectral parameters
    ! For more details about naming and HITRAN - see app/processParFile.f90
    real(kind=DP) :: lineWV ! [cm-1] -- current spectral line, transition wavenumber
    real :: refLineIntensity ! [cm-1/(molecule*cm-2)] -- spectral line intensity at refTemperature=296 K
    real :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM at refTemperature=296 K
    real :: gammaSelf ! [cm-1/atm] -- self-broadened component of Lorentz HWHM
    real :: lineLowerState ! [cm-1] -- lower state energy of the transition
    real :: foreignTempCoeff ! [dimensionless] (coefficient for temperature dependence of gammaForeign)
    integer :: jointMolIso ! [dimensionless] custom variable: joined reference to Molecule number (MOL) and Isotopologue number (ISO)
    real :: deltaForeign ! [cm-1/atm] (pressure shift of the line position at 296 K and 1 atm)
    
    real :: molarMass ! [g/mol] -- current species molar mass
    integer :: molType ! to define type of the molecule : 2 -CO2, 1 -H2O, 0 -other
    real, allocatable :: TIPS(:,:) ! TIPS array

contains

    real(kind=DP) function shiftedLinePosition(lineWVParameter, pressureParameter)
        ! calculation of shifted position of a transition wavenumber
        ! see https://hitran.org/docs/definitions-and-units/ formula (7)
        
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: pressureParameter

        shiftedLinePosition = lineWVParameter + deltaForeign*pressureParameter
    end function shiftedLinePosition

    real function dopplerHWHM(lineWVParameter, temperatureParameter, molarMassParameter)

        ! calculation of a doppler half-width
        ! see https://hitran.org/docs/definitions-and-units/ formula (5)

        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: temperatureParameter
        real, intent(in) :: molarMassParameter

        ! dopplerHWHM = dopplerCONST * abs(shiftedLinePosition(lineWVParameter, pressure)) * &
        !         sqrt(temperatureParameter / molarMassParameter)

        ! more correct way 
        dopplerHWHM = dopplerCONST * lineWVParameter * &
                sqrt(temperatureParameter / molarMassParameter)

    end function dopplerHWHM

    real function lorentzHWHM(pressureParameter, partialPressureParameter, temperatureParameter)
        ! calculation of a Lorentz half-width: pressure- and temperature-dependent
        ! see https://hitran.org/docs/definitions-and-units/ formula (6)
        
        ! this function is called inside the LBL_LOOP where spectroscopic data: foreignTempCoeff,
        ! gammaForeign, gammaSelf are fixed, so that is why these parameters are not treated as
        ! as inputs for this function
        
        real, intent(in) :: pressureParameter
        real, intent(in) :: partialPressureParameter
        real, intent(in) :: temperatureParameter

        lorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * &
                        (gammaForeign * (pressureParameter - partialPressureParameter) + &
                        gammaSelf * partialPressureParameter)

    end function lorentzHWHM

    real function TIPSofT(temperatureParameter)
        ! establishing a function which accepts temperature as an input
        ! and returning temperature-dependent TIPS as an output
        real, intent(in) :: temperatureParameter
        
        ! TODO: refactor when dealing with Gamache TIPS programs (python or Fortran)
        integer :: NTAB_G
        real :: C_G1, C_G2
        real :: t_G1
        
        isotopeNum = jointMolIso / 100
        NTAB_G = (temperatureParameter - 20.0) / 2 + 1
        t_G1 = NTAB_G * 2.0 + 18.
        C_G2 = (temperatureParameter - t_G1)/2.
        C_G1 = 1. - C_G2
        TIPSOfT = C_G1 * TIPS(isotopeNum, NTAB_G) + C_G2 * TIPS(isotopeNum, NTAB_G+1)
    end function TIPSOfT

    real function intensityOfT(temperatureParameter)
        ! calculation of intensity as a function of temperature
        ! see: https://hitran.org/docs/definitions-and-units/, formula (4)
        
        real, intent(in) :: temperatureParameter ! [K] -- temperature at the current atmospheric level

        real(kind=DP) :: shiftedLineWV
        real :: TIPSFactor, boltzmannFactor, emissionFactor
        real :: TIPSOfRefT

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)

        isotopeNum = jointMolIso / 100
        TIPSOfRefT = TIPS(isotopeNum, 139)

        TIPSFactor = TIPSOfRefT / TIPSOfT(temperatureParameter)
        boltzmannFactor = exp(-C2*lineLowerState/temperatureParameter) / exp(-C2*lineLowerState/refTemperature)
        emissionFactor = (1 - exp(-C2*lineWV/temperatureParameter)) / (1 - exp(-C2*lineWV/refTemperature))

        intensityOfT = refLineIntensity * TIPSFactor * boltzmannFactor * emissionFactor
    end function intensityOfT

    ! ----------------------------------------------------------------------------------------------------------------!

    real function parameterizedLorentzHWHM(pressureParameter, includeGammaSelf, partialPressureParameter, & 
                                includeTemperature, temperatureParameter)
        ! use this function for calculation of Lorentz half-width if temperature is not known (reference temperature
        ! will be set) or self-broadening is not known
        
        ! TODO:(?) add check if includeTemperature=true but not passed as an argument
        real, intent(in) :: pressureParameter
        logical, optional, intent(in) :: includeGammaSelf, includeTemperature
        real, optional, intent(in) :: partialPressureParameter
        real, optional, intent(in) :: temperatureParameter
        
        logical :: isIncludeGammaSelf, isIncludeTemperature
        
        ! defaults:
        isIncludeGammaSelf = .false.   ! do not count p_self, and gamma_self
        isIncludeTemperature = .false. ! no temperature dependency: temperature is set to 296 K
        
        if (present(includeTemperature)) isincludeTemperature = includeTemperature
        if (present(includeGammaSelf)) isIncludeGammaSelf = includeGammaSelf

        if (.not. isIncludeGammaSelf .and. .not. isIncludeTemperature) then
            ! temperature is set to 296 K and partial pressure is not counted
            parameterizedLorentzHWHM = gammaForeign * pressureParameter
        end if
        
        if (isIncludeGammaSelf .and. .not. isIncludeTemperature) then
            ! temperature is set to 296 K and partial pressure included
            parameterizedLorentzHWHM = gammaForeign * (pressureParameter - partialPressureParameter) + &
                            gammaSelf * partialPressureParameter
        end if

        if (.not. isIncludeGammaSelf .and. isIncludeTemperature) then
            ! temperature dependence is present, but partial pressure not included
            parameterizedLorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * (gammaForeign * pressureParameter)
        end if

        if (isIncludeGammaSelf .and. isIncludeTemperature) then
            ! full formula (6) from HITRAN docs 
            parameterizedLorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * &
                            (gammaForeign * (pressureParameter - partialPressureParameter) + &
                            gammaSelf * partialPressureParameter)
        end if
    end function parameterizedLorentzHWHM

end module Spectroscopy
