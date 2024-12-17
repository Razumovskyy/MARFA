module Spectroscopy
    ! This module contains calculation of spectroscopic features:
    ! pressure-induced shifted line positions, half-widths for various broadening mechanisms,
    ! establishing TIPS as a function of temperature and calculation of temprerature-dependent
    ! line intensities, simple line shapes, wing-correction functions and asymptotics for Voigt
    ! functions

    ! Line shape functions for self-consistent calculations are located in Shapes.f90 module
    
    use Constants
    use Atmosphere
    use ChiFactors
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
    integer :: moleculeIntCode ! 1 - H2O, 2 - CO2,  ...
    real, allocatable :: TIPS(:,:) ! TIPS array

    real :: dopHWHM, lorHWHM
    real :: lineIntensity

contains

    ! Add custom spectroscopic functions here in this module after the `contains` statement !

    pure function shiftedLinePosition(lineWVParameter, pressureParameter) result (position)
        ! calculation of shifted position of a transition wavenumber of a line shape
        ! see https://hitran.org/docs/definitions-and-units/ formula (7)
        
        implicit none
        
        real(kind=DP) :: position
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: pressureParameter

        position = lineWVParameter + deltaForeign*pressureParameter
    end function shiftedLinePosition

    
    pure function dopplerHWHM(lineWVParameter, temperatureParameter, molarMassParameter) result(dopplerWidth)
        ! calculation of a doppler half-width
        ! see https://hitran.org/docs/definitions-and-units/ formula (5)

        implicit none
        
        real :: dopplerWidth
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: temperatureParameter
        real, intent(in) :: molarMassParameter

        dopplerWidth = dopplerCONST * lineWVParameter * &
                sqrt(temperatureParameter / molarMassParameter)

    end function dopplerHWHM

    
    pure function lorentzHWHM(pressureParameter, partialPressureParameter, temperatureParameter) result(lorentzWidth)
        ! calculation of a Lorentz half-width: pressure- and temperature-dependent
        ! see https://hitran.org/docs/definitions-and-units/ formula (6)
        
        ! this function is called inside the LBL_LOOP where spectroscopic data: foreignTempCoeff,
        ! gammaForeign, gammaSelf are fixed, so that is why these parameters are not treated as
        ! as inputs for this function
        
        implicit none    
        
        real :: lorentzWidth
        real, intent(in) :: pressureParameter
        real, intent(in) :: partialPressureParameter
        real, intent(in) :: temperatureParameter

        lorentzWidth = ((refTemperature / temperatureParameter)**foreignTempCoeff) * &
                        (gammaForeign * (pressureParameter - partialPressureParameter) + &
                        gammaSelf * partialPressureParameter)

    end function lorentzHWHM

    
    real function TIPSofT(temperatureParameter)
        ! TODO: refactor when dealing with Gamache TIPS programs (python or Fortran)
        ! set isotope as a parameter for clarity to achieve pure function

        ! establishing a function which accepts temperature as an input
        ! and returning temperature-dependent TIPS as an output
        ! isotope -- is an external parameter from MolarMasses.f90 module

        ! in the TIPS array temperatures go from 20 K with 2 degrees step
        ! to find TIPS at random temperature, simple linear interpolation is employed
        
        implicit none
        
        real, intent(in) :: temperatureParameter
        integer :: lowerTempIndex ! Index of the closest temperature (available in TIPS array) below the input temperature
        real :: lowerTempValue ! Value of the closest temperature (available in TIPS array) below the input temperature
        real :: lowerWeight, upperWeight ! linear interpolation weights

        isotopeNum = jointMolIso / 100
        lowerTempIndex = (temperatureParameter - 20.0) / 2 + 1
        lowerTempValue = lowerTempIndex * 2.0 + 18.
        upperWeight = (temperatureParameter - lowerTempValue)/2.
        lowerWeight = 1. - upperWeight
        TIPSOfT = lowerWeight * TIPS(isotopeNum, lowerTempIndex) + upperWeight * TIPS(isotopeNum, lowerTempIndex+1)
    end function TIPSOfT

    
    real function intensityOfT(temperatureParameter)
        ! TODO: refactor when dealing with Gamache TIPS programs (python or Fortran)
        ! set isotope as a parameter for clarity
    
        ! calculation of intensity as a function of temperature
        ! see: https://hitran.org/docs/definitions-and-units/, formula (4)

        ! isotope -- is an external parameter from MolarMasses.f90 module
        
        implicit none
        
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


    ! real function parameterizedLorentzHWHM(pressureParameter, includeGammaSelf, partialPressureParameter, & 
    !                             includeTemperature, temperatureParameter)
    !     ! use this function for calculation Lorentz half-width if temperature is not known (reference temperature
    !     ! will be set) or self-broadening is not known (partial pressure will be set to zero)
        
    !     implicit none
        
    !     real, intent(in) :: pressureParameter
    !     logical, optional, intent(in) :: includeGammaSelf, includeTemperature
    !     real, optional, intent(in) :: partialPressureParameter
    !     real, optional, intent(in) :: temperatureParameter
        
    !     logical :: isIncludeGammaSelf, isIncludeTemperature
        
    !     ! defaults:
    !     isIncludeGammaSelf = .false.   ! do not count p_self, and gamma_self
    !     isIncludeTemperature = .false. ! no temperature dependency: temperature is set to 296 K
        
    !     if (present(includeTemperature)) isincludeTemperature = includeTemperature
    !     if (present(includeGammaSelf)) isIncludeGammaSelf = includeGammaSelf

    !     if (.not. isIncludeGammaSelf .and. .not. isIncludeTemperature) then
    !         ! temperature is set to 296 K and partial pressure is not counted
    !         parameterizedLorentzHWHM = gammaForeign * pressureParameter
    !     end if
        
    !     if (isIncludeGammaSelf .and. .not. isIncludeTemperature) then
    !         ! temperature is set to 296 K and partial pressure included
    !         parameterizedLorentzHWHM = gammaForeign * (pressureParameter - partialPressureParameter) + &
    !                         gammaSelf * partialPressureParameter
    !     end if

    !     if (.not. isIncludeGammaSelf .and. isIncludeTemperature) then
    !         ! temperature dependence is present, but partial pressure not included
    !         parameterizedLorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * (gammaForeign * pressureParameter)
    !     end if

    !     if (isIncludeGammaSelf .and. isIncludeTemperature) then
    !         ! full formula (6) from HITRAN docs 
    !         parameterizedLorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * &
    !                         (gammaForeign * (pressureParameter - partialPressureParameter) + &
    !                         gammaSelf * partialPressureParameter)
    !     end if
    ! end function parameterizedLorentzHWHM

    
    pure function lorentz(X, lorHWHMParameter) result(lorentzShape)
        implicit none

        real :: lorentzShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real, intent(in) :: lorHWHMParameter ! doppler HWHM

        lorentzShape = lorHWHMParameter / (pi*(X**2 + lorHWHMParameter**2)) * lineIntensity
    end function lorentz

    
    pure function doppler(X, dopHWHMParameter) result(dopplerShape)
        implicit none

        real :: dopplerShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real, intent(in) :: dopHWHMParameter ! doppler HWHM 

        dopplerShape = sqln2 / (sqrt(pi) * dopHWHMParameter) * exp(-(X/dopHWHMParameter)**2 * log(2.)) * lineIntensity
    end function doppler

    
    pure function chiCorrectedLorentz(X, lorHWHMParameter) result(chiCorrectedLorentzShape)
        ! Lorentz line shape with a χ-corrected wing
        implicit none

        real :: chiCorrectedLorentzShape
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        real, intent(in) :: lorHWHMParameter ! Lorentz HWHM

        chiCorrectedLorentzShape = lorentz(X, lorHWHMParameter) * chiFactorFuncPtr(X, moleculeIntCode)

    end function chiCorrectedLorentz

    
    pure function voigtAsymptotic1(X, lorHWHMParameter, VX) result(voigtAsymptotic1Value)
        ! Lorentz leading, Doppler-influenced asymptotic correction for the Voigt function (y<<1, x>>1)
        implicit none

        real :: voigtAsymptotic1Value
        real(kind=DP), intent(in) :: X
        real, intent(in) :: lorHWHMParameter
        real, intent(in) :: VX ! Voigt function K(x,y): x parameter
        
        voigtAsymptotic1Value = chiCorrectedLorentz(X, lorHWHMParameter) * (1 + 1.5/VX**2)
    end function voigtAsymptotic1


    pure function voigtAsymptotic2(X, lorHWHMParameter, VX, dopHWHMParameter) result(voigtAsymptotic2Value)
        ! Voigt ≈ Doppler + Lorentz + correction (y<<1, x>1)
        implicit none 

        real :: voigtAsymptotic2Value
        real(kind=DP), intent(in) :: X
        real, intent(in) :: lorHWHMParameter, dopHWHMParameter
        real, intent(in) :: VX ! Voigt function K(x,y): x parameter
        real, parameter :: U(9) = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
        real, parameter :: W(9) = [-0.688, 0.2667, 0.6338, 0.4405, 0.2529, 0.1601, 0.1131, 0.0853, 0.068]
        real :: F
        integer :: I

        I = VX/0.5 - 1.00001
        F = 2. * (W(I)*(U(I+1)-VX) + W(I+1)*(VX-U(I)))
        voigtAsymptotic2Value = doppler(X, dopHWHMParameter) + (lorHWHMParameter/(pi*X**2) * (1.+F)) * lineIntensity
    end function voigtAsymptotic2

end module Spectroscopy
