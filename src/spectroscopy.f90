module Spectroscopy
    use Constants
    use Atmosphere
    use MolarMasses
    implicit none
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
    real, allocatable :: TIPS(:,:) ! TIPS array (Total internal partition sums)

contains

    ! molarMass = WISO(isotopeNum)
    
    real(kind=DP) function shiftedLinePosition(lineWVParameter, pressureParameter)
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: pressureParameter

        shiftedLinePosition = lineWVParameter + deltaForeign * pressureParameter
    end function shiftedLinePosition

    real function dopplerHWHM(lineWVParameter, temperatureParameter, molarMassParameter)
        real(kind=DP), intent(in) :: lineWVParameter
        real, intent(in) :: temperatureParameter
        real, intent(in) :: molarMassParameter

        dopplerHWHM = dopplerCONST * shiftedLinePosition(lineWVParameter, pressure) * &
                sqrt(temperatureParameter / molarMassParameter)

    end function dopplerHWHM

    real function lorentzHWHM(pressureParameter, includeGammaSelf, partialPressureParameter, & 
                                includeTemperature, temperatureParameter)
        ! TODO: add check if includeTemperature=true but not passed as an argument
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

        ! TODO: rewrite using select case: use integer indicator which is mapped to two bool variables:
        ! if (isIncludeGammaSelf) then
        !     caseIndex = caseIndex + 2
        ! endif
        ! if (isIncludeTemperature) then
        !     caseIndex = caseIndex + 1
        ! endif

        if (.not. isIncludeGammaSelf .and. .not. isIncludeTemperature) then
            ! temperature is set to 296 K and partial pressure is not counted
            lorentzHWHM = gammaForeign * pressureParameter
        end if
        
        if (isIncludeGammaSelf .and. .not. isIncludeTemperature) then
            ! temperature is set to 296 K and partial pressure included
            lorentzHWHM = gammaForeign * (pressureParameter - partialPressureParameter) + &
                            gammaSelf * partialPressureParameter
        end if

        if (.not. isIncludeGammaSelf .and. isIncludeTemperature) then
            ! temperature dependence is present, but partial pressure not included
            lorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * (gammaForeign * pressureParameter)
        end if

        if (isIncludeGammaSelf .and. isIncludeTemperature) then
            ! full formula (6) from HITRAN docs 
            lorentzHWHM = ((refTemperature / temperatureParameter)**foreignTempCoeff) * &
                            (gammaForeign * (pressureParameter - partialPressureParameter) + &
                            gammaSelf * partialPressureParameter)
        end if
    end function lorentzHWHM

    real function intensityOfT(temperatureParameter)
        real, intent(in) :: temperatureParameter ! [K] -- temperature at the current atmospheric level
        ! real(kind=DP), intent(in) :: nu ! [cm-1],  gridWV -- spectral point where absorption coefficitent will be calculated
        ! procedure(shape), pointer, intent(in) :: lineShape ! [cm] -- normalized line shape function from the MyShapes module
        ! --------------------------------------------------- !
        real(kind=DP) :: shiftedLineWV
        real :: intensity ! [cm-1/(molecule*cm-2)] (refLineIntensity) -- the spectral line intensity for a single molecule per unit volume.
        real :: TIPSFactor, boltzmannFactor, emissionFactor
        real :: TIPSOfT
        real :: TIPSOfRefT
        integer :: NTAB_G
        real :: C_G1, C_G2
        real :: t_G1
        integer :: isotopeNum

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)

        isotopeNum = jointMolIso / 100
        NTAB_G = (temperatureParameter - 20.0) / 2 + 1
        t_G1 = NTAB_G * 2.0 + 18.
        C_G2 = (temperatureParameter - t_G1)/2.
        C_G1 = 1. - C_G2
        TIPSOfT = C_G1 * TIPS(isotopeNum, NTAB_G) + C_G2 * TIPS(isotopeNum, NTAB_G+1)
        TIPSOfRefT = TIPS(isotopeNum, 139)

        TIPSFactor = TIPSOfRefT / TIPSOfT
        boltzmannFactor = exp(-C2*lineLowerState/temperatureParameter) / exp(-C2*lineLowerState/refTemperature)
        emissionFactor = (1 - exp(-C2*lineWV/temperatureParameter)) / (1 - exp(-C2*lineWV/refTemperature))

        intensity = refLineIntensity * TIPSFactor * boltzmannFactor * emissionFactor
    end function intensityOfT

end module Spectroscopy
