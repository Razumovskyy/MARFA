module Spectroscopy
    use Kinds
    use Constants, only: C2, dopplerCONST, refTemperature
    implicit none
    ! HITRAN spectral data variables
    ! Not 100% sure about the explanations !!
    real(kind=DP) :: lineWV ! The wavenumber of the spectral line transition (cm-1) in vacuum
    real :: refLineIntensity  ! intensity in cm-1/(molec * cm-2) at 296 Kelvin
    
    ! Air-broadened half-width -- the half-width of the spectral line at half-maximum (HWHM)
    ! due to broadening by air (primarily nitrogen and oxygen)
    real :: gammaForeign ! ALFA_I ! The air-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm
    real :: gammaSelf ! ALFAS_I ! The self-broadened half width at half maximum (HWHM) (cm−1/atm) at Tref=296K and reference pressure pref=1atm
    
    real :: lineLowerState  ! E_I_ ! The lower-state energy of the transition (cm-1)
    
    ! A factor used to scale the line intensity, often related to isotopic abundance 
    ! or other factors that affect the overall intensity of the line.
    real :: foreignTempCoeff ! FACT_I_ ! The coefficient of the temperature dependence of the air-broadened half width
    
    integer :: jointMolIso! NISO_I_ ! joined reference to Molecule number (MOL) and Isotopologue number (ISO)
    
    ! The shift in the line position due to pressure. 
    ! It represents the change in the central frequency of the line under different pressure conditions.
    real :: deltaForeign ! SHIFT_I_ ! The pressure shift (cm−1/atm) at Tref=296K and pref=1atm of the line position with respect to the vacuum transition wavenumber νij

    ! ------- ^^^^ do not change the Kinds there ^^^^ ------ !
contains

    real function LorentzianHWHM(p, t, pSelf, tRef, n, gammaS, gammaF)
    ! see HITRAN docs, section 'Temperature and pressure dependence of the line width', formula 6
        real :: p, pSelf
        real :: t, tRef
        real :: n
        real :: gammaS, gammaF

        LorentzianHWHM = ((tRef/t) ** n) * (gammaF * (p - pSelf) + gammaS * pSelf)
    end function LorentzianHWHM 

    real function pressureShiftedWV(nu, deltaF, p)
        real(kind=DP) :: nu
        real :: deltaF
        real :: p
        
        pressureShiftedWV = nu + deltaF * p
    end function pressureShiftedWV

    real function DopplerHWHM(nu, t, molarMass)
        real(kind=DP) :: nu
        real :: t
        real :: molarMass

        DopplerHWHM = dopplerCONST * nu * sqrt(t/molarMass) 
    end function

    real function lineIntensityofT(t, refLineIntensity, isotopeNum, partFunc, lowerState, transitionWV)
        real, intent(in) :: t
        real, intent(in) :: refLineIntensity
        integer, intent(in) :: isotopeNum
        real, intent(in) :: partFunc(:,:)
        real, intent(in) :: lowerState 
        real(kind=DP), intent(in) :: transitionWV

        integer :: NTAB_G
        real :: C_G1, C_G2
        real :: t_G1
        real :: partFuncOfT, partFuncOfRefT
        real :: partFuncFactor
        real :: bolFactor
        real :: stEmissFactor

        NTAB_G = (t - 20.0)/2 + 1
        t_G1 = NTAB_G * 2.0 + 18.
        C_G2 = (t - t_G1)/2.
        C_G1 = 1. - C_G2
        partFuncOfT = C_G1 * partFunc(isotopeNum, NTAB_G) + C_G2 * partFunc(isotopeNum, NTAB_G+1)
        partFuncOfRefT = partFunc(isotopeNum, 139)
        partFuncFactor = partFuncOfRefT / partFuncOfT
        bolFactor = exp(-C2*lowerState/t) / exp(-C2*lowerState/refTemperature) 
        stEmissFactor = (1 - exp(-C2*transitionWV/t)) / (1 - exp(-C2*transitionWV/refTemperature))
        
        ! TODO consider to use stEmissFactor in Longwave, found it was 0.4 for 50 - 80 cm-1 for H2O
        ! lineIntensityofT = refLineIntensity * partFuncFactor * bolFactor * stEmissFactor
        lineIntensityofT = refLineIntensity * partFuncFactor * bolFactor
    end function lineIntensityofT
end module Spectroscopy
