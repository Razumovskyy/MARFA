module spectroscopy
    use kinds
    use constants, only: C2
    implicit none
    ! HITRAN spectral data variables
    ! Not 100% sure about the explanations !!
    real(kind=DP) :: lineWV ! The wavenumber of the spectral line transition (cm-1) in vacuum
    real :: lineIntensity  ! intensity in cm-1/(molec * cm-2) at 296 Kelvin
    
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

    ! ------- ^^^^ do not change the kinds there ^^^^ ------ !
contains
    real function VOIGT(X)
        real, intent(in) :: X
    end function VOIGT

    real function LORENTZ(X)
        real, intent(in) :: X
    end function LORENTZ

    real function DOPLER(X)
        real, intent(in) :: X
    end function DOPLER

    real function VAN_VLE(T, V)
        real :: T, V
    end function VAN_VLE

    real function R_factor(t, nu)
        real :: T
        real(kind=DP) :: nu
       
        R_factor = nu * tanh(C2*nu/(2*T))
    end function R_factor
end module spectroscopy

module van_fleck_huber_vars
    use kinds
    implicit none
    integer :: JM1 ! likely just a loop index
    real(kind=DP) :: lineShapeWV ! ** ! VIVI ! likely the grid line shape line position for grid calculations
    real :: numeratorCloughFactor ! numerator of the Clough Factor see works by Clough 1989.
    real :: cloughFactorExpTerm ! EVV_ ! ** ! auxilary variable for exponential factor in the R-fucntion
end module van_fleck_huber_vars
