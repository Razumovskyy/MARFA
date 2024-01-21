module spectroscopy
    use kinds
    use constants, only: C2
    implicit none
    ! HITRAN spectral data variables
    ! Not 100% sure about the explanations !!
    real(kind=DP) :: V_I_ ! line position in cm^-1
    real :: S_I_  ! line intensity (typically for given temperature, often 296 K)
    
    ! Air-broadened half-width -- the half-width of the spectral line at half-maximum (HWHM)
    ! due to broadening by air (primarily nitrogen and oxygen)
    real :: ALFA_I_ ! HWHM due to air
    real :: ALFAS_I_ ! HWHM due to self-collisions
    
    real :: E_I_  ! Lower-state energy (energy of the lower state of the transition)
    
    ! A factor used to scale the line intensity, often related to isotopic abundance 
    ! or other factors that affect the overall intensity of the line.
    real :: FACT_I_ ! Intensity scaling factor
    
    integer :: NISO_I_ ! Isotopolouge number
    
    ! The shift in the line position due to pressure. 
    ! It represents the change in the central frequency of the line under different pressure conditions.
    real :: SHIFT_I_ ! Pressure-induced line shift

    ! ------- ^^^^ do not change the kinds there ^^^^ ------ !
contains
    real function VOIGT(X)
        real, intent(in) :: X
    end function VOIGT

    real function VV_LOR(X)
        real, intent(in) :: X
    end function VV_LOR

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
