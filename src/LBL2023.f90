module line_shapes
    use globals
    implicit none
contains
    real function VOIGT(XXX)
        real(kind=DP), intent(in) :: XXX
    end function VOIGT

    real(kind=DP) function VV_LOR(X)
        real(kind=DP), intent(in) :: X
    end function VV_LOR

    real(kind=DP) function DOPLER(X)
        real(kind=DP), intent(in) :: X
    end function DOPLER

    real(kind=DP) function VAN_VLE(T, V)
        real(kind=DP) :: T, V
    end function VAN_VLE
end module line_shapes

module LBL
    use globals
    use shared_vars_main
    use MESH1
    use molecule_vars
    use line_shapes
    implicit none
contains
    
    subroutine LBL2023(MO_E, LINBEG, VAA, VFISH, NLIN)
        ! consider avoid providing arguments in subroutine.
        ! if they are declared in the external module -- why to declare them twice ?

        integer :: MO_E ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
        integer :: LINBEG
        integer :: NLIN ! number of lines in one file
        ! consider removing it from the subroutine argument list, otherwise warning 
        ! that this variable is already defined in the parent scope
        real(kind=DP) :: VAA, VFISH

    end subroutine
end module LBL
