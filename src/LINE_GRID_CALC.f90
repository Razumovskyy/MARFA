module LINE_GRID_CALC
    use globals
    use MESH1
    use shape_functions
    implicit none
contains
    subroutine LEFTLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape_func), pointer :: FSHAPE
    end subroutine LEFTLBL
    
    subroutine CENTLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape_func), pointer :: FSHAPE
    end subroutine CENTLBL

    subroutine RIGHTLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape_func), pointer :: FSHAPE
    end subroutine RIGHTLBL

end module LINE_GRID_CALC
