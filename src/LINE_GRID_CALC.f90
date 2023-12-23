! we need to pass function (line shape) itself to the subroutines arguments (LEFTLBL, CENTLBL, RIGHTLBL)
! in order to do that consider creating abstract interface that defines the procedure's signature (arguments and return type)
! sample code below: 

! module shape_functions
!     implicit none
!     abstract interface
!         real function shape_func(x)
!             real, intent(in) :: x
!         end function shape_func
!     end interface
! end module shape_functions


module LINE_GRID_CALC
    use globals
    use MESH1
    use spectroscopy, only: VOIGT, VV_LOR, DOPLER
    implicit none
contains
    subroutine LEFTLBL(FREQ, UL, FSHAPE, EPS)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real :: EPS
        real :: FSHAPE ! DOPLER, VOIGT or Lorentz !!! remove later
        ! procedure(shape_func), pointer :: ShapeFunction
        ! Use ShapeFunction as a function pointer
    end subroutine LEFTLBL
    
    subroutine CENTLBL(FREQ, UL, FSHAPE, EPS)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real :: EPS
        real :: FSHAPE ! DOPLER, VOIGT or Lorentz !!! remove later
        ! procedure(shape_func), pointer :: ShapeFunction
        ! Use ShapeFunction as a function pointer
    end subroutine CENTLBL

    subroutine RIGHTLBL(FREQ, UL, FSHAPE, EPS)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real :: EPS
        real :: FSHAPE ! DOPLER, VOIGT or Lorentz !!! remove later
        ! procedure(shape_func), pointer :: ShapeFunction
        ! Use ShapeFunction as a function pointer
    end subroutine RIGHTLBL

end module LINE_GRID_CALC