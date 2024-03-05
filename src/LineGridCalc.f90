module LineGridCalc
    use Kinds
    use Mesh
    use ShapeFuncInterface
    implicit none
contains
    subroutine leftLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape), pointer :: FSHAPE
    end subroutine leftLBL
    
    subroutine centerLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape), pointer :: FSHAPE
    end subroutine centerLBL

    subroutine rightLBL(FREQ, UL, FSHAPE)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL
        real(kind=DP) :: EPS
        procedure(shape), pointer :: FSHAPE
    end subroutine rightLBL

end module LineGridCalc
