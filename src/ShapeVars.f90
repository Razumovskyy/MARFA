module ShapeVars
    use Kinds
    implicit none
    integer :: molType ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
    real(kind=DP) :: VI
    real :: SLL
    real(kind=DP) :: gammaPT ! AL ! Lorentz HWHM -- T and P dependent !
end module ShapeVars
