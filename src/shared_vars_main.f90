module shared_vars_main
    use globals
    implicit none
    real(kind=DP) :: VSTART, VFINISH, WVA, WVB
    real(kind=DP) :: T, P, RO ! physical parameters in the atmosphere
    real(kind=DP) :: EPS ! line cut-off
end module shared_vars_main

module molecule_vars
    use globals
    implicit none
    integer :: MOTYPE
    integer :: NLIN
    character(len=17) :: LINE_PATH
end module molecule_vars