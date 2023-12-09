module shared_vars_main
    use globals
    implicit none
    real(kind=DP) :: VSTART, VFINISH, WVA, WVB
    real(kind=DP) :: T, P, RO ! physical parameters in the atmosphere
    real(kind=DP) :: EPS ! line cut-off
end module shared_vars_main