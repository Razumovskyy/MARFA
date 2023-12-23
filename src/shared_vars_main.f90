module shared_vars_main
    use globals
    implicit none
    real(kind=DP) :: VSTART, VFINISH, WVA, WVB
    real :: T, P, RO ! physical parameters in the atmosphere
    real(kind=DP) :: EPS ! line cut-off
    real, allocatable :: QofT(:,:)
end module shared_vars_main
