module shared_vars_main
    use kinds
    implicit none
    real(kind=DP) :: startWV, endWV
    real(kind=DP) :: extStartWV, extEndWV
    real(kind=DP) :: startDeltaWV, endDeltaWV ! starting and ending wavenumbers inside deltaWV interval *** ! VSTART, VFINISH !
    real(kind=DP) :: extStartDeltaWV, extEndDeltaWV ! extended deltaWV intervals including cut-offs ! WVA, WVB !
    real :: T, P, RO ! physical parameters of the atmosphere on the current height level
    real(kind=DP) :: EPS ! ??? epsilon ???
    real, allocatable :: QofT(:,:)
end module shared_vars_main
