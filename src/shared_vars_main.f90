module shared_vars_main
    use globals
    implicit none
    real(kind=DP) :: VSTART, VFINISH, WVA, WVB
    real(kind=DP) :: T, P, RO ! physical parameters in the atmosphere
    real(kind=DP) :: EPS ! line cut-off
end module shared_vars_main

module molecule_vars
    ! variables related to the molecules
    use globals
    implicit none
    integer :: MOTYPE ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
    integer :: NLIN ! number of lines in one HITRAN file (e.g. H16.01)
    character(len=20) :: LINE_PATH ! file names with HITRAN raw data (one file per molecule)
end module molecule_vars

module van_fleck_huber_vars
    implicit none
    integer :: JM1 ! likely just a loop index
    real :: VIVI ! modified wavenumber within the loop
    real :: EVV, EVV_ ! variables related to the exponential terms in the Van Vleck-Weisskopf-Huber factor calculations
    real :: FACTV ! actual Van-Vleck factor applied
end module van_fleck_huber_vars
