module MyConstants
    implicit none
    integer, parameter :: DP = selected_real_kind(15,  307)
    real, parameter :: pi = 3.14159

    real, parameter :: SPL = 2.99792458e10
    real, parameter :: BOL = 1.3806503e-16
    real, parameter :: AVOGADRO = 6.02214076e23
    real, parameter :: dopplerCONST = sqrt(2*AVOGADRO*BOL*log(2.)) / SPL

    real, parameter :: refTemperature = 296. ! in K
end module MyConstants