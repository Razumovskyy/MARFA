module Constants
    implicit none
    integer, parameter :: DP = selected_real_kind(15,  307) ! selected kind for double precision
    
    real, parameter :: PI = 3.1415926
    real, parameter :: sqln2 = sqrt(log(2.))

    ! CGS system is used
    real, parameter :: PLANCK = 6.626070e-27 ! [erg*s]
    real, parameter :: SPL = 2.99792458e10 ! [cm/s] -- speed of light
    real, parameter :: BOL = 1.3806503e-16 ! [erg/K] -- Boltzmann constant
    real, parameter :: C2 = (PLANCK * SPL / BOL) ! 1.438769 [cm * K] -- second radiation constant
    real, parameter :: stPressure = 1. ! [atm]
    real, parameter :: stTemperature = 273.15 ! [K]
    real, parameter :: LOSCHMIDT = 2.6867811e25 ! [1/m^3]
    real, parameter :: AVOGADRO = 6.02214076e23 ! [1/mol]
    real, parameter :: gasCONST = AVOGADRO * BOL ! [erg/(mol*K)]
    real, parameter :: dopplerCONST = sqrt(2*AVOGADRO*BOL*log(2.)) / SPL

    real, parameter :: refTemperature = 296. ! [K] -- refrence temperature for HITRAN data 

    ! ------------------------------------------ !
    real, parameter :: C1 = 2 * PI * PLANCK * (SPL**2) ! SOMETIMES PI is omitted
    real, parameter :: redPLANCK = PLANCK / (2 * PI) ! reduced planck constant

end module Constants
