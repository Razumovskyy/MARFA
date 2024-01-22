module constants
    implicit none
    real, parameter :: PI = 3.1415926
    real, parameter :: PLANCK = 6.626070e-34
    real, parameter :: redPLANCK = PLANCK / (2 * PI)
    real, parameter :: SPL = 2.99792458e8
    real, parameter :: BOL = 1.3806503e-23
    real, parameter :: C1 = 2 * PI * PLANCK * (SPL**2)
    real, parameter :: C2 = (PLANCK * SPL / BOL) * 100 ! 1.438769 cm * K
    real, parameter :: stPressure = 1. ! in atmospheres
    real, parameter :: stTemperature = 273.15 ! in K
    real, parameter :: LOSCHMIDT = 2.6867811e25
    real, parameter :: AVOGADRO = 6.02214076e23
    real, parameter :: gasCONST = AVOGADRO * BOL
    real, parameter :: dopplerCONST = sqrt(2*AVOGADRO*BOL*log(2.)) / SPL

    real, parameter :: refTemperature = 296. ! in K

end module constants
