module constants
    implicit none
    real, parameter :: PI = 3.1415926
    real, parameter :: PLANCK = 6.626070e-34
    real, parameter :: redPLANCK = PLANCK / (2 * PI)
    real, parameter :: SPL = 2.99792458e8
    real, parameter :: BOL = 1.3806503e-23
    real, parameter :: C1 = 2 * PI * PLANCK * (SPL**2)
    real, parameter :: C2 = (PLANCK * SPL / BOL) * 100 ! 1.438769 cm * K
end module constants
