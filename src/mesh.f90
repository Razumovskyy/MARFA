module MESH1
    implicit none
    
    ! Parameters
    integer, parameter :: NT0 = 10
    integer, parameter :: NT1 = NT0 * 2 !  = 20
    integer, parameter :: NT2 = NT1 * 2 !  = 40
    integer, parameter :: NT3 = NT2 * 2 !  = 80
    integer, parameter :: NT4 = NT3 * 2 !  = 160
    integer, parameter :: NT5 = NT4 * 2 !  = 320
    integer, parameter :: NT6 = NT5 * 2 !  = 640
    integer, parameter :: NT7 = NT6 * 2 !  = 1280
    integer, parameter :: NT8 = NT7 * 2 !  = 2560
    integer, parameter :: NT9 = NT8 * 2 !  = 5120
    integer, parameter :: NT = NT9 * 4 + 1 !  = 20481
    real, parameter :: OBR25 = 250.0 ! consider moving to the line_vars module
    real, parameter :: DELTA = 10.0
    real, parameter :: STEP = 1.0
    integer, parameter :: NINT = 10

    ! Arrays likely for various grids to cover the spectral shape
    ! consider setting double precision here !
    real :: RK(NT) ! array that likely holds calculated spectral data
    real :: RK0(NT0), RK0L(NT0), RK0P(NT0)
    real :: RK1(NT1), RK1L(NT1), RK1P(NT1)
    real :: RK2(NT2), RK2L(NT2), RK2P(NT2)
    real :: RK3(NT3), RK3L(NT3), RK3P(NT3)
    real :: RK4(NT4), RK4L(NT4), RK4P(NT4)
    real :: RK5(NT5), RK5L(NT5), RK5P(NT5)
    real :: RK6(NT6), RK6L(NT6), RK6P(NT6)
    real :: RK7(NT7), RK7L(NT7), RK7P(NT7)
    real :: RK8(NT8), RK8L(NT8), RK8P(NT8)
    real :: RK9(NT9), RK9L(NT9), RK9P(NT9)

    ! Additional variables
    real :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
    
end module MESH1