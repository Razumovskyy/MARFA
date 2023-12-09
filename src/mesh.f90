module MESH1
    implicit none
    
    ! Parameters
    integer, parameter :: NT0 = 10
    integer, parameter :: NT1 = NT0 * 2
    integer, parameter :: NT2 = NT1 * 2
    integer, parameter :: NT3 = NT2 * 2
    integer, parameter :: NT4 = NT3 * 2
    integer, parameter :: NT5 = NT4 * 2
    integer, parameter :: NT6 = NT5 * 2
    integer, parameter :: NT7 = NT6 * 2
    integer, parameter :: NT8 = NT7 * 2
    integer, parameter :: NT9 = NT8 * 2
    integer, parameter :: NT = NT9 * 4 + 1
    real, parameter :: OBR25 = 250.0
    real, parameter :: DELTA = 10.0
    real, parameter :: STEP = 1.0
    integer, parameter :: NINT = 10

    ! Arrays
    real :: RK(NT)
    real :: RK0(NT0)
    real :: RK0L(NT0)
    real :: RK0P(NT0)
    real :: RK1(NT1)
    real :: RK1L(NT1)
    real :: RK1P(NT1)
    real :: RK2(NT2)
    real :: RK2L(NT2)
    real :: RK2P(NT2)
    real :: RK3(NT3)
    real :: RK3L(NT3)
    real :: RK3P(NT3)
    real :: RK4(NT4)
    real :: RK4L(NT4)
    real :: RK4P(NT4)
    real :: RK5(NT5)
    real :: RK5L(NT5)
    real :: RK5P(NT5)
    real :: RK6(NT6)
    real :: RK6L(NT6)
    real :: RK6P(NT6)
    real :: RK7(NT7)
    real :: RK7L(NT7)
    real :: RK7P(NT7)
    real :: RK8(NT8)
    real :: RK8L(NT8)
    real :: RK8P(NT8)
    real :: RK9(NT9)
    real :: RK9L(NT9)
    real :: RK9P(NT9)

    ! Additional variables
    real :: H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H
    
end module MESH1