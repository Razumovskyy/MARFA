module Mesh
    implicit none
    
    ! Parameters: number of grid points for each grid:
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
    
    ! TODO:(!!): move it out from parameters, it might be an input value, or
    ! a dictionary value based on the input
    real, parameter :: deltaWV = 10.0 ! resolution in cm-1
    
    ! TODO:(!!) deal with it during refactor of the grid calculation
    real, parameter :: STEP = 1.0

    real :: cutOff ! cutOff condition in cm-1

    ! Arrays for various grids to cover the spectral shape !

    ! Naming example explanation:
    ! RK2L -- calculated contribution from left parts of lines on 2nd grid
    ! RK2R -- calculated contribution from right parts of lines on 2nd grid
    ! RK2 -- calculated contribution from central parts of lines on 2nd grid
    ! These values are used for sequenced interpolation to finally get the RK values.
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

    real :: RK(NT) ! array that holds final calculated spectral data

    ! Grids relative resolutions
    real, parameter :: H0 = STEP
    real, parameter :: H1 = H0 / 2.0
    real, parameter :: H2 = H1 / 2.0
    real, parameter :: H3 = H2 / 2.0
    real, parameter :: H4 = H3 / 2.0
    real, parameter :: H5 = H4 / 2.0
    real, parameter :: H6 = H5 / 2.0
    real, parameter :: H7 = H6 / 2.0
    real, parameter :: H8 = H7 / 2.0
    real, parameter :: H9 = H8 / 2.0
    real, parameter :: H = H9 / 4.0

contains

    subroutine resetAbsorptionGridValues
        implicit none
        RK = 0.0
        RK0 = 0.0; RK0L = 0.0; RK0P = 0.0
        RK1 = 0.0; RK1L = 0.0; RK1P = 0.0
        RK2 = 0.0; RK2L = 0.0; RK2P = 0.0
        RK3 = 0.0; RK3L = 0.0; RK3P = 0.0
        RK4 = 0.0; RK4L = 0.0; RK4P = 0.0
        RK5 = 0.0; RK5L = 0.0; RK5P = 0.0
        RK6 = 0.0; RK6L = 0.0; RK6P = 0.0
        RK7 = 0.0; RK7L = 0.0; RK7P = 0.0
        RK8 = 0.0; RK8L = 0.0; RK8P = 0.0
        RK9 = 0.0; RK9L = 0.0; RK9P = 0.0
    end subroutine resetAbsorptionGridValues

end module Mesh
