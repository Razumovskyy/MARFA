module Grids
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
    
    ! TODO: move it out from parameters, it might be an input value, or
    ! a dictionary value based on the input
    real, parameter :: deltaWV = 10.0 ! resolution in cm-1
    
    ! TODO: deal with it during refactor of the grid calculation
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


    subroutine cascadeInterpolation
        implicit none
        integer :: I, J
    
        ! Perform cascading interpolation steps from RK0 -> RK1, RK1 -> RK2, ... RK8 -> RK9
        call interpolationStep(NT0, RK0P, RK0, RK0L, RK1P, RK1, RK1L)
        call interpolationStep(NT1, RK1P, RK1, RK1L, RK2P, RK2, RK2L)
        call interpolationStep(NT2, RK2P, RK2, RK2L, RK3P, RK3, RK3L)
        call interpolationStep(NT3, RK3P, RK3, RK3L, RK4P, RK4, RK4L)
        call interpolationStep(NT4, RK4P, RK4, RK4L, RK5P, RK5, RK5L)
        call interpolationStep(NT5, RK5P, RK5, RK5L, RK6P, RK6, RK6L)
        call interpolationStep(NT6, RK6P, RK6, RK6L, RK7P, RK7, RK7L)
        call interpolationStep(NT7, RK7P, RK7, RK7L, RK8P, RK8, RK8L)
        call interpolationStep(NT8, RK8P, RK8, RK8L, RK9P, RK9, RK9L)
    
        ! Handle final step: RK9 -> RK
        I=1
        do J = 1, NT9
            I = I + 1
            RK(I) = RK(I) + (RK9P(J)*0.375 + RK9(J)*0.75 - RK9L(J)*0.125)
            I = I + 1
            RK(I) = RK(I) + RK9(J)
            I = I + 1
            RK(I) = RK(I) + (RK9L(J)*0.375 + RK9(J)*0.75 - RK9P(J)*0.125)
            I = I + 1
            RK(I) = RK(I) + RK9L(J)
        end do
    
        ! Remove negative outliers
        do J = 1, NT
            if (RK(J) < 0.0) RK(J)=0.0
        end do
    
    end subroutine cascadeInterpolation
    

    subroutine interpolationStep(NTsize, prevP, prev, prevL, nextP, next, nextL)
        !--------------------------------------------------------------------
        ! Helper subroutine to perform a single interpolation step (simple parabolic)
        ! It takes arrays from "previous" more coarse grid (prevP, prev, prevL) and
        ! updates the "next" more fine grid arrays (nextP, next, nextL).
        !--------------------------------------------------------------------
        implicit none
        integer, intent(in) :: NTsize
        real, intent(in) :: prevP(:), prev(:), prevL(:)
        real, intent(inout) :: nextP(:), next(:), nextL(:)
    
        integer :: J, I, M
    
        do J = 1, NTsize
            I = J*2 - 1
            nextP(I) = nextP(I) + prevP(J)
            next(I)  = next(I)  + prevP(J)*0.375 + prev(J)*0.75 - prevL(J)*0.125
            nextL(I) = nextL(I) + prev(J)
            M = I + 1
            nextP(M) = nextP(M) + prev(J)
            next(M)  = next(M)  + prevL(J)*0.375 + prev(J)*0.75 - prevP(J)*0.125
            nextL(M) = nextL(M) + prevL(J)
        end do
    
    end subroutine interpolationStep
end module Grids
