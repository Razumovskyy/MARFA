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
    real :: A

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
        RK = 0.0; A = 0.0
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
        integer :: I, J, M

        !*						*** SUMMARISING		***
			DO J = 1,NT0
                I=J*2-1
                RK1P(I)=RK1P(I)+RK0P(J)
                !### RK1(I) =RK1(I)+RK0P(J)*0.375+RK0(J)*0.75-RK0L(J)*0.125
    A=RK0P(J)*0.375+RK0(J)*0.75-RK0L(J)*0.125
    IF(A<RK0P(J).OR.A<RK0(J))A=(RK0P(J)+RK0(J))*0.5
    RK1(I) =RK1(I)+A
                RK1L(I)=RK1L(I)+RK0(J)
                M=I+1
                RK1P(M)=RK1P(M)+RK0(J)
                !### RK1(M) =RK1(M)+RK0L(J)*0.375+RK0(J)*0.75-RK0P(J)*0.125
    A=RK0L(J)*0.375+RK0(J)*0.75-RK0P(J)*0.125
    IF(A<RK0L(J).OR.A<RK0(J))A=(RK0L(J)+RK0(J))*0.5
    RK1(M) =RK1(M)+A
                RK1L(M)=RK1L(M)+RK0L(J)
                END DO
    !*
                DO J = 1,NT1
                I=J*2-1
                RK2P(I)=RK2P(I)+RK1P(J)
                !### RK2(I) =RK2(I)+RK1P(J)*0.375+RK1(J)*0.75-RK1L(J)*0.125
    A=RK1P(J)*0.375+RK1(J)*0.75-RK1L(J)*0.125
    IF(A<RK1P(J).OR.A<RK1(J))A=(RK1P(J)+RK1(J))*0.5
    RK2(I) =RK2(I)+A
                RK2L(I)=RK2L(I)+RK1(J)
                M=I+1
                RK2P(M)=RK2P(M)+RK1(J)
                !###RK2(M) =RK2(M)+RK1L(J)*0.375+RK1(J)*0.75-RK1P(J)*0.125
    A=RK1L(J)*0.375+RK1(J)*0.75-RK1P(J)*0.125
    IF(A<RK1L(J).OR.A<RK1(J))A=(RK1L(J)+RK1(J))*0.5
    RK2(M) =RK2(M)+A
                RK2L(M)=RK2L(M)+RK1L(J)
                    END DO
    !*
                DO J = 1,NT2
                I=J*2-1
                RK3P(I)=RK3P(I)+RK2P(J)
                !### RK3(I) =RK3(I)+RK2P(J)*0.375+RK2(J)*0.75-RK2L(J)*0.125
    A=RK2P(J)*0.375+RK2(J)*0.75-RK2L(J)*0.125
    IF(A<RK2P(J).OR.A<RK2(J))A=(RK2P(J)+RK2(J))*0.5
    RK3(I) =RK3(I)+A
                RK3L(I)=RK3L(I)+RK2(J)
                M=I+1
                RK3P(M)=RK3P(M)+RK2(J)
                !### RK3(M) =RK3(M)+RK2L(J)*0.375+RK2(J)*0.75-RK2P(J)*0.125
    A=RK2L(J)*0.375+RK2(J)*0.75-RK2P(J)*0.125
    IF(A<RK2L(J).OR.A<RK2(J))A=(RK2L(J)+RK2(J))*0.5
    RK3(M) =RK3(M)+A
                RK3L(M)=RK3L(M)+RK2L(J)
                    END DO
    !*
                DO J = 1,NT3
                I=J*2-1
                RK4P(I)=RK4P(I)+RK3P(J)
                !### RK4(I) =RK4(I)+RK3P(J)*0.375+RK3(J)*0.75-RK3L(J)*0.125
    A=RK3P(J)*0.375+RK3(J)*0.75-RK3L(J)*0.125
    IF(A<RK3P(J).OR.A<RK3(J))A=(RK3P(J)+RK3(J))*0.5
    RK4(I) =RK4(I)+A
                RK4L(I)=RK4L(I)+RK3(J)
                M=I+1
                RK4P(M)=RK4P(M)+RK3(J)
                !### RK4(M) =RK4(M)+RK3L(J)*0.375+RK3(J)*0.75-RK3P(J)*0.125
    A=RK3L(J)*0.375+RK3(J)*0.75-RK3P(J)*0.125
    IF(A<RK3L(J).OR.A<RK3(J))A=(RK3L(J)+RK3(J))*0.5
    RK4(M) =RK4(M)+A
                RK4L(M)=RK4L(M)+RK3L(J)
                    END DO
    !*
                DO J = 1,NT4
                I=J*2-1
                RK5P(I)=RK5P(I)+RK4P(J)
                !### RK5(I) =RK5(I)+RK4P(J)*0.375+RK4(J)*0.75-RK4L(J)*0.125
    A=RK4P(J)*0.375+RK4(J)*0.75-RK4L(J)*0.125
    IF(A<RK4P(J).OR.A<RK4(J))A=(RK4P(J)+RK4(J))*0.5
    RK5(I) =RK5(I)+A
                RK5L(I)=RK5L(I)+RK4(J)
                M=I+1
                RK5P(M)=RK5P(M)+RK4(J)
                !### RK5(M) =RK5(M)+RK4L(J)*0.375+RK4(J)*0.75-RK4P(J)*0.125
    A=RK4L(J)*0.375+RK4(J)*0.75-RK4P(J)*0.125
    IF(A<RK4L(J).OR.A<RK4(J))A=(RK4L(J)+RK4(J))*0.5
    RK5(M) =RK5(M)+A
                RK5L(M)=RK5L(M)+RK4L(J)
                    END DO
    !*
                DO J = 1,NT5
                I=J*2-1
                RK6P(I)=RK6P(I)+RK5P(J)
                !###  RK6(I) =RK6(I)+RK5P(J)*0.375+RK5(J)*0.75-RK5L(J)*0.125
    A=RK5P(J)*0.375+RK5(J)*0.75-RK5L(J)*0.125
    IF(A<RK5P(J).OR.A<RK5(J))A=(RK5P(J)+RK5(J))*0.5
    RK6(I) =RK6(I)+A
                RK6L(I)=RK6L(I)+RK5(J)
                M=I+1
                RK6P(M)=RK6P(M)+RK5(J)
                !###  RK6(M) =RK6(M)+RK5L(J)*0.375+RK5(J)*0.75-RK5P(J)*0.125
    A=RK5L(J)*0.375+RK5(J)*0.75-RK5P(J)*0.125
    IF(A<RK5L(J).OR.A<RK5(J))A=(RK5L(J)+RK5(J))*0.5
    RK6(M) =RK6(M)+A
                RK6L(M)=RK6L(M)+RK5L(J)
                    END DO
                DO J = 1,NT6
                I=J*2-1
                RK7P(I)=RK7P(I)+RK6P(J)
                !### RK7(I) =RK7(I)+RK6P(J)*0.375+RK6(J)*0.75-RK6L(J)*0.125
    A=RK6P(J)*0.375+RK6(J)*0.75-RK6L(J)*0.125
    IF(A<RK6P(J).OR.A<RK6(J))A=(RK6P(J)+RK6(J))*0.5
    RK7(I) =RK7(I)+A
                RK7L(I)=RK7L(I)+RK6(J)
                M=I+1
                RK7P(M)=RK7P(M)+RK6(J)
                !### RK7(M) =RK7(M)+RK6L(J)*0.375+RK6(J)*0.75-RK6P(J)*0.125
    A=RK6L(J)*0.375+RK6(J)*0.75-RK6P(J)*0.125
    IF(A<RK6L(J).OR.A<RK6(J))A=(RK6L(J)+RK6(J))*0.5
    RK7(M) =RK7(M)+A
                RK7L(M)=RK7L(M)+RK6L(J)
                    END DO
    !*
                DO J = 1,NT7
                I=J*2-1
                RK8P(I)=RK8P(I)+RK7P(J)
                !### RK8(I) =RK8(I)+RK7P(J)*0.375+RK7(J)*0.75-RK7L(J)*0.125
    A=RK7P(J)*0.375+RK7(J)*0.75-RK7L(J)*0.125
    IF(A<RK7P(J).OR.A<RK7(J))A=(RK7P(J)+RK7(J))*0.5
    RK8(I) =RK8(I)+A
                RK8L(I)=RK8L(I)+RK7(J)
                M=I+1
                RK8P(M)=RK8P(M)+RK7(J)
                !### RK8(M) =RK8(M)+RK7L(J)*0.375+RK7(J)*0.75-RK7P(J)*0.125
    A=RK7L(J)*0.375+RK7(J)*0.75-RK7P(J)*0.125
    IF(A<RK7L(J).OR.A<RK7(J))A=(RK7L(J)+RK7(J))*0.5
    RK8(M) =RK8(M)+A
                RK8L(M)=RK8L(M)+RK7L(J)
                    END DO
    !*
                DO J = 1,NT8
                I=J*2-1
                RK9P(I)=RK9P(I)+RK8P(J)
                !### RK9(I) =RK9(I)+RK8P(J)*0.375+RK8(J)*0.75-RK8L(J)*0.125
    A=RK8P(J)*0.375+RK8(J)*0.75-RK8L(J)*0.125
    IF(A<RK8P(J).OR.A<RK8(J))A=(RK8P(J)+RK8(J))*0.5
    RK9(I) =RK9(I)+A
                RK9L(I)=RK9L(I)+RK8(J)
                M=I+1
                RK9P(M)=RK9P(M)+RK8(J)
                !### RK9(M) =RK9(M)+RK8L(J)*0.375+RK8(J)*0.75-RK8P(J)*0.125
    A=RK8L(J)*0.375+RK8(J)*0.75-RK8P(J)*0.125
    IF(A<RK8L(J).OR.A<RK8(J))A=(RK8L(J)+RK8(J))*0.5
    RK9(M) =RK9(M)+A
                RK9L(M)=RK9L(M)+RK8L(J)
                    END DO
    !*
                I=1
                DO J = 1,NT9
                I=I+1
        !### RK(I) =RK(I)+(RK9P(J)*0.375+RK9(J)*0.75-RK9L(J)*0.125)
    A=RK9P(J)*0.375+RK9(J)*0.75-RK9L(J)*0.125
    IF(A<RK9P(J).OR.A<RK9(J))A=(RK9P(J)+RK9(J))*0.5
    RK(I) =RK(I)+A
                I=I+1
        RK(I)=RK(I)+RK9(J)
                I=I+1
        !### RK(I) =RK(I)+(RK9L(J)*0.375+RK9(J)*0.75-RK9P(J)*0.125)
    A=RK9L(J)*0.375+RK9(J)*0.75-RK9P(J)*0.125
    IF(A<RK9L(J).OR.A<RK9(J))A=(RK9L(J)+RK9(J))*0.5
    RK(I) =RK(I)+A
                I=I+1
        RK(I)=RK(I)+RK9L(J)
                    END DO
    end subroutine cascadeInterpolation
end module Grids
