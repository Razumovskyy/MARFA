module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy
    implicit none
contains

    ! all the shapes here are represented as (line shape function) x (temperature-dependent intensity)
    ! all the shape functions must satisfy abstract interface presented in `Interfaces.f90`
    
    ! add custom line shapes functions here !

    real function voigt(X)
        ! Implementation is based on the Humlicek method, improved by Kuntz:
        ! M. Kuntz. “A new implementation of the Humlicek algorithm for the calculation of the Voigt profile
        ! function”. In: Journal of Quantitative Spectroscopy and Radiative Transfer 57.6 (1997), pp. 819–824.

        ! Additional features and deviations:
        ! - scheme is not recursive
        ! - add region 0 where pure Lorentz is used
        ! - the boundaries between regions 3 and 4 are deviated from Kuntz's initial study.
        ! - in region 4 asymptotical analytical approximations are used for increasing speed
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !

        real :: VX, VY ! x and y parameters in the K(x,y) function
        
        real :: VXsquared ! x**2
        real :: Y1=0, Y2=0, Y3=0
        real :: Y_2
        real :: A1, B1, A2, B2, A3, B3, C3, D3, A4, B4, C4, D4, A5, B5, C5, D5, E5, &
                            A6, B6, C6, D6, E6
        
        ! TODO: figure out if this save can be removed (do when applying atmospheric level parallelization)
        save A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, C3, C4, C5, C6, D3, D4, D5, D6, E5, E6
        
        VX = abs(sqln2 * X / dopHWHM)
        VY = sqln2 * lorHWHM / dopHWHM

        ! REGION 0: Lorentz domination: pure Lorentz with χ-corrected wing (if set)
        if (VX >= 15.) then
            voigt = chiCorrectedLorentz(X, lorHWHM)
            return
        end if

        VXsquared = VX ** 2

        ! REGION 1: Voigt rational approximation 1
        if (VX + VY >= 15.0) then 

            if (VY /= Y1) then
                Y1 = VY
                Y_2 = Y1 ** 2
                A1 = (0.2820948 + 0.5641896*Y_2) * Y1
                B1 = 0.5641896 * Y1
                A2 = 0.25 + Y_2 + Y_2**2
                B2 = Y_2 + Y_2 - 1.
            end if
            ! rational approximation for Voigt using A1, A2, B1, B2 coefficients
            voigt = (A1+B1*VXsquared) / (A2+B2*VXsquared+VXsquared**2) / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 2: Voigt rational approximation 2       
        else if (VX + VY >= 5.5) then
            if (VY /= Y2) then
                Y2 = VY
                Y_2 = Y2**2
                A3 = Y2 * (((0.56419*Y_2+3.10304)*Y_2+4.65456)*Y_2+1.05786)
                B3 = Y2 * ((1.69257*Y_2+0.56419)*Y_2+2.962)
                C3 = Y2 * (1.69257*Y_2-2.53885)
                D3 = Y2*0.56419
                A4 = (((Y_2+6.0)*Y_2+10.5)*Y_2+4.5)*Y_2+0.5625
                B4 = ((4.0*Y_2+6.0)*Y_2+9.0)*Y_2-4.5
                C4 = 10.5+6.0*(Y_2-1.0)*Y_2
                D4 = 4.0*Y_2-6.0
            end if 

            voigt = (((D3*VXsquared+C3)*VXsquared+B3)*VXsquared+A3) / ((((VXsquared+D4)*VXsquared+C4)*VXsquared+B4)*VXsquared+A4) &
                            / sqrt(pi) / (dopHWHM/sqln2) * lineIntensity
        
        ! REGION 3: Voigt rational approximation 3
        else if (VX <= 1.0 .OR. VY >= 0.02) then 
            if (VY /= Y3) then
                Y3 = VY
                A5 = ((((((((0.564224*Y3+7.55895)*Y3+49.5213)*Y3+204.510)*Y3+	&
                    581.746)*Y3+1174.8)*Y3+1678.33)*Y3+1629.76)*Y3+973.778)*Y3+272.102
                B5 = ((((((2.25689*Y3+22.6778)*Y3+100.705)*Y3+247.198)*Y3+336.364)*	&
                    Y3+220.843)*Y3-2.34403)*Y3-60.5644
                C5 = ((((3.38534*Y3+22.6798)*Y3+52.8454)*Y3+42.5683)*Y3+18.546)*Y3+	&
                    4.58029
                D5 = ((2.25689*Y3+7.56186)*Y3+1.66203)*Y3-0.128922
                E5 = 0.971457E-3+0.564224*Y3
                A6 = (((((((((Y3+13.3988)*Y3+88.2674)*Y3+369.199)*Y3+1074.41)*Y3+	&
                    2256.98)*Y3+3447.63)*Y3+3764.97)*Y3+2802.87)*Y3+1280.83)*Y3+	&
                    272.102
                B6 = (((((((5.*Y3+53.5952)*Y3+266.299)*Y3+793.427)*Y3+1549.68)*Y3+	&
                    2037.31)*Y3+1758.34)*Y3+902.306)*Y3+211.678
                C6 = (((((10.*Y3+80.3928)*Y3+269.292)*Y3+479.258)*Y3+497.302)*Y3+	&
                    308.186)*Y3+78.866
                D6 = (((10.*Y3+53.5952)*Y3+92.7586)*Y3+55.0293)*Y3+22.0353
                E6 = (5.0*Y3+13.3988)*Y3+1.49645
            end if
            
            voigt = ((((E5*VXsquared+D5)*VXsquared+C5)*VXsquared+B5)*VXsquared+A5)/	&
                    (((((VXsquared+E6)*VXsquared+D6)*VXsquared+C6)*VXsquared+B6)*VXsquared+A6) / sqrt(pi) &
                            / (dopHWHM/sqln2) * lineIntensity
    
        else
            ! REGION 4: ASYMPTOTIC REGION: (1. < VX < 5.5 and y < 0.02)
            ! where instead of approximation of Voigt function, 
            ! can be used analytical asymptotical expressions based on series expansions 
            ! of Lorentz and Doppler shape functions	
            if (VX > 5.) then
                voigt = voigtAsymptotic1(X, lorHWHM, VX)
            else if (VX > sqrt(1.4)) then
                voigt = voigtAsymptotic2(X, lorHWHM, VX, dopHWHM)
            else 
                voigt = doppler(X, dopHWHM)
            end if

        end if
    end function voigt

end module Shapes
