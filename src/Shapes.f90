module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy
    use ChiFactors
    implicit none
contains

    ! all the shapes here are represented as (line shape function) x (temperature-dependent intensity)
    ! all the shape functions must satisfy abstract interface presented in `Interfaces.f90`
    
    ! add custom line shapes functions here !

    real function lorentz(X)

        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressureParameter=pressure, partialPressureParameter=pSelf, & 
                            temperatureParameter=temperature)
        lorentz = HWHM / (pi*(X**2 + HWHM**2))
        lorentz = lorentz * intensityOfT(temperature) 
    end function lorentz

    real function doppler(X)
        
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        real(kind=DP) :: HWHM ! [cm-1] -- Doppler HWHM

        HWHM = dopplerHWHM(lineWV, temperature, molarMass)      
        doppler = sqrt(log(2.) / (pi * HWHM**2)) * exp(-(X/HWHM)**2 * log(2.))
        doppler = doppler * intensityOfT(temperature)
    end function doppler

    real function chiCorrectedLorentz(X)
        
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X

        chiCorrectedLorentz = lorentz(X) * chiFactorFuncPtr(X)

    end function chiCorrectedLorentz

    ! TODO:(!) check if the return statements in else-if blocks affect calculation speed 

    real function voigt(X)
        
        ! ADD = dopHWHM / sqrt(ln2)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: dopHWHM, lorHWHM
        real(kind=DP) :: XX, YY, X2
        real(kind=DP) :: Y1=0, Y2=0, Y3=0
        real(kind=DP) :: Y_2
        real(kind=DP) :: A1, B1, A2, B2, A3, B3, C3, D3, A4, B4, C4, D4, A5, B5, C5, D5, E5, &
                            A6, B6, C6, D6, E6
        
        ! TODO: figure out if this save can be removed (do when applying atmospheric level parallelization)
        save A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, C3, C4, C5, C6, D3, D4, D5, D6, E5, E6
        
        dopHWHM = dopplerHWHM(lineWV, temperature, molarMass)
        lorHWHM = lorentzHWHM(pressureParameter=pressure, partialPressureParameter=pSelf, &
                                temperatureParameter=temperature)
        
        XX = abs(X / (dopHWHM/sqln2))        
        YY = lorHWHM / (dopHWHM/sqln2)

        ! REGION 0: Lorentz domination: pure Lorentz with χ-correction
        if (XX > 15.) then
            voigt = chiCorrectedLorentz(X)
            return
        end if

        X2 = XX ** 2

        ! REGION 1: Voigt rational approximation 1
        if (XX + YY >= 15.0) then 

            if (YY /= Y1) then
                Y1 = YY
                Y_2 = Y1 ** 2
                A1 = (0.2820948 + 0.5641896*Y_2) * Y1
                B1 = 0.5641896 * Y1
                A2 = 0.25 + Y_2 + Y_2**2
                B2 = Y_2 + Y_2 - 1.
            end if
            ! rational approximation for Voigt using A1, A2, B1, B2 coefficients
            voigt = (A1+B1*X2) / (A2+B2*X2+X2**2) / sqrt(pi) / (dopHWHM/sqln2) * &
                            intensityOfT(temperature)
            return
        
        ! REGION 2: Voigt rational approximation 2       
        else if (XX + YY >= 5.5) then
            if (YY /= Y2) then
                Y2 = YY
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

            voigt = (((D3*X2+C3)*X2+B3)*X2+A3) / ((((X2+D4)*X2+C4)*X2+B4)*X2+A4) / sqrt(pi) / (dopHWHM/sqln2) * &
                            intensityOfT(temperature)
            return
        
        ! REGION 3: Voigt rational approximation 3
        else if (XX <= 1.0 .OR. YY >= 0.02) then 
            if (YY /= Y3) then
                Y3 = YY
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
            
            voigt = ((((E5*X2+D5)*X2+C5)*X2+B5)*X2+A5)/	&
                    (((((X2+E6)*X2+D6)*X2+C6)*X2+B6)*X2+A6) / sqrt(pi) / (dopHWHM/sqln2) * &
                            intensityOfT(temperature)
            return
    
        ! REGION 4: ASYMPTOTIC REGION: (1. < XX < 5.5 and y < 0.02)
        ! where instead of estimating voigt function, 
        ! analytical asymptotical expressions based on Lorentz and Doppler shape functions might be used
        else					
            voigt = asymptoticDopplerLorentz(X)
        end if
    end function voigt

    real function asymptoticDopplerLorentz(X)
        
        ! used only for the Voigt function estimation
        
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        real, parameter :: U(9) = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
        real, parameter :: W(9) = [-0.688, 0.2667, 0.6338, 0.4405, 0.2529, 0.1601, 0.1131, 0.0853, 0.068]

        real(kind=DP) :: dopHWHM, lorHWHM
        real(kind=DP) :: XX, X2
        real :: XI
        real :: F
        integer :: I

        dopHWHM = dopplerHWHM(lineWV, temperature, molarMass)
        lorHWHM = lorentzHWHM(pressureParameter=pressure, partialPressureParameter=pSelf, &
                                temperatureParameter=temperature)
        XX = abs(X / (dopHWHM/sqln2))
        X2 = XX ** 2
        
        ! REGION 4a: Lorentz leading with Doppler-influenced asymptotic correction
        if (X2 >= 25.) then
            ! asymptotic approximation 1
            asymptoticDopplerLorentz = chiCorrectedLorentz(X) * (1.+1.5/XX)
            return
        
        ! REGION 4b: Doppler leading with Lorentz-influenced asymptotic correction
        else if (X2 >= 1.4) then
            ! asymptotic approximation 2
            XI = abs(X / (dopHWHM/sqln2))
            I = XI/0.5 - 1.00001
            F = 2. * (W(I)*(U(I+1)-XI) + W(I+1)*(XI-U(I)))
            asymptoticDopplerLorentz = doppler(X) + (lorHWHM/(pi*X**2) * (1.+F)) * intensityOfT(temperature)
            return
        
        ! REGION 4c: Dominant Doppler, Lorentz negligible
        else
            asymptoticDopplerLorentz = doppler(X)
            return
        end if
    end function asymptoticDopplerLorentz
end module Shapes
