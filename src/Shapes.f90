module Shapes
    use Constants
    use Atmosphere
    use Spectroscopy

    implicit none
contains
    ! TODO: consider removing various "lorentz" and keep one but with parameters
    ! real function simpleLorentz(X)
    !     !! 
    !     ! Simplistic Lorentz line shape function in which there is no account for self-broadening and no temperature dependence
    !     !!
    !     ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
    !     real(kind=DP), intent(in) :: X
    !     ! -------------------------------------------------------- !
    !     real(kind=DP) :: HWHM

    !     HWHM = lorentzHWHM(pressure, includeGammaSelf=.false., includeTemperature=.false.)
    !     simpleLorentz = HWHM / (pi*(X**2 + HWHM**2))
    ! end function simpleLorentz

    ! real function selfSimpleLorentz(X)
    !     !! 
    !     ! Lorentz line shape function in which there is no temperature dependence (T = 296 K), but self-broadening is taken into account
    !     !!
    !     ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
    !     real(kind=DP), intent(in) :: X
    !     ! -------------------------------------------------------- !
    !     real(kind=DP) :: HWHM

    !     HWHM = lorentzHWHM(pressure, includeGammaSelf=.true., partialPressureParameter=pSelf, includeTemperature=.false.)
    !     selfSimpleLorentz = HWHM / (pi*(X**2 + HWHM**2))
    ! end function selfSimpleLorentz

    ! real function noSelfLorentz(X)
    !     !! 
    !     ! Temperature-dependent Lorentz line shape with no account for self-broadening
    !     !!
    !     ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
    !     real(kind=DP), intent(in) :: X
    !     ! -------------------------------------------------------- !
    !     real(kind=DP) :: HWHM

    !     HWHM = lorentzHWHM(pressure, includeGammaSelf=.false., includeTemperature=.true., temperatureParameter=temperature)
    !     noSelfLorentz = HWHM / (pi*(X**2 + HWHM**2))
    ! end function noSelfLorentz

    real function lorentz(X)
        !! 
        ! Full-scale Lorentz shape calculation
        ! Temperature-dependent Lorentz line shape with self-broadening included
        !!
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM

        HWHM = lorentzHWHM(pressure, includeGammaSelf=.true., partialPressureParameter=pSelf, & 
                            includeTemperature=.true., temperatureParameter=temperature)
        lorentz = HWHM / (pi*(X**2 + HWHM**2))
        lorentz = lorentz * density * intensityOfT(temperature) 
    end function lorentz

    real function doppler(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: HWHM ! [cm-1] -- Doppler HWHM

        HWHM = dopplerHWHM(lineWV, temperature, molarMass)      
        doppler = sqrt(log(2.) / (pi * HWHM**2)) * exp(-(X/HWHM)**2 * log(2.))
        doppler = doppler * density * intensityOfT(temperature)
    end function doppler

    real function tonkov(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: shiftedLineWV

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        tonkov = lorentz(X)
        if (molType == 2) then
            if (shiftedLineWV > 3750 .and. shiftedLineWV < 4700. .and. abs(X) > 3.) then
                if (abs(X) <= 150.) then
                    tonkov = lorentz(X) * 1.084 * exp(-0.027*abs(X))
                else
                    tonkov = lorentz(X) * 0.208 * exp(-0.016*abs(X))
                end if
            end if
            return
        end if
    end function tonkov

    real function correctedDoppler(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        real, parameter :: bound = 12.5 * 12.5
        real, parameter :: U(9) = [1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
        real, parameter :: W(9) = [-0.688, 0.2667, 0.6338, 0.4405, 0.2529, 0.1601, 0.1131, 0.0853, 0.068]

        real(kind=DP) :: dopHWHM, lorHWHM
        real(kind=DP) :: XX
        real :: XI
        real :: F
        integer :: I
        !* ---------------------------------------------------------------*
        !*		For DOPLER are used :								      *
        !* X - distance from line center VI ( X = V - VI cm**-1 )	      *
        !* ADD - (Doppler half width)/(Ln(2))**0.5					      *
        !* SL - (intensity*density*VVH_factor)/pi (see LBL93.f lbf93.f)   *
        !* ---------------------------------------------------------------*
        dopHWHM = dopplerHWHM(lineWV, temperature, molarMass)
        lorHWHM = lorentzHWHM(pressure, includeGammaSelf=.true., partialPressureParameter=pSelf, &
                                includeTemperature=.true., temperatureParameter=temperature)
        XX = (X / (dopHWHM/sqln2))**2
        if (XX < bound) then
            !* Lorentz correction *
            if (XX >= 25.) then
                correctedDoppler = lorentz(X) * (1.+1.5/XX)
                return
            end if
            !* Pure Doppler *! -- the same is my doppler function
            correctedDoppler = (exp(-XX) / sqrt(pi) / (dopHWHM/sqln2)) * density * intensityOfT(temperature)
            if (XX <= 1.4) return
            !* Doppler shape begins to transform to Lorentz shape *
            XI = abs(X / (dopHWHM/sqln2))
            I = XI/0.5 - 1.00001
            F = 2. * (W(I)*(U(I+1)-XI) + W(I+1)*(XI-U(I)))
            correctedDoppler = correctedDoppler + lorHWHM/(pi*X**2) * (1.+F)
        else
            !* Doppler must be changed to Lorentz *
            correctedDoppler = lorentz(X)
            return
        end if
    end function correctedDoppler

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
        save A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, C3, C4, C5, C6, D3, D4, D5, D6, E5, E6
        
        dopHWHM = dopplerHWHM(lineWV, temperature, molarMass)
        lorHWHM = lorentzHWHM(pressure, includeGammaSelf=.true., partialPressureParameter=pSelf, &
                                includeTemperature=.true., temperatureParameter=temperature)
        
        YY = lorHWHM / (dopHWHM/sqln2)
        XX = abs(X / (dopHWHM/sqln2))

        if (XX > 15.) then	! Kuntz
            voigt = lorentz(X)
            return
        end if

        X2 = XX ** 2
        
        if (XX + YY >= 15.0) then 
            ! Region 1
            if (YY /= Y1) then
                Y1 = YY
                Y_2 = Y1 ** 2
                A1 = (0.2820948 + 0.5641896*Y_2) * Y1
                B1 = 0.5641896 * Y1
                A2 = 0.25 + Y_2 + Y_2**2
                B2 = Y_2 + Y_2 - 1.
            end if
            voigt = (A1+B1*X2) / (A2+B2*X2+X2**2) / sqrt(pi) / (dopHWHM/sqln2) * density * intensityOfT(temperature)
        else
            ! Region 2
            if (XX + YY >= 5.5) then
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
                    density * intensityOfT(temperature)
            else
                ! IF(Y >= 0.195*X-0.176)THEN ! Region 3 in accordance with Kuntz
                if (XX <= 1.0 .OR. YY >= 0.02) then 
                    ! Region 3 - Fomin suggestion
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
                             density * intensityOfT(temperature)
                else					
                    ! Region 4
                    voigt = correctedDoppler(X)
                end if
            end if
        end if
    end function voigt
end module Shapes

