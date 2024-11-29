module ChiFactors
    use Constants
    use Interfaces
    use Atmosphere
    use Spectroscopy
    implicit none
    procedure(shape), pointer :: chiFactorFuncPtr
    character(len=30) :: chiFactorFuncName ! name of the chi-factor function
contains

    subroutine fetchChiFactorFunction()
        ! TODO: add flow for the incorrect chi factor function input
        select case(trim(adjustl(chiFactorFuncName)))
        case ('none')
            chiFactorFuncPtr => noneChi
        case ('tonkov')
            chiFactorFuncPtr => tonkov
        case ('pollack')
            chiFactorFuncPtr => pollack
        case ('perrin')
            chiFactorFuncPtr => perrin
        end select
    end subroutine fetchChiFactorFunction
    
    real function noneChi(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        noneChi = 1.
    end function noneChi

    real function tonkov(X)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: shiftedLineWV
        
        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        tonkov = 1. ! default value of the chi-factor
        if (molType == 2) then ! CO2
            if (abs(X) > 3.) then
                if (abs(X) <= 150.) then
                    tonkov = 1.084 * exp(-0.027*abs(X))
                else
                    tonkov = 0.208 * exp(-0.016*abs(X))
                end if
            end if
            return
        end if
    end function tonkov

    real function pollack(X)
        ! Pollack(1993):  recommended cutoff is 125 cm-1 or 160 cm-1 for the 1.2 um window complex (see p. 5.4 in the paper)
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: shiftedLineWV

        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        pollack = 1. ! default value of the chi-factor

        if (molType == 2) then
            if (abs(X) >= 3.) then
                if (abs(X) < 10.) then
                    pollack = 1.35 * exp(-abs(X)/10.)
                else
                    pollack = 0.614 * exp(-abs(X)/47.)
                    ! pollack = 0.614 * exp(abs(X)/90.) ! for the 1.18 um window
                end if
            end if
        end if
    end function pollack 

    real function perrin(X)
        ! Perrin and Hartmann (1989): temperature dependent chi-factor; CO2-CO2 broadening
        ! X - [cm-1] -- distance from the shifted line center to the spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: X
        ! -------------------------------------------------------- !
        real(kind=DP) :: shiftedLineWV
        real :: B1, B2, B3, S1, S2, S3
        shiftedLineWV = shiftedLinePosition(lineWV, pressure)

        perrin = 1. ! default value for the chi-factor

        ! temperature dependence
        B1 = 0.0888 - 0.16*exp(-0.00410*temperature)
        B2 = 0.0526 * exp(-0.00152*temperature)
        B3 = 0.0232

        ! boundaries
        S1 = 3.
        S2 = 30.
        S3 = 120.

        if (molType == 2) then 
            if (abs(X) < S2) then
                perrin = exp(-B1 * (abs(X)-S1))
            else 
                if (abs(X) < S3) then
                    perrin = exp(-B1*(S2-S1) - B2*(abs(X)-S2))
                else
                    perrin = exp(-B1*(S2-S1) - B2*(S3-S2) - B3*(abs(X)-S3))
                end if
            end if 
        end if
    end function perrin
end module ChiFactors