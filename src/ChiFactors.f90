module ChiFactors
    use Constants
    use Interfaces
    use Atmosphere
    implicit none
    
    procedure(chifactor), pointer :: chiFactorFuncPtr
    character(len=30) :: chiFactorFuncName ! name of the chi-factor function
contains

    subroutine setChiFactorFunction()
        select case(trim(adjustl(chiFactorFuncName)))
        case ('none')
            chiFactorFuncPtr => noneChi
        case ('tonkov')
            chiFactorFuncPtr => tonkov
        case ('pollack')
            chiFactorFuncPtr => pollack
        case ('perrin')
            chiFactorFuncPtr => perrin
        case default
            print *, "unknown chi-factor function: ", trim(adjustl(chiFactorFuncName))
            stop 10
        end select
    end subroutine setChiFactorFunction
    
    real function noneChi(X, moleculeIntCode)
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !
        noneChi = 1.
    end function noneChi

    real function tonkov(X, moleculeIntCode)
        ! MV Tonkov et al. “Measurements and empirical modeling of pure CO2 absorption in the 2.3-μm
        ! region at room temperature: far wings, allowed and collision-induced bands”. In: Applied optics 35.24
        ! (1996), pp. 4863–4870.

        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        tonkov = 1. ! default value of the chi-factor
        if (moleculeIntCode == 2) then ! CO2
            if (abs(X) > 3.) then
                if (abs(X) <= 150.) then
                    tonkov = 1.084 * exp(-0.027*abs(X))
                else
                    tonkov = 0.208 * exp(-0.016*abs(X))
                end if
            end if
        end if
    end function tonkov

    real function pollack(X, moleculeIntCode)
        ! James B. Pollack et al. “Near-Infrared Light from Venus’ Nightside: A Spectroscopic Analysis”. In:
        ! Icarus 103 (1 1993), pp. 1–42. ISSN: 10902643. DOI: 10.1006/icar.1993.1055
        
        ! Pollack(1993):  recommended cutoff is 125 cm-1 or 160 cm-1 for the 1.2 um window complex (see p. 5.4 in the paper)
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !

        pollack = 1. ! default value of the chi-factor
        
        if (moleculeIntCode == 2) then
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

    real function perrin(X, moleculeIntCode)
        ! MY Perrin and JM Hartmann. “Temperature-dependent measurements and modeling of absorption
        ! by CO2-N2 mixtures in the far line-wings of the 4.3 μm CO2 band”. In: Journal of Quantitative
        ! Spectroscopy and Radiative Transfer 42.4 (1989), pp. 311–317.
        
        implicit none
        ! X - [cm-1] -- distance from the shifted line center to the spectral point of function evaluation
        real(kind=DP), intent(in) :: X
        integer, intent(in) :: moleculeIntCode
        ! -------------------------------------------------------- !
        real :: B1, B2, B3, S1, S2, S3

        perrin = 1. ! default value for the chi-factor

        ! temperature dependence
        B1 = 0.0888 - 0.16*exp(-0.00410*temperature)
        B2 = 0.0526 * exp(-0.00152*temperature)
        B3 = 0.0232

        ! boundaries
        S1 = 3.
        S2 = 30.
        S3 = 120.

        if (moleculeIntCode == 2) then 
            ! CO2
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
