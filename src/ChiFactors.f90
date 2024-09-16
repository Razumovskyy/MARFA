module ChiFactors
    use Constants
    use ShapeFuncInterface
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
        ! case ('pollack')
        !     shapeFuncPtr => pollack
        ! case ('bezard')
        !     shapeFuncPtr => bezard
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
        tonkov = 1. ! default value of chi-factor
        if (molType == 2) then
            if (shiftedLineWV > 3750 .and. shiftedLineWV < 4700. .and. abs(X) > 3.) then
                if (abs(X) <= 150.) then
                    tonkov = 1.084 * exp(-0.027*abs(X))
                else
                    tonkov = 0.208 * exp(-0.016*abs(X))
                end if
            end if
            return
        end if
    end function tonkov
end module ChiFactors