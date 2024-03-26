module Shapes
    use Kinds
    use Constants
    use ShapeVars, only: VI, SLL, gammaPT, molType
    implicit none
contains
    real function lorentz(X)
        real, intent(in) :: X
        !*---------------------------------------------------------------------*
        !* For this shape the following variables are used:			           *
        !* X	- distance from line center VI ( X = V - VI cm**-1 ),	       *
        !* ALAL (gammaPT) - (Lorentz half width)**2,								       *
        !* SLL	- (intensity*density*half_width*VVH_factor)/pi (see LBL93),    *
        !* MOTYPE - type of the molecule: 2 - CO2, 1 - H2O, 0 - other,         *
        !* --------------------------------------------------------------------*
        lorentz = SLL / (X*X + gammaPT*gammaPT)
        !* CO2 far wing correction *
        select case (molType)
        case(2)
            if (VI>3750 .and. VI<4700. .and. abs(X)>3.0) THEN  ! <==== Tonkov et al. 3800-4700 cm-1 
                if (abs(X) <= 150.) then
                    lorentz = lorentz * 1.084 * exp(-0.027*abs(X))
                else
                    lorentz = lorentz * 0.208 * exp(-0.016*abs(X))
                end if
            end if
        end select
    end function lorentz

    real function doppler(X)
        real, intent(in) :: X
    end function doppler

    real function voigt(X)
        real, intent(in) :: X
    end function voigt
end module Shapes
