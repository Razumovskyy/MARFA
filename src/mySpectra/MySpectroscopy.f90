module MySpectroscopy
    use MyConstants, only: DP, pi, dopplerCONST
    use MyAtmosphere
    implicit none
    real(kind=DP) :: shiftedLinePosition ! cm-1, (lineWV + deltaForeign * pressure) -- shifted line center
    real :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM
    real :: deltaForeign ! cm^-1/atm (pressure shift of the line position at 296 K and 1 atm)
    real(kind=DP) :: lineWV ! [cm-1] -- current spectral line
contains

    real(kind=DP) function shiftedLinePositionFunc()

        shiftedLinePositionFunc = lineWV + deltaForeign * pressure 
    end function shiftedLinePositionFunc

    real function DopplerHWHM()

        DopplerHWHM = dopplerCONST * shiftedLinePosition * sqrt(temperature/molarMass) 
    end function DopplerHWHM

end module MySpectroscopy
