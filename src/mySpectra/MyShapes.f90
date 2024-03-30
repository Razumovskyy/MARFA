module MyShapes
    use MyConstants
    use MyAtmosphere
    use MySpectroscopy
    implicit none 
contains
    
    real function simpleLorentz(nu)
        real(kind=DP), intent(in) :: nu ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated

        simpleLorentz = (gammaForeign*pressure) / (pi*((nu-shiftedLinePosition)**2 + (gammaForeign*pressure)**2))
    end function simpleLorentz

    real function simpleDoppler(nu)
        real(kind=DP), intent(in) :: nu ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated

        real :: refAlphaD ! cm-1 -- Doppler HWHM at reference temperature of 296 K (refTemperature)
        real(kind=DP) :: shiftedLineWV
        
        shiftedLineWV = shiftedLinePositionFunc()
        refAlphaD = DopplerHWHM()
        simpleDoppler = sqrt(log(2.) / (pi * refAlphaD**2)) * exp(-((nu - shiftedLinePosition)**2 * log(2.)) / (refAlphaD**2))
    end function simpleDoppler

end module MyShapes