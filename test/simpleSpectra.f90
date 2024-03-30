program simpleSpectra
    use SimpleShapes, only: simpleLorentz, shiftedLinePositionFunc
    implicit none
    
    integer, parameter :: DP = selected_real_kind(15, 307)
    integer, parameter :: hitranFileUnit = 7777
    integer, parameter :: outputFileUnit = 7778
    character(len=20), parameter :: hitranFile = 'data/HITRAN16/H16.01'
    character(len=27), parameter :: outputFile = 'test/calc/simpleSpectra.dat'

    real(kind=DP) :: startWV, endWV, step
    real(kind=DP) :: lineWV ! cm-1 (wavenumber of the transition)
    real :: refLineIntensity ! cm-1/(molecule*cm-2) (spectral line intensity at 296 K)
    real :: gammaForeign, gammaSelf ! cm-1/atm (air- and self-broadened HWHM at 296 K)
    real :: lineLowerState ! cm-1 (lower state energy of the transition)
    real :: foreignTempCoeff ! dimensionless (coeff for temp dependence of gammaForeign)
    integer :: jointMolIso ! dimensionless (joined reference to Molecule number (MOL) and Isotopologue number (ISO))
    real :: deltaForeign ! cm^-1/atm (pressure shift of the line position at 296 K and 1 atm)

    real :: pressure ! atm
    real :: density ! molecules/(cm^2*km)
    real :: lineCutOff ! 

    real(kind=DP), allocatable :: spectra(:,:) ! 2D array: first item is grid wavenumber in cm^-1 and the second is the absorption coefficient in km^-1
    real(kind=DP) :: gridWV ! more human-readable variable for the x-values of the spectra array
    real(kind=DP) :: gridAbsorptionCoeff ! more human readable variable for the y-values of the spectra array
    integer :: len ! len of the spectra array
    integer :: i, j, jStart, k ! loop variables

    ! ---------------------------------------------------------------------------- !
    ! INPUT PARAMETERS !
    pressure = 1. ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
    density = 0.6994E+20 ! ~ Venus 50 km level (see data/Atmospheres/H2O_117.dat)
    startWV = 100.
    endWV = 110.
    step = 0.01
    lineCutOff = 25.
    ! ---------------------------------------------------------------------------- !   
    len = int((endWV-startWV) / step) + 1

    allocate(spectra(len, 2))

    lineWV = startWV

    jStart = 1
    do
        open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
        read(hitranFileUnit, rec=jStart) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                        lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign
        if (lineWV > startWV - lineCutOff) exit
        jStart = jStart + 1
    end do
    

    do i = 1, len
        write(*,*) i, ' of ', len, ' is processed'
        spectra(i, 1) = startWV + (i-1) * step
        spectra(i, 2) = 0.
        gridWV = spectra(i, 1)
        j = jStart
        do
            open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
            read(hitranFileUnit, rec=j) lineWV, refLineIntensity, gammaForeign, gammaSelf, &
                                    lineLowerState, foreignTempCoeff, jointMolIso, deltaForeign                      
            if (lineWV >= endWV + lineCutOff) exit
            if (abs(lineWV - gridWV) < lineCutOff) then
                gridAbsorptionCoeff = gridAbsorptionCoeff + &
                                    simpleLorentz(gridWV, refLineIntensity, lineWV, gammaForeign, pressure, density)
            end if
            j = j + 1
        spectra(i, 1) = gridWV
        spectra(i, 2) = gridAbsorptionCoeff
        end do
    end do

    open(outputFileUnit, file=outputFile, status='replace', action='write')

    do k = 1, len
        write(outputFileUnit, '(F8.3, ", ", E10.3)') spectra(k,1), spectra(k, 2)
    end do
    
    close(hitranFileUnit)
    close(outputFileUnit)
    deallocate(spectra)

end program simpleSpectra

module SimpleShapes
    implicit none
    integer, parameter :: DP = selected_real_kind(15, 307)
    integer, parameter :: pi = 3.14159
contains
    real(kind=DP) function shiftedLinePositionFunc(lineWV, deltaForeign, pressure)
        real(kind=DP), intent(in) :: lineWV ! cm-1 (lineWV) -- wavenumber of the transition (read from HITRAN or HITEMP)
        real, intent(in) :: deltaForeign ! cm-1/atm (deltaForeign) -- pressure shift coefficient (read from HITRAN or HITEMP)
        real, intent(in) :: pressure  ! [atm] -- pressure at the given atmospheric level 

        shiftedLinePositionFunc = lineWV + deltaForeign * pressure 
    end function shiftedLinePositionFunc

    real function simpleLorentz(nu, intensity, shiftedLinePosition, gammaForeign, pressure, density)
        real(kind=DP), intent(in) :: nu ! cm-1, (gridWV) -- spectral point in which the total contribution from lines is calculated
        real(kind=DP), intent(in) :: shiftedLinePosition ! cm-1, (lineWV + deltaForeign * pressure) -- shifted line center
        real, intent(in) :: intensity ! [cm-1/(molecule*cm-2)] (refLineIntensity) -- the spectral line intensity for a single molecule per unit volume.
        real, intent(in) :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM
        real, intent(in) :: pressure ! [atm] -- pressure at the given atmospheric level 
        real, intent(in) :: density ! [molecule/(cm^2 * km)] -- such density units are needed for having absorption coefficient in km-1
        
        simpleLorentz = (intensity*gammaForeign*pressure*density) / (pi*((nu-shiftedLinePosition)**2 + (gammaForeign*pressure)**2))
    end function simpleLorentz
end module SimpleShapes