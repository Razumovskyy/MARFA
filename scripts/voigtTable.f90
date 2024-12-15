program voigtTable
    ! draft script for creating a Voigt function calculation table
    use Constants
    use Atmosphere
    use Interfaces
    use ChiFactors
    use Spectroscopy
    use Shapes
    implicit none

    integer, parameter :: voigtTableUnit = 344
    character, parameter :: voigtTableFile = 'voigtTable.dat'
    real(kind=DP), dimension(10), parameter :: vxValues = (/1D-3, 1D-2, 1D-1, 1.1D0, 3.0D0, 5.2D0, 10.0D0, 20.0D0, 1D2, 1D3/)

    character(len=300) :: header
    real(kind=DP) :: shiftedLineWV ! VI ! pressure-induced shifted line transition wavenumber
    real(kind=DP), dimension(10) :: xValues = (/0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/)
    real :: VY
    real, dimension(10) :: voigtValues = (/0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
    integer :: i

    ! line parameters:
    lineWV =    252.32362100000000
    gammaForeign =    6.47000000E-02
    gammaSelf =   7.10000023E-02
    foreignTempCoeff =  0.730000019
    deltaForeign =  2.05999997E-04
    jointMolIso = 1002
    reflineIntensity = 2.49100005E-30
    lineLowerState = 3456.05225
    
    atmProfileFile = 'VenusCO2.dat'
    call readAtmosphericParameters()
    call readTips

    header = 'HEIGHT        DOPHWHM         LORHWHM                VY                VX1       VX2     VX3        VX4        &
                           VX5        VX6         VX7        VX8            VX9     VX10'
    open(voigtTableUnit, file='voigtTable.dat', status='replace', action='write')
    write(voigtTableUnit, *) trim(adjustl(header))
    
    chiFactorFuncPtr => noneChi
    
    isotopeNum = jointMolIso/100
    molarMass = WISO(isotopeNum)
    do levelsIdx = 1, levels
        pressure = pressureArray(levelsIdx)
        density = densityArray(levelsIdx)
        temperature = temperatureArray(levelsIdx)
        lineIntensity = intensityOfT(temperature)
        print *, 'SL: ', lineIntensity * density / 3.14
        pause
        pSelf = density * 10. / LOSCHMIDT * temperature/stTemperature
        pForeign = pressure - pSelf
        shiftedLineWV = shiftedLinePosition(lineWV, pressure)
        dopHWHM = dopplerHWHM(shiftedLineWV, temperature, molarMass)
        lorHWHM = lorentzHWHM(pressure, pSelf, temperature)
        VY = lorHWHM * sqln2 / dopHWHM
        xValues = vxValues * dopHWHM / sqln2
        print *, xValues
        print *, levelsIdx
        do i = 1, 10
            voigtValues(i) = voigt(xValues(i))
        end do
        write(voigtTableUnit,  '(F10.4, 4X, E15.8E3, 4X, E15.8E3, 4X, E15.8E3, 4X, 10(4X, E15.8E3))') heightArray(levelsIdx), dopHWHM, lorHWHM, VY, voigtValues
    end do
    close(voigtTableUnit)
    print *, 'READY!'

end program
