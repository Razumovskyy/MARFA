module LBL
    use Constants
    use Interfaces
    use Atmosphere
    use Mesh
    use IO
    use Spectroscopy
    use MolarMasses
    use Shapes
    use LineGridCalc
    implicit none
    integer, parameter :: databaseFileUnit = 7777
contains

    subroutine modernLBL(lineIdxParameter, capWVParameter)
        
        integer, intent(inout) :: lineIdxParameter ! integer line label used for locating record in direct access file
        real(kind=DP), intent(inout) :: capWVParameter ! inout parameter for capWV. For definition see main.f90 file

        real(kind=DP) :: shiftedLineWV ! VI ! pressure-induced shifted line transition wavenumber

        ! techincal variables
        integer :: I ! loop variable for accessing the record in direct access file
        integer :: ios2

    !------------------------------------------------------------------------------------------------------------------!

        ! line-by-line loop: reading spectral data for the subinterval and performing summation
        I = lineIdxParameter - 1
        ios2 = 0
        LBL_LOOP: do while (.not. is_iostat_end(ios))
            I = I + 1
            read(databaseFileUnit, rec=I, iostat=ios2) lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, & 
                                                    foreignTempCoeff, jointMolIso, deltaForeign
     
            ! DEBUG SECTION !
            ! print *, 'lineWV: ', lineWV
            ! print *, 'refLineIntensity: ', refLineIntensity
            
            if (ios2 > 0) then
                print *, 'ERROR: when reading file with spectral data.'
                stop 9
            end if

            ! exit the loop when extended subinterval (endDeltaWV + cutOff) ends
            ! extEndDeltaWV = endDeltaWV + cutOff
            if  (lineWV >= extEndDeltaWV) exit

            ! setting the lineIdxParameter for the next call -- for corrct switching to the next subinterval
            ! capWVParameter = endDeltaWV - cutOff, so capWVParameter < extEndDeltaWV
            if  (lineWV <= capWVParameter) lineIdxParameter = I

            ! Setting the molar weight
            isotopeNum = jointMolIso/100
            molarMass = WISO(isotopeNum)

            ! calculation of pressure-induced shift
            shiftedLineWV = shiftedLinePosition(lineWV, pressure)
            dopHWHM = dopplerHWHM(shiftedLineWV, temperature, molarMass)
            lorHWHM = lorentzHWHM(pressureParameter=pressure, partialPressureParameter=pSelf, &
                                    temperatureParameter=temperature)

            lineIntensity = intensityOfT(temperature)

            ! DEBUG SECTION !
            ! print *, 'dopHWHM: ', dopHWHM
            ! print *, 'lorHWHM: ', lorHWHM
            ! print *, 'VY: ', lorHWHM * sqln2 /dopHWHM
            ! pause
            
            if (shiftedLineWV < startDeltaWV) then
                ! in this case current LBL_LOOP spectral line falls into the interval: 
                ! [startDeltaWV - cutOff; startDeltaWV] and to find its contribution 
                ! on the [startDeltaWV; endDeltaWV] interval, only RIGHT WING of this line
                ! must be accounted for. It is done in leftLBL subroutine:
                call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr) 
            else if (shiftedLineWV > endDeltaWV) then
                ! in this case current LBL_LOOP spectral line falls into the interval: 
                ! [endDeltaWV; endDeltaWV + cutOff] and to find its contribution 
                ! on the [startDeltaWV; endDeltaWV] interval, only LEFT WING of this line
                ! must be accounted for. It is done in rightLBL subroutine:
                call rightLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
            else
                ! in this case current LBL_LOOP spectral line falls into the subinterval itself
                ! BOTH WINGS of this line must be accounted for.
                ! It is implemented in the centerLBL subroutine:
                call centerLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
            end if
        end do LBL_LOOP
    end subroutine modernLBL
end module LBL
