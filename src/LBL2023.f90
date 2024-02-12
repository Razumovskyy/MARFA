module LBL
    use kinds
    use constants, only: PI, refTemperature, stPressure, stTemperature, LOSCHMIDT, C2
    use shared_vars_main
    use mesh
    use molecule_vars
    use spectroscopy, only: lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, foreignTempCoeff, &
                            jointMolIso, deltaForeign, lineIntensityofT
    use spectroscopy, only: cloughRFunc, LorentzianHWHM, DopplerHWHM, VOIGT, LORENTZ, DOPLER, pressureShiftedWV
    use molar_masses, only: WISO
    use shape_functions
    use LINE_GRID_CALC
    implicit none
contains

    subroutine LBL2023(molType, LINBEG, capWV, totalLines, loopLevel) ! VFISH is from shared_vars_main (extEndDeltaWV)
        
        procedure(shape_func), pointer :: ShapeFuncPtr  

        integer :: molType ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
        integer :: LINBEG ! integer line label used for locating record in direct access file
        integer :: I ! loop variable for accessing the record in direct access file
        integer :: totalLines ! number of lines in one file
        integer :: loopLevel
        integer, save :: currentLevel

        ! this variable serves as a cap in the loop to stop increasing the lineIdx
        ! it is the value as when there is the next step over the DeltaWV, hitran reading will from the correct record number
        ! the definition and the value assignment is in K_HTRAN module
        real(kind=DP) :: capWV

        real, parameter :: BOUNDL = 10. ! some boundary parameter
        real, parameter :: BOUNDD = 0.01 ! likely some boundary value related to Doppler broadening, given its small value

        real, save :: unitlessT ! unitless temperature (refTemperature/T)

        real :: lineIntensity ! QDQ * PI / RO ! temperature and pressure dependent line intensity
        
        real, save :: pSelf ! self-broadening pressure component
        real, save :: pForeign ! foreign gas broadening component
        
        real(kind=DP) :: gammaPT ! AL ! Lorentz HWHM -- T and P dependent !
        
        ! Appears in the denominator of the Lorentzian line profile L(\nu)
        real(kind=DP) :: shiftedLineWV ! VI ! shifted line position under the current atmospheric pressure

        integer :: isotopeNum ! N_MOLIS ! for categorizing isotopolouges in more broader group
        
        real :: alphaT ! ADD ! The half-width at half-maximum (HWHM) of the Doppler-broadened component

        ! it characterizes the relative contributions of Lorentzian and Doppler broadening to the overall shape of the spectral line
        real(kind=DP) :: shapePrevailFactor ! ALAD ! ratio of is the ratio of the Lorentz HWHM (AL) to the Doppler width (ADD).

        real(kind=DP) :: SL ! temperature-dependent line intensity
        
        ! write(*,*) 'currentLevel: ', currentLevel
        ! write(*,*) 'loopLevel: ', loopLevel
        ! write(*,*) '**************'
        
        if (currentLevel /= loopLevel) then
            currentLevel = loopLevel
            unitlessT = refTemperature/T
            
            ! write(*,*) 'unitlessT: ', unitlessT
            ! pause

            pSelf = RO * 10. / LOSCHMIDT * T/stTemperature
            pForeign = P - pSelf

            ! write(*,*) 'pSelf: ', pSelf
            ! write(*,*) 'pForeign: ', pForeign
            ! pause
        end if

        ! -------- Line-by-line loop (iteration over records in HITRAN file) ------ !

        I = LINBEG - 1 ! MBR !
        ! TODO: change this loop to 'DO WHILE'
        do I = LINBEG, totalLines
            read(7777, rec=I) lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, foreignTempCoeff, &
                                jointMolIso, deltaForeign

            ! write(*,*) 'lineWV: ', lineWV
            ! write(*,*) 'refLineIntensity: ', refLineIntensity
            ! write(*,*) 'gammaForeign: ', gammaForeign
            ! write(*,*) 'gammaSelf: ', gammaSelf
            ! write(*,*) 'lineLowerState: ', lineLowerState
            ! write(*,*) 'foreignTempCoeff: ', foreignTempCoeff
            ! write(*,*) 'jointMolIso: ', jointMolIso
            ! write(*,*) 'deltaForegin: ', deltaForeign
            ! pause

            if  (lineWV >= extEndDeltaWV) exit
            
            if  (lineWV <= capWV) LINBEG = I

            isotopeNum = jointMolIso/100 
            
            ! write(*,*) 'isotopeNum: ', isotopeNum
            ! pause

            gammaPT = LorentzianHWHM(P, T, pSelf, refTemperature, foreignTempCoeff, gammaSelf, gammaForeign)
            ! write(*,*) gammaPT
            ! pause
            shiftedLineWV = pressureShiftedWV(lineWV, deltaForeign, P)
            alphaT = DopplerHWHM(shiftedLineWV, T, WISO(isotopeNum)) ! DIFFERENCE WITH FOMIN on sqrt(log(2)), check in subsequent sections
            ! pause

            ! write(*,*) 'gammaPT (AL): ', gammaPT
            ! write(*,*) 'shiftedLineWV (VI): ', shiftedLineWV
            ! write(*,*) 'alphaT (ADD): ', alphaT
            ! pause
            
            ! calculation of line intensity dependent on temperature based on reference line intensity from HITRAN
            ! write(*,*) 'STR3: ', STR3
            ! write(*,*) 'STS3: ', '***'
            ! write(*,*) 'RO: ', RO
            ! write(*,*) 'PI: ', PI
            lineIntensity = lineIntensityofT(T, refLineIntensity, isotopeNum, QofT, lineLowerState, lineWV)
            ! write(*,*) 'lineIntensity * RO/PI (QDQ): ', lineIntensity * RO / PI
            ! pause

            shapePrevailFactor = gammaPT / alphaT ! <----- ratio to see which effect (Doppler or Lorentz prevails)
            ! write(*,*) 'alphaT (ADD): ', alphaT 
            ! write(*,*) 'shapePrevailFactor (ALAD): ', shapePrevailFactor
            ! pause

            lineIntensity = lineIntensity * RO/PI
            write(*,*) 'lineIntensity before Clough: ', lineIntensity
            ! write(*,*) 'FVVHSL: ', 1/cloughRFunc(T, shiftedLineWV)
            ! pause
            ! lineIntensity = lineIntensity / cloughRFunc(T, shiftedLineWV)
            write(*,*) 'SL with Clough: ', lineIntensity
            pause
            ! SLL = SL * gammaPT --- > it seems that this is needed only in shape functions

            if (shapePrevailFactor > BOUNDL) then
                if (shiftedLineWV < startDeltaWV) then
                    ShapeFuncPtr => Lorentz
                    call LEFTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr) 
                else 
                    if (shiftedLineWV >= endDeltaWV) then
                        ShapeFuncPtr => Lorentz
                        call RIGHTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                    else 
                        ShapeFuncPtr => Lorentz
                        call CENTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                    end if
                end if
            else
                if (shapePrevailFactor > BOUNDD) then
                    if (shiftedLineWV < startDeltaWV) then
                        ShapeFuncPtr => VOIGT
                        call LEFTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                    else
                        if (shiftedLineWV >= endDeltaWV) then
                            ShapeFuncPtr => VOIGT
                            call RIGHTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                        else 
                            ShapeFuncPtr => VOIGT
                            call CENTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                        end if
                    end if
                else 
                    if ( shiftedLineWV < startDeltaWV ) then
                        ShapeFuncPtr => DOPLER
                        call LEFTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                    else 
                        if ( shiftedLineWV >= endDeltaWV) then
                            ShapeFuncPtr => DOPLER
                            call RIGHTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                        else 
                            ShapeFuncPtr => DOPLER
                            call CENTLBL(startDeltaWV, shiftedLineWV, ShapeFuncPtr)
                        end if
                    end if
                end if
            end if 
        end do
        ! ------------------------------------------------------ !

        ! -------- End of line-by-line loop (iteration over records in HITRAN file) --!
    end subroutine LBL2023
end module LBL
