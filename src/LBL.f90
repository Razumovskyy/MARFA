module LBL
    use Constants
    use ShapeFuncInterface
    use Mesh
    use IO
    use Spectroscopy
    use MolarMasses, only: WISO
    use Shapes
    use LineGridCalc
    implicit none
contains

    subroutine modernLBL(molTypeArg, LINBEG, capWV, totalLines, loopLevel)
        
        integer :: molTypeArg
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

        real :: lineIntensity ! temperature and pressure dependent line intensity
        
        real, save :: pSelf ! self-broadening pressure component
        real, save :: pForeign ! foreign gas broadening component
        
        ! Appears in the denominator of the Lorentzian line profile L(\nu)
        real(kind=DP) :: shiftedLineWV ! VI ! shifted line position under the current atmospheric pressure

        integer :: isotopeNum ! N_MOLIS ! for categorizing isotopolouges in more broader group
        
        real :: alphaT ! ADD ! The half-width at half-maximum (HWHM) of the Doppler-broadened component

        ! it characterizes the relative contributions of Lorentzian and Doppler broadening to the overall shape of the spectral line
        real(kind=DP) :: shapePrevailFactor ! ALAD ! ratio of is the ratio of the Lorentz HWHM (AL) to the Doppler width (ADD).

        ! -------- Line-by-line loop (iteration over records in HITRAN file) ------ !

        do I = LINBEG, totalLines
            read(7777, rec=I) lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, foreignTempCoeff, &
                                jointMolIso, deltaForeign

            if  (lineWV >= extEndDeltaWV) exit
            
            if  (lineWV <= capWV) LINBEG = I

            isotopeNum = jointMolIso/100

            !gammaPT = LorentzianHWHM(P, T, pSelf, refTemperature, foreignTempCoeff, gammaSelf, gammaForeign)
            
            !shiftedLineWV = pressureShiftedWV(lineWV, deltaForeign, P)

            !alphaT = DopplerHWHM(shiftedLineWV, T, WISO(isotopeNum))
            
            !lineIntensity = lineIntensityofT(T, refLineIntensity, isotopeNum, QofT, lineLowerState, lineWV)

            !shapePrevailFactor = gammaPT / alphaT ! <----- ratio to see which effect (Doppler or Lorentz prevails)

            ! SLL = gammaPT * 

            if (shapePrevailFactor > BOUNDL) then
                if (shiftedLineWV < startDeltaWV) then
                    ! shapeFuncPtr => lorentz
                    ! call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr) 
                else 
                    if (shiftedLineWV >= endDeltaWV) then
                        ! shapeFuncPtr => lorentz
                        ! call rightLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                    else 
                        ! shapeFuncPtr => lorentz
                        ! call centerLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                    end if
                end if
            else
                if (shapePrevailFactor > BOUNDD) then
                    if (shiftedLineWV < startDeltaWV) then
                        ! shapeFuncPtr => voigt
                        ! call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                    else
                        if (shiftedLineWV >= endDeltaWV) then
                            ! shapeFuncPtr => voigt
                            ! call rightLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                        else 
                            ! shapeFuncPtr => voigt
                            ! call centerLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                        end if
                    end if
                else 
                    if ( shiftedLineWV < startDeltaWV ) then
                        ! shapeFuncPtr => doppler
                        ! call leftLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                    else 
                        if ( shiftedLineWV >= endDeltaWV) then
                            ! shapeFuncPtr => doppler
                            ! call rightLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                        else 
                            ! shapeFuncPtr => doppler
                            ! call centerLBL(startDeltaWV, shiftedLineWV, shapeFuncPtr)
                        end if
                    end if
                end if
            end if 
        end do
        ! -------- End of line-by-line loop (iteration over records in HITRAN file) --!
    end subroutine modernLBL
end module LBL
