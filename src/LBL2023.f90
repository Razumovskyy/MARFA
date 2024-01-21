module LBL
    use kinds
    use constants, only: PI, refTemperature, stPressure, stTemperature, LOSCHMIDT
    use shared_vars_main
    use mesh
    use molecule_vars
    use spectroscopy
    use molar_masses, only: WISO
    use shape_functions
    use LINE_GRID_CALC
    implicit none
contains

    !subroutine LBL2023(MO_E, LINBEG, VAA, VFISH, totalLines)
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

        real :: TOLD, ISOOLD ! likely old values for save block

        real, parameter :: BOUNDL = 10. ! some boundary parameter
        real, parameter :: BOUNDD = 0.01 ! likely some boundary value related to Doppler broadening, given its small value

        real, save :: unitlessT ! unitless temperature (div by standard temperature)
        
        real, save :: EK ! some temperature dependent factor for energy-related calculations
        ! ---------- second radiation constant is used for EK calculations (C2)
        
        ! variables likely related to temperature interpolation or indexing within
        ! temperature-realted data  (e.g. partition functions)
        integer, save :: NTAB_G ! likely an index or pointer
        real, save :: T_G1 ! likely temperature index corresponding to the specific index
        real, save :: C_G1, C_G2 ! likely coefficients used for linear interpolation between table entries
        
        ! variables related to pressure calculations
        real, save :: pSelf ! self-broadening pressure component
        real, save :: pForeign ! foreign gas broadening component
        
        real(kind=DP) :: SLSS ! replica of lineIntensity (line intensity)
        real(kind=DP) :: APALS, APALF ! scaled self- and foreign- broadened HWMW
        real(kind=DP) :: AL ! Lorentz HWHM -- T and P dependent !
        
        ! Appears in the denominator of the Lorentzian line profile L(\nu)
        real(kind=DP) :: ALAL ! squared Lorentz half-width
        real(kind=DP) :: VI ! likely shifted line position under the current atmospheric pressure

        integer :: ISO ! replica of jointMolIso isotopolouge number
        
        !  The division by 100 could be a way to obtain a molecular identifier 
        ! from the isotopologue number. For example, if isotopologues are numbered such
        ! that the first two digits represent the molecule and the remaining digits 
        ! represent different isotopes of that molecule, 
        ! then this division would extract the molecular part of the identifier.
        ! N_MOLIS treated as an index to WISO array of molar masses
        integer :: N_MOLIS ! for categorizing isotopolouges in more broader group
        
        real :: DOPCON ! likely the Doppler broadening constant

        real(kind=DP) :: STS3 ! likely temperature-dependent partition function for a specific isotopologue
        real(kind=DP) :: STR3 ! likely the partition function value at the standard temperature
        real(kind=DP) :: QDQ ! scaled partition function 

        real(kind=DP) :: ADD ! Doppler width
        
        ! it characterizes the relative contributions of Lorentzian and Doppler broadening to the overall shape of the spectral line
        real(kind=DP) :: ALAD ! ratio of is the ratio of the Lorentz HWHM (AL) to the Doppler width (ADD).

        real(kind=DP) :: SL ! temperature-dependent line intensity
        
        real(kind=DP) :: BETVITS, BETVIT ! scaled version of the line position
        
        real(kind=DP) :: EXPVV, EXPVVS ! exponential terms used in the calculation of Van Vleck-Weiss-Huber factor
        real(kind=DP) :: FVVHSL ! VV-W-H factor (calculation in two ways: for small BETVITS value approximation is used)
        real(kind=DP) :: SLL ! can be interpreted as a measure of the contribution of the Lorentzian component to the overall line profile.


        ! CHANGE THIS IN LEGACY CODE !!! T, P, RO are fetched from the module shared_vars_main !
        ! T = 0. 
        ! P = 0.
        ! RO = 0.
        
        ! TOLD = -5.
        ! ISOOLD = -5.

        ! VS = 1.0-04
        ! endDeltaWV = 1.0-04

        ! if ( VS /= startDeltaWV ) then 
        !     VS = startDeltaWV
        !     endDeltaWV = VS + deltaWV
        ! end if

        ! DEBUG SECTION
        ! write(*,*) 'T in LBL2023.f90: ', T
        ! pause

        if (currentLevel /= loopLevel) then
            currentLevel = loopLevel
            unitlessT = refTemperature/T
            EK = C2 * (T - refTemperature)/(T * refTemperature)
            NTAB_G = (T - 20.0)/2 + 1
            T_G1 = NTAB_G * 2.0 + 18.
            C_G2 = (T - T_G1)/2.
            C_G1 = 1. - C_G2
            
            !###  296K -> N_TAB=139 ###TIPS=C_G1*QofT(N_MOLIS,NTAB_G)+C_G2*QofT(N_MOLIS,NTAB_G+1) ! <<--- WTF is it ??
            
            ! DEBUG SECTION
            ! write(*,*) 'T: ', T
            ! write(*,*) 'NTAB_G', NTAB_G 
            ! pause

            pSelf = RO * 10. / LOSCHMIDT * T/stTemperature
            ! pSelf = RO * 10 * BOL * T ! more simpler ! 10 factor is because concentration data is in molecules/(cm^2*km)
            pForeign = P - pSelf
        end if

        ! -------- Line-by-line loop (iteration over records in HITRAN file) ------ !

        I = LINBEG - 1 ! MBR !
        do I = LINBEG, totalLines

            !DEBUG SECTION
            !write(*,*) '! ---------- Reading and preprpocessing the HITRAN data ---------------------- !'
            !write(*,*) 'I= ', I
            !write(*,*) 'totalLines= ', totalLines
            !pause

            read(7777, rec=I) lineWV, lineIntensity, gammaForeign, gammaSelf, lineLowerState, foreignTempCoeff, &
                                jointMolIso, deltaForeign

            ! DEBUG SECTION
            ! write(*,*) '!------ read HITRAN data ------ !'
            ! write(*,*)  lineWV: ', lineWV
            ! write(*,*) 'lineIntensity: ', lineIntensity
            ! write(*,*) 'ALFA_I: ', gammaForeign
            ! write(*,*) 'gammaSelf: ', gammaSelf
            ! write(*,*) 'lineLowerState: ', lineLowerState
            ! write(*,*) 'foreignTempCoeff: ', foreignTempCoeff
            ! write(*,*) 'jointMolIso: ', NISO_I
            ! write(*,*) 'SHIFT_I: ', deltaForeign
            ! pause
            
            if  (lineWV >= extEndDeltaWV) exit
            
            if  (lineWV <= capWV) LINBEG = I

            SLSS = lineIntensity
            APALF = (pForeign / stPressure) * gammaForeign
            APALS = (pSelf / stPressure) * gammaSelf

            ! Lorentz half-width
            AL = unitlessT ** foreignTempCoeff * (APALF * (pForeign/P) + APALS * (pSelf/P))

            ! squared Lorentz half-width
            ALAL = AL * AL

            VI = lineWV + deltaForeign * P ! <------- FORMULA !!
            
            ISO = jointMolIso
            N_MOLIS = jointMolIso/100 

            !!! WHY ?
            if (ISO == 222) ISO = 221 ! Attention: the second N2 isotopolouge is treated as the first !!!
            
            if (ISO /= ISOOLD .OR. T /= TOLD) then
                ISOOLD = ISO
                TOLD = T ! <--- CHANGED LEGACY CODE !!! likely a tipo (was TOOLD)
                
                DOPCON = 4.29e-7 * sqrt(T/WISO(N_MOLIS)) ! <--------- FORMULA !!

                !* Total internal partition function -- statistics sum ? *!
                
                !DEBUG SECTION
                ! write(*,*) 'QofT dimensions:', size(QofT, 1), size(QofT, 2)
                ! pause
                ! ! write(*,*) 'N_MOLIS: ', N_MOLIS
                ! write(*,*) 'NTAB_G', NTAB_G

                STS3 = C_G1 * QofT(N_MOLIS, NTAB_G) + C_G2 * QofT(N_MOLIS, NTAB_G+1) ! <---- INTERPOLATION FORMULA between two values of the partition function
                STR3 = QofT(N_MOLIS, 139) ! At 296K
                
                ! ration of the partition function at the standard temperature to the interpolated
                QDQ = STR3/STS3 * RO/PI ! <-------- FORMULA (scaled partition function)
            end if

            ! Line 166 from the LEGACY CODE
            ADD = DOPCON * VI ! <----- FORMULA

            ALAD = AL / ADD ! <----- ratio

            SL = SLSS * QDQ

            if (STR3 >= 0.) SL = SL * exp(EK * lineLowerState)  ! <-- apply correction if st sum > 0

            ! ----- LEGACY CODE COMMENTARY --- V-W-H factor for LTE and non-LTE ---------- !
            ! LTE:	Pure radiation  (van Vleck-Weisskopf-Huber) factor!!! 26 Feb.,2009  !!!	*
            ! b=C2=1.438786 the second radiation constant 
            !  Fi(V)=S*[1/FACTOR(VI)]*FACTOR(V)] =
            !  Sref*(...)*[(1+exp(-Vi*b/T)/[Vi*(1-exp(-Vi*b/Tref))] x       ! see below  
            ! *** ATTENTON=> there are T in the NUMENATOR    and  Tref in the DENUMERATOR *** LTE-case !!!     
            !       x [(1-exp(-V*b/T)/(1+exp(-V*b/T))]*V        ! see K_COEFF
            ! <<< in this non-LTE program only T is used (not referenced !!!)
            ! ----------------------------------------------------------------------------- !

            ! This expressions calculate scaled versions of the line position at the standard temperature and current temperature
            ! may be used in temperature-dependent ajustments
            BETVITS = C2 * VI / refTemperature
            BETVIT = C2 * VI / T

            if (BETVITS > 1E-5) then
                EXPVV = exp(-BETVIT)
                EXPVVS = exp(-BETVITS)
                FVVHSL = (1. + EXPVV) / (1. - EXPVVS) / VI
            else 
                FVVHSL = (2. + BETVIT) / BETVITS / VI
            end if

            SL = FVVHSL * SL ! applying VV-W-H factor on line intensity
            SLL = SL * AL

            if (ALAD > BOUNDL) then
                if (VI < startDeltaWV) then
                    ShapeFuncPtr => Lorentz
                    call LEFTLBL(startDeltaWV, VI, ShapeFuncPtr) 
                else 
                    if (VI >= endDeltaWV) then
                        ShapeFuncPtr => Lorentz
                        call RIGHTLBL(startDeltaWV, VI, ShapeFuncPtr)
                    else 
                        ShapeFuncPtr => Lorentz
                        call CENTLBL(startDeltaWV, VI, ShapeFuncPtr)
                    end if
                end if
            else
                if (ALAD > BOUNDD) then
                    if (VI < startDeltaWV) then
                        ShapeFuncPtr => VOIGT
                        call LEFTLBL(startDeltaWV, VI, ShapeFuncPtr)
                    else
                        if (VI >= endDeltaWV) then
                            ShapeFuncPtr => VOIGT
                            call RIGHTLBL(startDeltaWV, VI, ShapeFuncPtr)
                        else 
                            ShapeFuncPtr => VOIGT
                            call CENTLBL(startDeltaWV, VI, ShapeFuncPtr)
                        end if
                    end if
                else 
                    if ( VI < startDeltaWV ) then
                        ShapeFuncPtr => DOPLER
                        call LEFTLBL(startDeltaWV, VI, ShapeFuncPtr)
                    else 
                        if ( VI >= endDeltaWV) then
                            ShapeFuncPtr => DOPLER
                            call RIGHTLBL(startDeltaWV, VI, ShapeFuncPtr)
                        else 
                            ShapeFuncPtr => DOPLER
                            call CENTLBL(startDeltaWV, VI, ShapeFuncPtr)
                        end if
                    end if
                end if
            end if 
        end do
        ! ------------------------------------------------------ !

        ! -------- End of line-by-line loop (iteration over records in HITRAN file) --!
    end subroutine LBL2023
end module LBL
