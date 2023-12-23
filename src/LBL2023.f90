module spectroscopy
    use globals
    implicit none
    ! HITRAN spectral data variables
    ! Not 100% sure about the explanations !!
    real(kind=DP) :: V_I_ ! line position in cm^-1
    real :: S_I_  ! line intensity (typically for given temperature, often 296 K)
    
    ! Air-broadened half-width -- the half-width of the spectral line at half-maximum (HWHM)
    ! due to broadening by air (primarily nitrogen and oxygen)
    real :: ALFA_I_ ! HWHM due to air
    real :: ALFAS_I_ ! HWHM due to self-collisions
    
    real :: E_I_  ! Lower-state energy (energy of the lower state of the transition)
    
    ! A factor used to scale the line intensity, often related to isotopic abundance 
    ! or other factors that affect the overall intensity of the line.
    real :: FACT_I_ ! Intensity scaling factor
    
    integer :: NISO_I_ ! Isotopolouge number
    
    ! The shift in the line position due to pressure. 
    ! It represents the change in the central frequency of the line under different pressure conditions.
    real :: SHIFT_I_ ! Pressure-induced line shift

    ! ------- ^^^^ do not change the kinds there ^^^^ ------ !
contains
    real function VOIGT(XXX)
        real(kind=DP), intent(in) :: XXX
    end function VOIGT

    real(kind=DP) function VV_LOR(X)
        real(kind=DP), intent(in) :: X
    end function VV_LOR

    real(kind=DP) function DOPLER(X)
        real(kind=DP), intent(in) :: X
    end function DOPLER

    real(kind=DP) function VAN_VLE(T, V)
        real(kind=DP) :: T, V
    end function VAN_VLE
end module spectroscopy

module LBL
    use globals
    use shared_vars_main
    use MESH1
    use molecule_vars
    use spectroscopy
    use molar_masses, only: WISO
    use LINE_GRID_CALC
    implicit none
contains
    
    subroutine LBL2023(MO_E, LINBEG, VAA, VFISH, NLIN)
        ! consider avoid providing arguments in subroutine.
        ! if they are declared in the external module -- why to declare them twice ?

        ! avoided SAVE blocks
        ! avoided COMMON blocks
        
        ! uncomment when ready to provide line-shape functions as arguments to the 
        ! LEFTLBL, CENTLBL, RIGHTLBL subroutines
        !Ensure all functions assigned to ShapeFuncPtr conform to the shape_func interface.
        ! procedure(shape_func), pointer :: ShapeFuncPtr

        ! ----- move this section below where LEFTLBL, CENTLBL and RIGHTLBL are called ---- !
        ! ! Assign the function pointer to VOIGT and call LEFTLBL
        ! ShapeFuncPtr => VOIGT
        ! call LEFTLBL(..., ShapeFuncPtr)
    
        ! ! Later reassign the function pointer to DOPLER and call RightLBL
        ! ShapeFuncPtr => DOPLER
        ! call RightLBL(..., ShapeFuncPtr)
        ! ----------------------------------------------------------------------------------!

        integer :: MO_E ! molecule type-integer: '1' - for SO2 and H2O, '2' -- for CO2
        integer :: LINBEG ! integer line label used for locating record in direct access file
        integer :: I ! loop variable for accessing the record in direct access file
        integer :: NLIN ! number of lines in one file
        ! consider removing it from the subroutine argument list, otherwise warning 
        ! that this variable is already defined in the parent scope

        real(kind=DP) :: VAA, VFISH

        real :: TOLD, ISOOLD ! likely old values for save block
        real(kind=DP) :: VS, VF ! MBR; similar to VSTART, VS, VFINISH

        real, parameter :: PI = 3.141593 ! pi
        real, parameter :: TS = 296. ! Standard temperature
        real, parameter :: PS = 1. ! Standard pressure in atmsopheres

        real, parameter :: BOUNDL = 10. ! some boundary parameter
        real, parameter :: BOUNDD = 0.01 ! likely some boundary value related to Doppler broadening, given its small value
        real, parameter :: BET = 1.438786  ! hc/k -- possibly second radiation constant in specific units (look for Planck law)

        real(kind=DP) :: DLT8 ! likely the wavenumber-like delta with DP for loop 

        real :: T1, P1, RO1 ! replicas of T, R, P -- can't find in legacy where they are initialised !
        
        real :: TT ! temperature
        real :: TST ! unitless temperature (div by standard temperature)
        
        real :: EK ! some temperature dependent factor for energy-related calculations
        ! ---------- second radiation constant is used for EK calculations (BET)
        
        ! variables likely related to temperature interpolation or indexing within
        ! temperature-realted data  (e.g. partition functions)
        integer :: NTAB_G ! likely an index or pointer
        real :: T_G1 ! likely temperature index corresponding to the specific index
        real :: C_G1, C_G2 ! likely coefficients used for linear interpolation between table entries
        
        ! variables related to pressure calculations
        real :: PSELF ! self-broadening pressure component
        real :: PFOREI ! foreign gas broadening component
        ! ^^^ likely calculated from the total pressure P, density RO and standard conditions 
        real :: APF ! likely ratio of foregin gas pressure to standard pressure
        real :: APF2 ! likely ratio of foregin gas pressure to total pressure
        real :: APS ! ratio of self-broadening pressure to standard pressure
        real :: APS2 ! ration of self-broadening pressure to total pressure
        
        real(kind=DP) :: SLSS ! replica of S_I_ (line intensity)
        real(kind=DP) :: APALS, APALF ! scaled self- and foreign- broadened HWMW
        real(kind=DP) :: AL ! Lorentz HWHM -- T and P dependent !
        
        ! Appears in the denominator of the Lorentzian line profile L(\nu)
        real(kind=DP) :: ALAL ! squared Lorentz half-width
        real(kind=DP) :: VI ! likely shifted line position under the current atmospheric pressure

        integer :: ISO ! replica of NISO_I_ isotopolouge number
        
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
        
        ! zero initialisation of pressure, temperature and density from shared_vars_main module
        ! MBR

        ! CHANGE THIS IN LEGACY CODE !!! T, P, RO are fetched from the module shared_vars_main !
        ! T = 0. 
        ! P = 0.
        ! RO = 0.  

        TOLD = -5.
        ISOOLD = -5.

        VS = 1.0-04
        VF = 1.0-04

        MOTYPE = MO_E ! might be confusing with `select case` block from the K_HITRAN subroutine

        if ( VS /= VSTART ) then 
            DLT8 = 10.0
            VS = VSTART
            VF = VS + DLT8
        end if

        ! DEBUG SECTION
        ! write(*,*) 'T in LBL2023.f90: ', T
        ! pause

        ! if ( P /= P1 .OR. T/=T1 .OR. RO /= RO1 ) then
        !     if ( T /= T1 ) then
                !T = T1 
                TT = T
                TST = TS/T
                EK = BET * (T - TS)/(T * TS)
                NTAB_G = (T - 20.0)/2 + 1
                T_G1 = NTAB_G * 2.0 + 18.
                C_G2 = (T - T_G1)/2.
                C_G1 = 1. - C_G2
                
                !###  296K -> N_TAB=139 ###TIPS=C_G1*QofT(N_MOLIS,NTAB_G)+C_G2*QofT(N_MOLIS,NTAB_G+1) ! <<--- WTF is it ??
                
                ! DEBUG SECTION
                ! write(*,*) 'T: ', T
                ! write(*,*) 'NTAB_G', NTAB_G 
                ! pause
                
            ! end if
            ! P = P1
            ! RO = RO1
            PSELF = RO * 10. / 2.6872E25 * T/273.15
            PFOREI = P - PSELF
            APF = PFOREI/PS
            APS = PSELF/PS
            APF2 = PFOREI/P
            APS2 = PSELF/P 
        ! end if

        ! -------- Line-by-line loop (iteration over records in HITRAN file) ------ !

        I = LINBEG - 1 ! MBR !
        do I = LINBEG, NLIN

            !DEBUG SECTION
            !write(*,*) '! ---------- Reading and preprpocessing the HITRAN data ---------------------- !'
            !write(*,*) 'I= ', I
            !write(*,*) 'NLIN= ', NLIN
            !pause
            if (I<=0) then
                !!!!! 10th level falls !!!
                write(*,*) '!!!!!!!', I, T;
                exit
            end if
            read(7777, rec=I) V_I_, S_I_, ALFA_I_, ALFAS_I_, E_I_, FACT_I_, NISO_I_, SHIFT_I_

            ! DEBUG SECTION
            ! write(*,*) '!------ read HITRAN data ------ !'
            ! write(*,*) 'V_I_: ', V_I_
            ! write(*,*) 'S_I_: ', S_I_
            ! write(*,*) 'ALFA_I: ', ALFA_I_
            ! write(*,*) 'ALFAS_I_: ', ALFAS_I_
            ! write(*,*) 'E_I_: ', E_I_
            ! write(*,*) 'FACT_I_: ', FACT_I_
            ! write(*,*) 'NISO_I_: ', NISO_I
            ! write(*,*) 'SHIFT_I: ', SHIFT_I_
            ! pause
            
            if (V_I_ >= VFISH) exit

            if (V_I_ <= VAA) then
                LINBEG = I ! MBR !
            end if 

            SLSS = S_I_
            APALF = APF * ALFA_I_
            APALS = APS * ALFAS_I_

            ! Lorentz half-width
            AL = TST ** FACT_I_ * (APALF * APF2 + APALS * APS2) !  <---- FORMULA !!

            ! squared Lorentz half-width
            ALAL = AL * AL

            VI = V_I_ + SHIFT_I_ * P ! <------- FORMULA !!
            
            ISO = NISO_I_
            N_MOLIS = NISO_I_/100 

            !!! WHY ?
            if (ISO == 222) ISO = 221 ! Attention: the second N2 isotopolouge is treated as the first !!!
            
            if (ISO /= ISOOLD .OR. T /= TOLD) then
                ISOOLD = ISO
                TOLD = T ! <--- CHANGED LEGACY CODE !!! likely a tipo (was TOOLD)
                
                DOPCON = 4.29e-7 * sqrt(T/WISO(N_MOLIS)) ! <--------- FORMULA !!

                !* Total internal partition function -- statistics sum ? *!
                
                !DEBUG SECTION
                ! write(*,*) 'QofT dimensions:', size(QofT, 1), size(QofT, 2)
                ! write(*,*) 'N_MOLIS: ', N_MOLIS
                ! write(*,*) 'NTAB_G', NTAB_G 
                ! pause

                STS3 = C_G1 * QofT(N_MOLIS, NTAB_G) + C_G2 * QofT(N_MOLIS, NTAB_G+1) ! <---- INTERPOLATION FORMULA between two values of the partition function
                STR3 = QofT(N_MOLIS, 139) ! At 296K
                
                ! ration of the partition function at the standard temperature to the interpolated
                QDQ = STR3/STS3 * RO/PI ! <-------- FORMULA (scaled partition function)
            end if

            ! Line 166 from the LEGACY CODE
            ADD = DOPCON * VI ! <----- FORMULA

            ALAD = AL / ADD ! <----- ratio

            SL = SLSS * QDQ

            if (STR3 >= 0.) SL = SL * exp(EK * E_I_)  ! <-- apply correction if st sum > 0

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
            BETVITS = BET * VI / TS
            BETVIT = BET * VI / T

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
                if (VI < VS) then
                    call LEFTLBL(VS, VI, VV_LOR, EPS) 
                else 
                    if (VI >= VF) then
                        call RIGHTLBL(VS, VI, VV_LOR, EPS)
                    else 
                        call CENTLBL(VS, VI, VV_LOR, EPS)
                    end if
                end if
            else
                if (ALAD > BOUNDD) then
                    if (VI < VS) then
                        call LEFTLBL(VS, VI, VOIGT, EPS)
                    else
                        if (VI >= VF) then
                            call RIGHTLBL(VS, VI, VOIGT, EPS)
                        else 
                            call CENTLBL(VS, VI, VOIGT, EPS)
                        end if
                    end if
                else 
                    if ( VI < VS ) then
                        call LEFTLBL(VS, VI, DOPLER, EPS)
                    else 
                        if ( VI >= VF) then
                            call RIGHTLBL(VS, VI, DOPLER, EPS)
                        else 
                            call CENTLBL(VS, VI, DOPLER, EPS)
                        end if
                    end if
                end if
            end if 
        end do

        ! ------------------------------------------------------ !

        ! -------- End of line-by-line loop (iteration over records in HITRAN file) --!
    end subroutine LBL2023
end module LBL
