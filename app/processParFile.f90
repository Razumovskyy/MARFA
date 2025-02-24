! -------------------------------------------------------------------------------------------------
! Each item is defined below, with its format shown in parenthesis.
! To run marfa: 9 parameters from .par files are needed ONLY:
! -------------------------------------------------------------------------------------------------
!      1          (I2)          MO = molecule number {I use MO_ISO= MO*10+ISO (I3)}
!      2          (A1)          ISO = isotope number {1 = most abundant, 2 = second, ..., 0 = 10, A = 11, ...}
!      3          (F12.6)       lineWV (V) = frequency of transition in wavenumbers (cm-1)
!      4          (E10.3)       refLineIntentsity (S) = intensity in cm-1/(molec * cm-2) at 296 Kelvin
!*  not needed    (E10.3)       A = Einstein A coefficient
!      5          (F5.4)        gammaForeign (AGAM) = air-broadened halfwidth  in cm-1/atm at 296 Kelvin
!      6          (F5.3)        gammasSelf (SGAM) = self-broadened halfwidth in cm-1/atm at 296 Kelvin
!      7          (F10.4)       lineLowerState (E) = lower state energy in wavenumbers (cm-1)
!      8          (F4.2)        foreignTempCoeff (N) = coeff. of temperat. dependence of air-broadened halfwidth
!      9          (F8.6)        deltaForeign (d) = shift of transition due to pressure (cm-1)
!*  not needed    (A15)         V1 = upper state global quanta 
!*  not needed    (A15)         V2 = lower state global quanta 
!*  not needed    (A15)         Q1 = upper state local quanta
!*  not needed    (A15)         Q2 = lower state local quanta
!*  not needed    (6I1)         IERR = accuracy indexes
!*  not needed    (6I2)         IREF = reference indexes
!*  not needed    (A1)          FLAG = flag for line coupling {(*) if "Yes"}
!*  not needed    (F7.1)        Gu = upper stat. weight
!*  not needed    (F7.1)        Gl = lower stat. weight
! -------------------------------------------------------------------------------------------------
!
! The molecule numbers are encoded as shown in the table below (https://hitran.org/docs/molec-meta/):
!
! -------------------------------------------------------------------------------------------------
!  1=   H2O   2=    CO2   3=   O3    4=   N2O   5=   CO   6=     CH4
!  7=   O2    8=    NO    9=   SO2   10=  NO2   11=  NH3  12=    HNO3   
! -------------------------------------------------------------------------------------------------
program processParFile
    implicit none
    ! selected kind for double precision
    integer, parameter :: DP = selected_real_kind(15,  307)
    
    ! spectral parameters which are needed:
    integer :: MO ! Molcelule ID
    character(len=1) ISO ! Isotopologue ID
    real(kind=DP) lineWV ! [cm-1] -- current spectral line, transition wavenumber
    real(kind=DP) refLineIntensityDP ! lineIntensity initialy read in double precision, than converted to single precision
    real :: refLineIntensity ! [cm-1/(molecule*cm-2)] -- spectral line intensity at refTemperature=296 K
    real :: gammaForeign ! [cm-1/atm] -- Lorentzian foreign-broadened Lorentz HWHM at refTemperature=296 K
    real :: gammaSelf ! [cm-1/atm] -- self-broadened component of Lorentz HWHM
    real(kind=DP) lineLowerStateDP ! lowerStateEnergy initialy read in double precision, than converted to single precision
    real :: lineLowerState ! [cm-1] -- lower state energy of the transition
    real :: foreignTempCoeff ! N ! [dimensionless] (coefficient for temperature dependence of gammaForeign)
    real :: deltaForeign ! [cm-1/atm] (pressure shift of the line position at 296 K and 1 atm)
    
    ! spectral parameters which are not needed
    real :: A ! Einstein A coefficient
    character(len=15) V1, V2 ! global upper and lower quanta
    character(len=15) Q1, Q2 ! local upper and lower quanta
    integer :: IERR(6), IREF(6) ! accuracy and reference indexes
    character(len=1) FLAG ! line mixing flag
    real(kind=DP) Gu, Gl ! The lower and upper state statistical weights

    ! command line arguments
    character(len=10) :: inputMolecule
    character(len=20) :: databaseSlug ! Slug for the output database file
    character(len=100) :: rawSpectralFile ! path to the file with raw Spectral data
    
    ! other parameters
    integer :: IS ! Isotopologue ID transformed to numeric format
    integer :: jointMolIso !  [dimensionless] custom variable: joined reference to Molecule number (MOL) and Isotopologue number (ISO)
    integer :: isotopeNumber
    character(len=2) nameSuffix
    integer :: i
    integer, parameter :: N_MOL_ISO(20,12) = reshape([ &
        ! 1 H2O
        1, 2, 3, 4, 5, 6, 7, 8, 9, [(0, i=1,11)], &
        ! 2 CO2
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, [(0, i=1,7)], &
        ! 3 O3
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, [(0, i=1,2)], &
        ! 4 N2O
        41, 42, 43, 44, 45, [(0, i=1,15)], &
        ! 5 CO
        46, 47, 48, 49, 50, 51, 52, 53, 54, [(0, i=1,11)], &
        ! 6 CH4
        55, 56, 57, 58, [(0, i=1,16)], &
        ! 7 O2
        59, 60, 61, 62, 63, 64, [(0, i=1,14)], &
        ! 8 NO
        65, 66, 67, [(0, i=1,17)], &
        ! 9 SO2
        68, 69, [(0, i=1,18)], &
        ! 10 NO2
        70, [(0, i=1,19)], &
        ! 11 NH3
        71, 72, [(0, i=1,18)], &
        ! 12 HNO3
        73, 74, [(0, i=1,18)] ], shape=[20,12])
    character(len=81) :: parFormat = '(I2, A1, F12.6, 2E10.3, 2F5.4, F10.4, F4.2, F8.6, 4A15, 6I1, 6I2, A1, F7.1, F7.2)'

    ! technical variables and parameters
    integer :: record_length
    integer :: argc
    integer, parameter :: rawSpectralFileUnit = 331
    integer, parameter :: outputFileUnit = 341
    integer, parameter :: outputHumanReadableFileUnit = 744
    integer :: l
    integer :: ios
    integer :: lineCounter
! -------------------------------------------------------------------------------------------------

    argc = command_argument_count()

    ! Check if the number of arguments is not equal to 8
    if (argc /= 3) then
        print *, 'Incorrect number of arguments.'
        print *, 'Expected 3 arguments, but received ', argc
        print *, 'Usage: program_name databaseSlug rawSpectralFile'
        stop 1
    end if

    do l = 1, command_argument_count()
        select case (l)
        case(1)
            call get_command_argument(l, inputMolecule)  
        case (2)
            call get_command_argument(l, databaseSlug)
        case (3)
            call get_command_argument(l, rawSpectralFile)
        end select
    end do

    nameSuffix = getSpeciesCode(inputMolecule)

    if (nameSuffix == "00") then
        print *, 'ERROR: unsupported molecule: ', trim(adjustl(inputMolecule))
        stop 2
    end if

    ! open .par file
    open(rawSpectralFileUnit, file=rawSpectralFile)

    ! set length of the record for output unformatted file where necessary spectral data would be stored
    ! and open this file:
    record_length = 36 ! 4 * 7 bytes + 8 bytes (for DP lineWV)
    open(outputFileUnit, access='DIRECT', form='UNFORMATTED', recl=record_length, &
            file="data"//"/"//"databases"//"/"//trim(adjustl(databaseSlug))//"."//nameSuffix)

    ! open output file where data would be stored in human-readable format
    open(outputHumanReadableFileUnit, file=trim(adjustl(databaseSlug))//'.'//nameSuffix//'.dat')
    
    lineCounter = 0
    ios = 0
    do while (.not. is_iostat_end(ios))
        lineCounter = lineCounter + 1
        read(unit=rawSpectralFileUnit, fmt=parFormat, iostat=ios) MO, ISO, lineWV, refLineIntensityDP, A, &
                            gammaForeign, gammaSelf, lineLowerStateDP, &
                            foreignTempCoeff, deltaForeign, V1, V2, Q1, Q2, IERR, IREF, FLAG, Gu, Gl
        
        if (ios > 0) then
            print *, 'ERROR: when reading file:', rawSpectralFile
            stop 3
        end if

        IS = convert_ISO(ISO)

        if (IS == 99) then
            print *, 'WARNING: unknown ISO value: ', ISO, 'skipping ...'
            cycle
        end if
        
        ! setting ordinary precision for lineIntensity and lower state energy
        refLineIntensity = refLineIntensityDP
        lineLowerState = lineLowerStateDP

        isotopeNumber = N_MOL_ISO(IS,MO) ! The isotop number in the internal database 
                                         ! (to define mol. weight and TIPS-17).

        jointMolIso = isotopeNumber*100 + MO    ! [isotop number*100 + molecular number]

        write(outputFileUnit, rec=lineCounter)  lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, &
                                                 foreignTempCoeff, jointMolIso, deltaForeign
        write(outputHumanReadableFileUnit, *)  lineWV, refLineIntensity, gammaForeign, gammaSelf, lineLowerState, & 
                                                 foreignTempCoeff, jointMolIso, deltaForeign
    end do
    
    close(outputFileUnit)
    close(outputHumanReadableFileUnit)
    print *, 'The file preprocession is finished !'
    print *, 'number of lines: ', lineCounter-1

contains

    function getSpeciesCode(species) result(code)
        character(len=*), intent(in) :: species
        character(len=2) :: code
   
        ! initialize code to 0 (not found)
        code = "00"
    
        select case (trim(adjustl(species)))
            case ("H2O")
                code = "01"
            case ("CO2")
                code = "02"
            case ("O3")
                code = "03"
            case ("N2O")
                code = "04"
            case ("CO")
                code = "05"
            case ("CH4")
                code = "06"
            case ("O2")
                code = "07"
            case ("NO")
                code = "08"
            case ("SO2")
                code = "09"
            case ("NO2")
                code = "10"
            case ("NH3")
                code = "11"
            case ("HNO3")
                code = "12"
            case default
                code = "00"  ! species not found
        end select
    end function getSpeciesCode

    function convert_ISO(ISOparameter) result(converted_ISO)
        character(len=1), intent(in) :: ISOparameter
        integer :: converted_ISO
        select case (ISOparameter)
            case ('1') 
                converted_ISO = 1
            case ('2') 
                converted_ISO = 2
            case ('3') 
                converted_ISO = 3
            case ('4') 
                converted_ISO = 4
            case ('5') 
                converted_ISO = 5
            case ('6') 
                converted_ISO = 6
            case ('7') 
                converted_ISO = 7
            case ('8') 
                converted_ISO = 8
            case ('9') 
                converted_ISO = 9
            case ('0') 
                converted_ISO = 10
            case ('A') 
                converted_ISO = 11
            case ('B') 
                converted_ISO = 12
            case ('C') 
                converted_ISO = 13
            case ('D') 
                converted_ISO = 14
            case ('E') 
                converted_ISO = 15
            case ('F') 
                converted_ISO = 16
            case ('G') 
                converted_ISO = 17
            case ('H') 
                converted_ISO = 18
            case ('I') 
                converted_ISO = 19
            case default
                converted_ISO = 99
            end select
    end function convert_ISO
end program processParFile
