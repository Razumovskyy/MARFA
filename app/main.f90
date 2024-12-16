program main
    use Constants
    use Atmosphere
    use Spectroscopy
    use LBL
    use Grids
    implicit none

    !------------------------------------------------------------------------------------------------------------------!
    ! Output files units !
    integer, parameter :: outputUnit = 47 ! output file unit where PT-tables are stored
    integer, parameter :: infoUnit = 67
    integer, parameter :: latestRunUnit = 87
    
    ! Record number in the output file
    integer :: outputRecNum 

    ! Input command line arguments and their trimmed values as strings !
    character(len=20) :: databaseSlug
    character(len=10) :: startWVcla, endWVcla, startWVclaTrimmed, endWVclaTrimmed
    character(len=10) :: cutOffcla, cutOffclaTrimmed
    character(len=30) :: targetValue, targetValuecla
    character(len=10) :: inputMolecule   
    
    ! Parameters and variables for constructing directories names !
    character(len=200), parameter :: rootDirName = 'output'
    character(len=300) :: subDirName
    character(len=100) :: formattedStartWV, formattedEndWV
    character(len=200) :: parentDir
    character(len=300) :: fullSubDirPath
    character(len=500) :: mkdirCommand
    character(len=300) :: infoFilePath
    character(len=300) :: latestRunFilePath
    
    ! Technical variables !
    integer :: argc ! number of the passed command line arguments
    integer :: dateTimeValues(8)
    integer :: year, month, day, hour, minute, second
    real(kind=DP) startTime, endTime ! times for measuring a machine time
    character(len=20) :: timestamp ! 
    integer :: status
    character(len=3) :: levelLabel ! unique identifier for each atmospheric level: ('1__','2__',...,'100',...)
    integer :: l ! loop variable

    ! database file extension based on the species title (input molecule), see `data/databases/` directory
    character(len=2) :: DBfileExtension
    
    ! after the atmospheroc levels loop steps on the next level, LBL starts from startingLineIdx
    ! Index (record number in database file) from which the whole calculation goes
    ! Determined by startWV and cutOff
    integer :: startingLineIdx
    
    ! Transition wavenumber of a spectral line at index `lineBegIdx ! VAA0
    ! Mostly it is used to determine capWV
    real(kind=DP) :: startingLineWV
    
    integer :: lineIdx ! record number in direct access file corresonding to specific line ! LINBEG 
    
    ! two extended neighbour subintervals: 
    ! [extStartDeltaWV1; extEndDeltaWV1] and [extStartDeltaWV2, extEndDeltaWV2] are highly overlap
    ! so when going through the first interval, it is needed to set at which lineIdx to start next subinterval.
    ! capWV is used to determine this, see the implementation of it in the line-by-line scheme
    real(kind=DP) :: capWV
    
    logical :: isFirstSubinterval

    !------------------------------------------------------------------------------------------------------------------!

    ! Start a timer to track the machine time
    call cpu_time(startTime)

    ! Get a number of command-line arguments
    argc = command_argument_count()

    ! Check if the number of arguments is not sufficient
    if (argc < 8) then
        print *, 'Insufficient number of arguments.'
        print *, 'Expected at least 8 arguments, but received ', argc
        print *, 'Usage: marfa Molecule StartWV EndWV DatabaseSlug & 
                    CutOff ChiFactorFuncName TargetValue AtmProfileFile'
        print *, 'Molecule: CO2, H2O'
        print *, 'StartWV: 10 - 20000'
        print *, 'EndWV: 10 - 20000, and greater than StartWV'
        print *, 'DatabaseSlug: basefilenames from the data/databases folder'
        print *, 'CutOff: integer value of a line cutoff condition in cm-1'
        print *, 'ChiFactorFuncName: name of the wing correction function from &
                     the ChiFactors.f90 module'
        print *, 'TargetValue: ACS (for cross-section in cm^2/molecule) or &
                     VAC (for volume absorption coefficient in km-1)'
        stop 1
    end if

    !------------------------------------------------------------------------------------------------------------------!

    ! Establishing input parameters based on command line arguments
    do l = 1, command_argument_count()
        select case (l)
        case (1)
            call get_command_argument(l, inputMolecule)
        case (2)
            call get_command_argument(l, startWVcla)
            startWVclaTrimmed = trim(startWVcla)
            read(startWVclaTrimmed, *) startWV
        case (3)
            call get_command_argument(l, endWVcla)
            endWVclaTrimmed = trim(endWVcla)
            read(endWVclaTrimmed, *) endWV
        case (4)
            call get_command_argument(l, databaseSlug)
        case (5)
            call get_command_argument(l, cutOffcla)
            cutOffclaTrimmed = trim(cutOffcla)
            read(cutOffclaTrimmed, *) cutOff
        case (6)
            call get_command_argument(l, chiFactorFuncName)
            call fetchChiFactorFunction
        case (7)
            call get_command_argument(l, targetValuecla)
            ! Set ACS or VAC with validation
            targetValue = trim(targetValuecla)
            select case (targetValue)
            case ('ACS', 'VAC')
                ! Valid targetValue; proceed as normal
            case default
                print *, 'Error: Invalid targetValue "', targetValue, '". & 
                            Must be either "ACS" (Absorption cross-section) or &
                             "VAC" (Volume absorption coefficient).'
                stop 2
            end select
        case (8)
            call get_command_argument(l, atmProfileFile)
            call readAtmosphericParameters
        end select
    end do

    ! reads and initializes the TIPS array
    call readTIPS()

    ! setting a line shape function (choose from Shapes.f90 module)
    shapeFuncPtr => voigt
    !------------------------------------------------------------------------------------------------------------------!

    ! Directory where the PT-tables are stored is generated automatcally and contains a timestamp in its name
    call date_and_time(values=dateTimeValues)
    year = dateTimeValues(1)
    month = dateTimeValues(2)
    day = dateTimeValues(3)
    hour = dateTimeValues(5)
    minute = dateTimeValues(6)
    second = dateTimeValues(7)
    write(timestamp, '(I4, I2.2, I2.2, "_", I2.2, I2.2, I2.2, I2.2)') & 
                        year, month, day, hour, minute, second
    
    ! Format numeric parameters to remove decimal points and convert to integers
    write(formattedStartWV, '(I0)') int(startWV)
    write(formattedEndWV, '(I0)') int(endWV)

    ! Assemble a subdirectory name: Molecule_StartWV-EndWV_Timestamp, absolute and relative paths
    subDirName = trim(inputMolecule) // "_" // trim(formattedStartWV) // "-" //&
                     trim(formattedEndWV) // "_" // trim(timestamp)
    parentDir = trim(adjustl(rootDirName)) // "/ptTables/"
    fullSubDirPath = trim(adjustl(parentDir)) // trim(adjustl(subDirName))
    
    ! Running a makedir command to create directories
    mkdirCommand = 'mkdir "' // trim(fullSubDirPath) // '"'
    call execute_command_line(mkdirCommand, wait=.true., exitstat=status)
    if (status /= 0) then
        print *, "Error: Failed to create directory ", trim(fullSubDirPath)
        stop 3
    else
        print *, "Directory for storing PT-tables is created: ", trim(fullSubDirPath)
    end if

    ! Establishing technical infoFile for storin information about the parameters of the run
    infoFilePath = trim(fullSubDirPath) // '/info.txt'
    open(infoUnit, file=infoFilePath, status='replace', action='write', iostat=status)

    if (status /= 0) then
        print *, "Error: Unable to create info file at ", trim(infoFilePath)
        stop 4
    end if

    ! Write command-line arguments to the info file
    write(infoUnit, '(A)') 'Command-Line Arguments:'
    write(infoUnit, '(A,A)') 'Input Molecule: ', trim(inputMolecule)
    write(infoUnit, '(A,A)') 'Start Wavenumber: ', trim(startWVclaTrimmed)
    write(infoUnit, '(A,A)') 'End Wavenumber: ', trim(endWVclaTrimmed)
    write(infoUnit, '(A,A)') 'Cut Off: ', trim(cutOffclaTrimmed)
    write(infoUnit, '(A,A)') 'Chi Factor Function Name: ', trim(chiFactorFuncName)
    write(infoUnit, '(A,A)') 'Target Value: ', trim(targetValue)
    write(infoUnit, '(A,A)') 'Atmospheric Profile File: ', trim(atmProfileFile)
    close(infoUnit)
    print *, "Info file created at: ", trim(infoFilePath)

    ! Storing information about the time of a last calculation (needed for python scripts to save the processed data and plots)
    latestRunFilePath = trim(adjustl(parentDir)) // 'latest_run.txt'
    open(unit=latestRunUnit, file=trim(latestRunFilePath), &
                status='replace', action='write', iostat=status)
    if (status /= 0) then
        print *, "Error: Unable to create latest run file at ", trim(latestRunFilePath)
        stop 5
    end if
    write(latestRunUnit, '(A)') subDirName
    close(latestRunUnit)
    
    ! DEBUG SECTION !
    ! print *, inputMolecule
    ! print *, startWV
    ! print *, endWV
    ! print *, cutOff
    ! print *, chiFactorFuncName
    ! print *, targetValue
    ! print *, atmProfileFile

    call getSpeciesCode(inputMolecule, moleculeIntCode, DBfileExtension)
    if (moleculeIntCode == 0) then
        print *, 'ERROR: unsupported molecule: ', trim(adjustl(inputMolecule))
        stop 8
    end if

    call determineStartingSpectralLine(databaseSlug, startingLineIdx, startingLineWV)
    
    ATMOSPHERIC_LEVELS_LOOP: do levelsIdx = 1, levels
        ! OUTER LOOP: over the atmospheric levels. After each iteration, PT-table file is generated for
        ! the levelsIdx level
        
        ! Initializing the pressure (total), temperature and density (of the species) at current atm level
        pressure = pressureArray(levelsIdx)
        temperature = temperatureArray(levelsIdx)
        density = densityArray(levelsIdx)
        
        ! Calculation of the self and foreign pressures based on the Loschmidt formula
        ! Alternative -- using Dalton's law and Mendeleev equation
        unitlessT = refTemperature / temperature
        pSelf = density * 10. / LOSCHMIDT * temperature/stTemperature
        pForeign = pressure - pSelf

        ! Construction of the extenstion for the PT-table output file, reflecting an atmospheric level
        levelLabel = '___'
        if ( levelsIdx < 10 ) then
            write(levelLabel(1:1), '(I1)') levelsIdx
        else
            if ( levelsIdx < 100 ) then
                write(levelLabel(1:2), '(I2)') levelsIdx
            else
                write(levelLabel(1:3), '(I3)') levelsIdx
            end if
        end if
        
        ! Direct access files section: be aware about the OS compatability
        open(outputUnit, access='DIRECT', form='UNFORMATTED', recl=NT*4, &
            file=trim(fullSubDirPath)//'/'//levelLabel//'.ptbin')
        ! RECL = NT for old Windows Fortrans !
        
        ! the whole interval [startWV, endWV] is separated into subintervals of the same length of deltaWV
        ! each record in the output file corresponds to the one subinterval
        startDeltaWV = startWV ! left boundary of the first subinterval
        endDeltaWV = startDeltaWV + deltaWV ! right boundary of the first subinterval

        isFirstSubinterval = .true.

        SUBINTERVALS_LOOP: do while (startDeltaWV < endWV)

            ! Relation between record number and left boundary of a subinterval
            ! TODO: rewrite when working on introducing the dynamic resolution and fixing files sizes issue
            outputRecNum = (startDeltaWV + 1.0) / 10.0 ! *** (0.0 -> 0 , 10.0 -> 1,..., 560.0 -> 56, etc.)

            ! Defining the boundaries of extended subinterval (add cutOff) 
            extStartDeltaWV = startDeltaWV - cutOff
            extEndDeltaWV = startDeltaWV + deltaWV + cutOff
            
            ! Proceed to calculation inside subinterval !
            if (isFirstSubinterval) then
                lineIdx = startingLineIdx
                isFirstSubinterval = .false.
            end if

            ! Defining capWV: wavenumber to determine from which spectral line to start calculation
            ! of the next subinterval
            ! TODO: introduce pointers logic for that -- do when fixing the issue with subintervals overlap
            capWV = extStartDeltaWV + deltaWV
            ! The same is:
            ! capWV = endDeltaWV - cutOff
            if (capWV <= startingLineWV) capWV = startingLineWV ! for consistency, because cutOff might be large
            
            ! Arrays for storing absorption values on all grids
            ! must be initialized to zero at the beginning of each calculation inside subinterval
            call resetAbsorptionGridValues()
    
            ! DEBUG SECTION !
            ! print *, lineIdx
            ! print *, capWV
            
            ! Proceed to this subroutine for reading spectral features line-by-line inside a subinterval
            call lblCalculation(lineIdx, capWV)
            
            ! DEBUG SECTION !
            ! print *, 'lineIdx after LBL: ', lineIdx
    
            call cascadeInterpolation()

            ! WRITE OPERATION OF THE ABSORPTION SIGNATURES TO THE OUTPUT FILE !
            if (targetValue == 'VAC') then
                ! Multiply on density if want to calculate a volume absorption coefficient
                write(outputUnit, rec=outputRecNum) density * RK
            elseif (targetValue == 'ACS') then
                write(outputUnit, rec=outputRecNum) RK
            end if

            ! switching to the next subinterval
            startDeltaWV = startDeltaWV + deltaWV
            endDeltaWV = startDeltaWV + deltaWV
            
        end do SUBINTERVALS_LOOP
        
        close(outputUnit)
        
        ! for real time tracking how many levels has been processed:
        print *, levelsIdx, ' of ', levels, ' atmospheric levels is being processed ...' 
    end do ATMOSPHERIC_LEVELS_LOOP
    
    print *, ' *** Congratulations! PT-tables have been calculated! ***'
    deallocate(heightArray)
    deallocate(pressureArray)
    deallocate(densityArray)
    deallocate(temperatureArray)
    deallocate(TIPS)
    
    call cpu_time(endTime)
    print *, "Took: ", endTime - startTime, " seconds"


contains


    subroutine determineStartingSpectralLine(DBName, lineBegIdx, lineWVBeg)
        implicit none
        ! database file basename (e.g. "HITRAN2020"), see `data/databases/` directory
        character(len=20), intent(in) :: DBName
        integer, intent(out) :: lineBegIdx
        real(kind=DP), intent(out) :: lineWVBeg
        
        real(kind=DP) :: iterLineWV 
        integer :: iterRecNum

        lineBegIdx = 1

        databaseFile = 'data'//'/'//'databases'//'/'//trim(adjustl(DBName))//'.'//DBfileExtension

        ! Reading the direct access files (Linux, MacOS, modern Windows(?)). 
        ! Comment out this line for old Windows
        open(databaseFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=trim(databaseFile))
        
        ! Uncomment if the direct access file was created under old Windows(?)
        ! open(7777, access='DIRECT', form='UNFORMATTED', recl=9, file=databaseFile)

        ! This section identifies a starting spectral line from which to proceed with line-by-line scheme
        ! based on the left boundary of the initial spectral interval [startWV, endWV]
        read(databaseFileUnit, rec=lineBegIdx) lineWVBeg

        ! This `if` block is needed to determine from which spectral line in database to start calculation !
        ! See the `lineBegIdx` and `lineWVBeg` variables
        if (startWV > cutOff) then
            ! if the startWV is e.g. 300 cm-1 and the cutOff is e.g. 25 cm-1, then
            ! the first spectral line to be counted must be the first line with transition 
            ! wavenumber bigger than extStartWV = 300 - 25 = 275 cm-1
            extStartWV = startWV - cutOff
            iterRecNum = lineBegIdx
            iterLineWV = lineWVBeg
            do while(iterLineWV <= extStartWV)
                iterRecNum = iterRecNum + 1
                read(databaseFileUnit, rec=iterRecNum) iterLineWV
            end do
            lineBegIdx = iterRecNum
            lineWVBeg = iterLineWV
        else 
            ! if the startWV is e.g. 20 cm-1, but cutOff is 125 cm-1, then
            ! proceed calculation from the first spectral line presented in the database file
            ! lineWVBeg and lineBegIdx are set to initial values
        end if

        ! simple check for consistency of the database file
        if ((abs(startWV - cutOff - lineWVBeg) > 25.) .and. (startWV > cutOff)) then
            print *, "ATTENTION: database file might be insufficient for input spectral interval:"
            print *, "Your input left boundary - cutOff condition: ", startWV - cutOff, " cm-1"
            print *, "Line-by-line scheme starts from: ", lineWVBeg, " cm-1"
        end if
    end subroutine determineStartingSpectralLine

    
    subroutine readTIPS()
        implicit none
        character(len=300), parameter :: TIPSFile = 'data/QofT_formatted.dat' ! path to the file with TIPS data
        integer, parameter :: TIPSUnit = 5467 ! unit for file with TIPS data
        integer :: nIsotopes, nTemperatures ! number of different isotopes and temperatures in the TIPS file
        integer :: stSumTIdx, stSumIsoIdx  ! loop indices for partition sums: temperatures and isotopes
        ! This subroutine opens a file containing TIPS data, allocates an array for storage, 
        ! and populates the array with data in a loop
        open(unit=TIPSUnit, file='data/TIPS/TIPS.dat', status='old', action='read')
        read(TIPSUnit, *) nIsotopes, nTemperatures
        allocate(TIPS(nIsotopes, nTemperatures))
        do stSumIsoIdx = 1, nIsotopes
            read(TIPSUnit, *) (TIPS(stSumIsoIdx, stSumTIdx), stSumTIdx=1, nTemperatures)
        end do
        close(TIPSUnit)
    end subroutine readTIPS
    
    
    subroutine getSpeciesCode(species, codeInt, codeStr)
        ! Subroutine to map species to both integer and string codes accordnig to the HITRAN coding system
        implicit none
        character(len=*), intent(in) :: species  ! molecule title as string
        integer, intent(out) :: codeInt ! output code as an integer
        character(len=2), intent(out) :: codeStr ! output code as a string

        select case (trim(adjustl(species)))
            case ("H2O")
                codeInt = 1
                codeStr = "01"
            case ("CO2")
                codeInt = 2
                codeStr = "02"
            case ("O3")
                codeInt = 3
                codeStr = "03"
            case ("N2O")
                codeInt = 4
                codeStr = "04"
            case ("CO")
                codeInt = 5
                codeStr = "05"
            case ("CH4")
                codeInt = 6
                codeStr = "06"
            case ("O2")
                codeInt = 7
                codeStr = "07"
            case ("NO")
                codeInt = 8
                codeStr = "08"
            case ("SO2")
                codeInt = 9
                codeStr = "09"
            case ("NO2")
                codeInt = 10
                codeStr = "10"
            case ("NH3")
                codeInt = 11
                codeStr = "11"
            case ("HNO3")
                codeInt = 12
                codeStr = "12"
            case default
                codeInt = 0
                codeStr = "00"  ! Species not found
        end select
    end subroutine getSpeciesCode


end program main
