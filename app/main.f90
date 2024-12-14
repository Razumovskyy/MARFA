program main
    use Constants
    use IO
    use Atmosphere
    use Spectroscopy
    use LBL
    use Mesh
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
    real(kind=DP) startTime, endTime ! times for measuring a machine time of the run
    character(len=20) :: timestamp ! 
    integer :: status
    character(len=3) :: levelLabel ! unique identifier for each atmospheric level: ('1__','2__',...,'100',...)
    integer :: l ! loop variable

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
    call readTIPS

    ! setting a line shape function
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

    !------------------------------------------------------------------------------------------------------------------!

    ! TODO:(!!) optimization proposition. Reduce reading operations because of highly overlaping intervals: 
    !       [extStartWV1, extEndWV1] and [extStartWV2, extEndWV2]
    
    ! DEBUG SECTION !
    ! print *, inputMolecule
    ! print *, startWV
    ! print *, endWV
    ! print *, cutOff
    ! print *, chiFactorFuncName
    ! print *, targetValue
    ! print *, atmProfileFile

    ! Defining the grid constants (see the Mesh.f90 module)
    ! TODO:(!) move out from the main.f90 file, it initialized only once
    EPS = 0.0
    H0 = STEP
    H1 = H0/2.
    H2 = H1/2.
    H3 = H2/2.
    H4 = H3/2.
    H5 = H4/2.
    H6 = H5/2.
    H7 = H6/2. 
    H8 = H7/2. 
    H9 = H8/2. 
    H = H9/4.

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
        SUBINTERVALS_LOOP: do while (startDeltaWV < endWV)
            ! loop over subintervals !

            ! Relation between record number and left boundary of a subinterval
            ! TODO:(!!) rewrite when working on introducing the dynamic resolution and fixing file sizes issue
            outputRecNum = (startDeltaWV + 1.0) / 10.0 ! *** (0.0 -> 0 , 10.0 -> 1,..., 560.0 -> 56, etc.)

            ! TODO:(!) that couldn't happen -- remove
            ! if ( extStartDeltaWV < startWV ) then
            !     extStartDeltaWV = startDeltaWV
            ! end if

            ! Proceed to calculation inside subinterval !
            call processSpectra(inputMolecule, databaseSlug, levelsIdx)

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

    subroutine processSpectra(molecule, slug, atmosphericLoopLevel)
        ! TODO: (!) remove `molecule` from input arguments (used only in first call) 
        implicit none
        character(len=5), intent(in) :: molecule ! inputMolecule as a string
        integer, intent(in) :: atmosphericLoopLevel ! atmospheric level number as passed from the outer loop
        
        ! database file basename (e.g. "HITRAN2020"), see `data/databases/` directory
        character(len=20), intent(in) :: slug
        
        ! ---------------------------------------------------------------------------- !
        ! database file extension based on the species title (input molecule), see `data/databases/` directory
        character(len=2) :: DBfileExtension

        ! this save is to decouple the logic run only for the first call
        logical, save :: isFirstCall = .true.
        
        ! this save is for the keeping the current level. 
        ! It is calculated on the first call and then updates if the level changes in the calling procedure
        ! current atmospheric level at previous subroutine call, used to fix when the atmospheric level changes
        integer, save :: currentLevel ! JOLD
        
        ! this save is for the loop over height levels
        ! after the step into the new level, LBL starts from startingLineIdx
        ! Index (record number in database file) from which the whole calculation goes
        ! Determined by startWV and cutOff
        integer, save :: startingLineIdx ! LINBEG0
        
        ! Transition wavenumber of a spectral line at index `startingLineIdx ! VAA0
        ! Mostly it is used to determine capWV
        real(kind=DP), save :: startingLineWV 

        ! this save is for the loop over subintervals
        ! after the step into the next delta interval, we need to keep the lineIdx from the LBL calc
        ! `lineIdx` is used as an iterator over spectral lines in SUBINTERVALS_LOOP
        integer, save :: lineIdx ! record number in direct access file corresonding to specific line ! LINBEG 

        ! two extended subsequent subintervals: 
        ! [extStartDeltaWV1; extEndDeltaWV1] and [extStartDeltaWV2, extEndDeltaWV2] are highly overlap
        ! so when going through the first interval, it is needed to set at which lineIdx the next interval starts.
        ! capWV is used to determine this, see the implementation of it in the modernLBL subroutine
        real(kind=DP), save :: capWV

        ! auxillary variables
        integer :: iterRecNum ! loop variable
        real(kind=DP) :: iterLineWV ! loop dummy variable
        integer :: I, J, M ! loop variables mainly for grid calculations
        
        ! TODO:(!) moving `ifFirstCall` part to the separate subroutine, because
        ! it is only for reading from the spectral file and must be done only once. Do it during
        ! implementing parallelization with OpenMP
         
        ! DEBUG SECTION !
        ! real(kind=DP) :: VFISH ! *** ! extEndDeltaWV
        ! real(kind=DP) :: VA ! *** ! extStartDeltaWV
        ! real(kind=DP) :: VS ! the same as startDeltaWV
        ! real(kind=DP) :: VR4 ! the same as startDeltaWV
        
        if (isFirstCall) then
            
            isFirstCall = .false.
            startingLineIdx = 1
            currentLevel = 0
            lineIdx = 0
            
            ! TODO: (!) move out from this subroutine (do it when moving the first call) 
            call get_species_code(molecule, moleculeIntCode, DBfileExtension)

            if (moleculeIntCode == 0) then
                print *, 'ERROR: unsupported molecule: ', trim(adjustl(molecule))
                stop 8
            end if

            databaseFile = 'data'//'/'//'databases'//'/'//trim(adjustl(slug))//'.'//DBfileExtension

            ! Reading the direct access files (Linux, MacOS, modern Windows(?)). 
            ! Comment out this line for old Windows
            open(databaseFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=trim(databaseFile))
            
            ! Uncomment if the direct access file was created under old Windows(?)
            ! open(7777, access='DIRECT', form='UNFORMATTED', recl=9, file=databaseFile)

            ! This section identifies a starting spectral line from which to proceed with line-by-line scheme
            ! based on the left boundary of the initial spectral interval [startWV, endWV]
            ! TODO:(!) optimize this block and startingLineIdx and startingLineWV variables
            read(databaseFileUnit, rec=startingLineIdx) startingLineWV

            
            ! This `if` block is needed to determine from which spectral line in database to start calculation !
            ! See the `startingLineIdx` and `startingLineWv` variables
            if (startWV > cutOff) then
                ! if the startWV is e.g. 300 cm-1 and the cutOff is e.g. 25 cm-1, then
                ! the first spectral line to be counted must be the first line with transition 
                ! wavenumber bigger than extStartWV = 300 - 25 = 275 cm-1
                extStartWV = startWV - cutOff
                iterRecNum = startingLineIdx
                iterLineWV = startingLineWV
                do while(iterLineWV <= extStartWV)
                    iterRecNum = iterRecNum + 1
                    read(databaseFileUnit, rec=iterRecNum) iterLineWV
                end do
                startingLineIdx = iterRecNum
                startingLineWV = iterLineWV
            else 
                ! if the startWV is e.g. 20 cm-1, but cutOff is 125 cm-1, then
                ! proceed calculation from the first spectral line presented in the database file
                ! startingLineWV and startingLineIdx are set to initial values
            end if

            ! simple check for consistency of the database file
            if ((abs(startWV - cutOff - startingLineWV) > 25.) .and. (startWV > cutOff)) then
                print *, "ATTENTION: database file might be insufficient for input spectral interval:"
                print *, "Your input left boundary - cutOff condition: ", startWV - cutOff, " cm-1"
                print *, "Line-by-line scheme starts from: ", startingLineWV, " cm-1"
            end if
        end if

        ! When the ATMOSPHERIC_LOOP steps onto the new atmospheric level, then the lineIdx must be reset 
        ! to the the startingLineIdx value calculated in the firstCall section, because when it is new
        ! atmospheric level, calculation goes again from startWV.
        if ( currentLevel /= atmosphericLoopLevel ) then
            ! DEBUG SECTION !
            ! write(*,*) 'Level change !'
            ! write(*,*) 'current Level: ', currentLevel
            ! write(*,*) 'loop Level: ', atmosphericLoopLevel
            currentLevel = atmosphericLoopLevel
            lineIdx = startingLineIdx
        end if

        ! TODO:(?) move to the main.f90 file in the inner loop
        
        ! Defining the boundaries of extended subinterval (add cutOff) 
        extStartDeltaWV = startDeltaWV - cutOff
        extEndDeltaWV = startDeltaWV + deltaWV + cutOff
        
        ! Defining capWV: wavenumber to determine from which spectral line to start calculation
        ! of the next subinterval
        ! TODO:(!!) introduce pointers logic for that -- do when fixing the issue with subintervals overlap
        capWV = extStartDeltaWV + deltaWV
        ! The same is:
        ! capWV = endDeltaWV - cutOff
        if (capWV <= startingLineWV) capWV = startingLineWV ! for consistency, because cutOff might be large

        !------------------------------------------------------------------------------------------------------------------!
        
        ! Note: Grid arrays for storing absorption values on all grids
        ! are initialized to zero at the beginning of each calculation inside subinterval.
        
        ! TODO:(!) consider moving this section to inner subroutines
        RK = 0.0
        RK0 = 0.0; RK0L = 0.0; RK0P = 0.0
        RK1 = 0.0; RK1L = 0.0; RK1P = 0.0
        RK2 = 0.0; RK2L = 0.0; RK2P = 0.0
        RK3 = 0.0; RK3L = 0.0; RK3P = 0.0
        RK4 = 0.0; RK4L = 0.0; RK4P = 0.0
        RK5 = 0.0; RK5L = 0.0; RK5P = 0.0
        RK6 = 0.0; RK6L = 0.0; RK6P = 0.0
        RK7 = 0.0; RK7L = 0.0; RK7P = 0.0
        RK8 = 0.0; RK8L = 0.0; RK8P = 0.0
        RK9 = 0.0; RK9L = 0.0; RK9P = 0.0

        !------------------------------------------------------------------------------------------------------------------!

        ! DEBUG SECTION !
        ! print *, lineIdx
        ! print *, capWV
        
        ! Proceed to this subroutine for reading of spctral features and line summation in the subinterval
        call modernLBL(lineIdx, capWV)
        
        ! DEBUG SECTION !
        ! write(*,*) 'lineIdx after moderLBL: ', lineIdx
        
        !------------------------------------------------------------------------------------------------------------------!
        
        ! Final interpolation scheme (?)

        ! TODO:(!) Create subroutine to remove repetitions
        ! (be careful that last loop differs from NT0-NT8) 
        
        ! TODO:(!) consider moving this section to inner subroutines
        ! TODO:(!) establish block and define loop variables there
        do J = 1, NT0
            I = J*2 - 1
            RK1P(I) = RK1P(I) + RK0P(J)
            RK1(I) = RK1(I) + RK0P(J)*0.375 + RK0(J)*0.75 - RK0L(J)*0.125
            RK1L(I) = RK1L(I) + RK0(J)
            M = I + 1
            RK1P(M) = RK1P(M) + RK0(J)
            RK1(M) = RK1(M) + RK0L(J)*0.375 + RK0(J)*0.75 - RK0P(J)*0.125
            RK1L(M) = RK1L(M) + RK0L(J)
        end do

        do J = 1, NT1
            I = J*2 - 1
            RK2P(I) = RK2P(I) + RK1P(J)
            RK2(I) = RK2(I) + RK1P(J)*0.375 + RK1(J)*0.75 - RK1L(J)*0.125
            RK2L(I) = RK2L(I) + RK1(J)
            M = I + 1
            RK2P(M) = RK2P(M) + RK1(J)
            RK2(M) = RK2(M) + RK1L(J)*0.375 + RK1(J)*0.75 - RK1P(J)*0.125
            RK2L(M) = RK2L(M) + RK1L(J)
        end do

        do J = 1, NT2
            I = J*2 - 1
            RK3P(I) = RK3P(I) + RK2P(J)
            RK3(I) = RK3(I) + RK2P(J)*0.375 + RK2(J)*0.75 - RK2L(J)*0.125
            RK3L(I) = RK3L(I) + RK2(J)
            M = I + 1
            RK3P(M) = RK3P(M) + RK2(J)
            RK3(M) = RK3(M) + RK2L(J)*0.375 + RK2(J)*0.75 - RK2P(J)*0.125
            RK3L(M) = RK3L(M) + RK2L(J)
        end do

        do J = 1, NT3
            I = J*2 - 1
            RK4P(I) = RK4P(I) + RK3P(J)
            RK4(I) = RK4(I) + RK3P(J)*0.375 + RK3(J)*0.75 - RK3L(J)*0.125
            RK4L(I) = RK4L(I) + RK3(J)
            M = I + 1
            RK4P(M) = RK4P(M) + RK3(J)
            RK4(M) = RK4(M) + RK3L(J)*0.375 + RK3(J)*0.75 - RK3P(J)*0.125
            RK4L(M) = RK4L(M) + RK3L(J)
        end do

        do J = 1, NT4
            I = J*2 - 1
            RK5P(I) = RK5P(I) + RK4P(J)
            RK5(I) = RK5(I) + RK4P(J)*0.375 + RK4(J)*0.75 - RK4L(J)*0.125
            RK5L(I) = RK5L(I) + RK4(J)
            M = I + 1
            RK5P(M) = RK5P(M) + RK4(J)
            RK5(M) = RK5(M) + RK4L(J)*0.375 + RK4(J)*0.75 - RK4P(J)*0.125
            RK5L(M) = RK5L(M) + RK4L(J)
        end do

        do J = 1, NT5
            I = J*2 - 1
            RK6P(I) = RK6P(I) + RK5P(J)
            RK6(I) = RK6(I) + RK5P(J)*0.375 + RK5(J)*0.75 - RK5L(J)*0.125
            RK6L(I) = RK6L(I) + RK5(J)
            M = I + 1
            RK6P(M) = RK6P(M) + RK5(J)
            RK6(M) = RK6(M) + RK5L(J)*0.375 + RK5(J)*0.75 - RK5P(J)*0.125
            RK6L(M) = RK6L(M) + RK5L(J)
        end do
        
        do J = 1, NT6
            I = J*2 - 1
            RK7P(I) = RK7P(I) + RK6P(J)
            RK7(I) = RK7(I) + RK6P(J)*0.375 + RK6(J)*0.75 - RK6L(J)*0.125
            RK7L(I) = RK7L(I) + RK6(J)
            M = I + 1
            RK7P(M) = RK7P(M) + RK6(J)
            RK7(M) = RK7(M) + RK6L(J)*0.375 + RK6(J)*0.75 - RK6P(J)*0.125
            RK7L(M) = RK7L(M) + RK6L(J)
        end do

        do J = 1, NT7
            I = J*2 - 1
            RK8P(I) = RK8P(I) + RK7P(J)
            RK8(I) = RK8(I) + RK7P(J)*0.375 + RK7(J)*0.75 - RK7L(J)*0.125
            RK8L(I) = RK8L(I) + RK7(J)
            M = I + 1
            RK8P(M) = RK8P(M) + RK7(J)
            RK8(M) = RK8(M) + RK7L(J)*0.375 + RK7(J)*0.75 - RK7P(J)*0.125
            RK8L(M) = RK8L(M) + RK7L(J)
        end do

        do J = 1, NT8
            I = J*2 - 1
            RK9P(I) = RK9P(I) + RK8P(J)
            RK9(I) = RK9(I) + RK8P(J)*0.375 + RK8(J)*0.75 - RK8L(J)*0.125
            RK9L(I)= RK9L(I) + RK8(J)
            M = I + 1
            RK9P(M) = RK9P(M) + RK8(J)
            RK9(M) = RK9(M) + RK8L(J)*0.375 + RK8(J)*0.75 - RK8P(J)*0.125
            RK9L(M)= RK9L(M) + RK8L(J)
        end do

        I=1
        do J = 1, NT9
            I = I + 1
            RK(I) = RK(I) + (RK9P(J)*0.375 + RK9(J)*0.75 - RK9L(J)*0.125)
            I = I + 1
            RK(I) = RK(I) + RK9(J)
            I = I + 1
            RK(I) = RK(I) + (RK9L(J)*0.375 + RK9(J)*0.75 - RK9P(J)*0.125)
            I = I + 1
            RK(I) = RK(I) + RK9L(J)
        end do

        do J = 1, NT
            if (RK(J) < 0.) RK(J)=0.
        end do
        !------------------------------------------------------------------------------------------------------------------!
    end subroutine processSpectra

    subroutine get_species_code(species, code_int, code_str)
        ! Subroutine to map species to both integer and string codes accordnig to the HITRAN coding system
        implicit none
        character(len=*), intent(in) :: species  ! molecule title as string
        integer, intent(out) :: code_int ! output code as an integer
        character(len=2), intent(out) :: code_str ! output code as a string

        select case (trim(adjustl(species)))
            case ("H2O")
                code_int = 1
                code_str = "01"
            case ("CO2")
                code_int = 2
                code_str = "02"
            case ("O3")
                code_int = 3
                code_str = "03"
            case ("N2O")
                code_int = 4
                code_str = "04"
            case ("CO")
                code_int = 5
                code_str = "05"
            case ("CH4")
                code_int = 6
                code_str = "06"
            case ("O2")
                code_int = 7
                code_str = "07"
            case ("NO")
                code_int = 8
                code_str = "08"
            case ("SO2")
                code_int = 9
                code_str = "09"
            case ("NO2")
                code_int = 10
                code_str = "10"
            case ("NH3")
                code_int = 11
                code_str = "11"
            case ("HNO3")
                code_int = 12
                code_str = "12"
            case default
                code_int = 0
                code_str = "00"  ! Species not found
        end select
    end subroutine get_species_code
end program main
