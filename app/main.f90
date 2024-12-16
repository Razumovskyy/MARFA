program main
    use Constants
    use Atmosphere
    use Spectroscopy
    use LBL
    use Grids
    use EnvironmentSetup
    implicit none

    ! Output files units !
    integer, parameter :: outputUnit = 47 ! output file unit where PT-tables are stored
    
    ! Record number in the output file
    integer :: outputRecNum 
    
    ! Technical variables !
    real(kind=DP) startTime, endTime ! times for measuring a machine time
    character(len=3) :: levelLabel ! unique identifier for each atmospheric level: ('1__','2__',...,'100',...)

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

    
    ! Start a timer to track the machine time
    call cpu_time(startTime)

    call readCommandLineArguments()
    call readAtmosphericParameters()
    
    call setShapeFunction()
    call setChiFactorFunction()

    ! reads and initializes the TIPS array
    call readTIPS()

    call initialiseDirectories()
    call initialiseLogFiles()
    
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
            else if (targetValue == 'ACS') then
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


end program main
