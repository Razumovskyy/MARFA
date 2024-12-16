module EnvironmentSetup
    use Constants
    use LBL
    use Atmosphere
    implicit none
    
    ! Input command line arguments and their trimmed values as strings !
    integer :: argc ! number of the passed command line arguments
    character(len=20) :: databaseSlug
    character(len=10) :: inputMolecule   
    character(len=10) :: startWVcla, endWVcla, startWVclaTrimmed, endWVclaTrimmed
    character(len=10) :: cutOffcla, cutOffclaTrimmed
    character(len=30) :: targetValue, targetValuecla

    ! Parameters and variables for constructing directories names !
    character(len=200), parameter :: rootDirName = 'output'
    character(len=300) :: subDirName
    character(len=100) :: formattedStartWV, formattedEndWV
    character(len=200) :: parentDir
    character(len=300) :: fullSubDirPath
    character(len=500) :: mkdirCommand
    character(len=300) :: infoFilePath
    character(len=300) :: latestRunFilePath
    integer, parameter :: infoUnit = 67
    integer, parameter :: latestRunUnit = 87
    character(len=20) :: timestamp ! 
    integer :: dateTimeValues(8)
    integer :: year, month, day, hour, minute, second
    integer :: status

contains


    subroutine readCommandLineArguments()
        implicit none
        integer :: l

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
            end select
        end do
    end subroutine readCommandLineArguments

    
    subroutine initialiseDirectories()
        implicit none
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
    end subroutine initialiseDirectories
        

    subroutine initialiseLogFiles
        implicit none
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
    end subroutine initialiseLogFiles


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
    
    
    pure subroutine getSpeciesCode(species, codeInt, codeStr)
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


end module EnvironmentSetup