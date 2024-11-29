module Atmosphere
    use Constants
    implicit none
    
    integer, parameter :: atmProfileUnit = 777 ! unit for the atmospheric file
    character(len=300) :: atmProfileFile ! title of the file with atmospheric data located in `data/Atmospheres`
    character(len=300) :: fullNameAtmProfile ! full path to the atmospheric file
    character(len=100) :: atmTitle ! specific atmosphere title (in the header of the atmospheric file)
    integer :: levels ! number of levels of the considered atmosphere ( located in a file header )
    real, allocatable :: heightArray(:), pressureArray(:), temperatureArray(:), densityArray(:) ! PPP, TTT, RORO (legacy)
    integer, parameter :: levelsThreshold = 200 

    ! macroscopic parameters on the current atmospheric level !
    ! temperature: temperauture on the current atmospheric level [K]
    ! pressure: total pressure on the current atmospheric level [atm]
    ! density: number density of the considered species: [(number of molecules)/(cm^2 * km)]
    !           such number density units are needed for having volume absorption coefficient in km-1
    real :: temperature, pressure, density

    real :: pSelf ! [atm] -- partial pressure of considered gaseous species
    real :: pForeign ! [atm] -- foreign pressure (p(total) - pSelf)
    real :: unitlessT ! see the main.f90 implementation
    
    ! Technical variables
    integer :: ios
    integer :: levelsIdx ! loop index

contains
    ! TODO:(?) consider moving loop over atmospheric levels here
    subroutine readAtmosphericParameters()
        
        fullNameAtmProfile = 'data/Atmospheres/'//atmProfileFile

        ! Reading the header of the atmospheric file: atmTitle, levels !
        open(atmProfileUnit, file=fullNameAtmProfile, status='old')
        read(atmProfileUnit, '(A20)') atmTitle
        read(atmProfileUnit, *) levels
        if (levels > 200) then
            write(*, '(A, I3, A)') 'WARNING: input number of atmospheric levels is &
                                    bigger than ', levelsThreshold
            stop 6
        end if

        ! Allocation and population of arrays for macroscopic parameters !
        allocate(heightArray(levels))
        allocate(pressureArray(levels))
        allocate(temperatureArray(levels))
        allocate(densityArray(levels))
        do levelsIdx = 1, levels
            read(atmProfileUnit, *, iostat=ios) heightArray(levelsIdx), pressureArray(levelsIdx), &
            temperatureArray(levelsIdx), densityArray(levelsIdx)
            if (ios /= 0) then
                print *, 'Error: Unable to read data from file "', trim(fullNameAtmProfile), '".'
                print *, 'Check the data format in the atmospheric file'
                close(atmProfileUnit)
                stop 7
            end if
        end do
        close(atmProfileUnit)

    end subroutine readAtmosphericParameters
end module Atmosphere
