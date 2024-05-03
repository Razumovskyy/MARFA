module Atmosphere
    use Constants
    implicit none
    character(len=20) :: atmProfileFile ! NEW ! ** ! ATM
    character(len=20) :: atmTitle
    integer, parameter :: atmProfileUnit = 777
    integer :: numGasSpecies ! NGS ! ** ! number of species in the analysis
    integer :: levels ! JMAX ! ** ! number of levels in the header of the atmospheric profile file
    real, allocatable :: heightArray(:), pressureArray(:), temperatureArray(:), densityArray(:) ! NEW ! ** ! PPP, TTT, RORO
    integer, parameter :: levelsThreshold = 200 ! NEW ! ** !
    character(len=7) :: inputMolecule
    integer :: levelsIdx ! JJJ ! ** ! index for loop over levels
    ! ----------------------------------------------------------- !
    ! real :: pressure ! [atm] -- pressure on the current atmospheric level
    real :: pSelf ! [atm] -- partial pressure of the considered gaseous species
    real :: pForeign ! [atm] -- foreign pressure ( p - pSelf)
    real :: unitlessT
    ! real :: temperature ! [K] -- temperauture on the current atmospheric level
    ! real :: density ! [molecule/(cm^2 * km)] -- such density units are needed for having absorption coefficient in km-1

    ! TODO: consider moving loop over atmospheric levels here
    real :: temperature, pressure, density ! physical parameters of the atmosphere on the current height level
contains
    subroutine readAtmosphericParameters()
        ! ---------- reading from the ATM file ------------------ !
        open(atmProfileUnit, file='data/Atmospheres/'//atmProfileFile, status='old')
        read(atmProfileUnit, '(A20)') atmTitle
        read(atmProfileUnit, *) numGasSpecies, levels
        if (levels > 200) then
            write(*, '(A, I3, A)') 'WARNING: input number of atmospheric levels is &
                                    bigger than ', levelsThreshold, '. WARNING message here.'
            stop
        endif
        read(atmProfileUnit, '(A7)') inputMolecule
        allocate(heightArray(levels))
        allocate(pressureArray(levels))
        allocate(temperatureArray(levels))
        allocate(densityArray(levels))
        do levelsIdx = 1, levels
            read(atmProfileUnit, *) heightArray(levelsIdx), pressureArray(levelsIdx), &
                                    temperatureArray(levelsIdx), densityArray(levelsIdx)
        end do
        close(atmProfileUnit)
    ! ------------------------------------------------------------- !
    end subroutine readAtmosphericParameters
end module Atmosphere