!*******************************************************************************
! main.f90
!*******************************************************************************
!
! Legacy code: Boris Fomin
! Rework and new features: Mikhail Razumovskii
!
! Program Description:
! This program calculates PT-TABLE (Pressure-Temperature Table) for Venus Atmosphere,
! specifically designed to handle temperatures up to 1000 K. The current version focuses
! on line-by-line spectroscopic analysis for a single gas species at a time,
! supporting H2O, CO2, or SO2. It does not incorporate continuum models.
!
! Principal Features:
! 1. 9 internal grids for detailed spectral analysis.
! 2. Adapted for high-resolution spectroscopic data from HITRAN database.
!
!*******************************************************************************

program main
    use kinds, only: DP
    use shared_vars_main
    use mesh
    use K_HITRAN
    implicit none

    ! MAIN OUTPUT FILE !
    integer, parameter :: outputUnit = 47 ! IOUT ! ** ! output file unit where spectral PT-tables will be stored for each atmospheric level
    integer :: outputRecNum ! NW ! ** ! record number in the output file

    ! INPUT FILES !
    integer, parameter :: configUnit = 675 ! NEW ! ** !
    integer, parameter :: atmProfileUnit = 676 ! NEW ! ** ! 
    integer, parameter :: stSumUnit = 677 ! NEW ! ** !

    ! CONTROL AND DEBUG FILES !
    integer, parameter :: atmControlUnit = 5555 ! NEW ! ** !

    ! LIMITATIONS !
    integer, parameter :: levelsThreshold = 200 ! NEW ! ** !

    ! Input config file !
    !real(kind=DP) :: startWV, endWV
    character(len=20) :: atmProfileFile ! NEW ! ** ! ATM
    !real(kind=DP) :: deltaWV ! NEW ! ** ! DLT8 ! also declared in LBL2023 !!! ! to be added to the input config file

    ! Atmospheric Profile Data !
    character(len=20) :: atmTitle
    integer :: numGasSpecies ! NGS ! ** ! number of species in the analysis
    integer :: levels ! JMAX ! ** ! number of levels in the header of the atmospheric profile file
    character(len=7) :: inputMolecule
    real, allocatable :: height(:), pressure(:), temperature(:), density(:) ! NEW ! ** ! PPP, TTT, RORO
    integer :: levelsIdx ! JJJ ! ** ! index for loop over levels

    ! Statistical sums file !
    integer :: nIsotopes, nTemperatures
    integer :: stSumTIdx, stSumIsoIdx  ! indices for statistical sums

    ! Other variables !
    character(len=3) :: reducedMoleculeName ! 
    character(len=3) :: levelLabel ! unique identifier for each Zj level: ('1__','2__',...,'100',...)


    ! TODO: optimization proposition. Reduce reading operations because of highly overlaping intervals: [extStartWV1, extEndWV1] and [extStartWV2, extEndWV2]
    
    EPS = 0.0
    H0=STEP
    H1=H0/2.
    H2=H1/2.
    H3=H2/2.
    H4=H3/2.
    H5=H4/2.
    H6=H5/2.
    H7=H6/2. 
    H8=H7/2. 
    H9=H8/2. 
    H=H9/4.

    ! --------- reading from the config file -------------- !
    open(configUnit, file='simConfig.ini', status='old')
    read(configUnit, *) startWV, endWV
    read(configUnit, '(A20)') atmProfileFile
    close(configUnit)
    ! ------------------------------------------------------ !

    ! ---------- reading from the ATM file ------------------ !

    open(atmProfileUnit, file='data/Atmospheres/'//atmProfileFile, status='old')
    read(atmProfileUnit, '(A20)') atmTitle
    read(atmProfileUnit, *) numGasSpecies, levels
    if (levels > 200) then
        write(*, '(A, I3, A)') 'WARNING: input number of atmospheric levels is &
                                bigger than ', levelsThreshold, '. WARNING message here.'
    endif
    read(atmProfileUnit, '(A7)') inputMolecule
    allocate(height(levels))
    allocate(pressure(levels))
    allocate(temperature(levels))
    allocate(density(levels))
    do levelsIdx = 1, levels
        read(atmProfileUnit, *) height(levelsIdx), pressure(levelsIdx), temperature(levelsIdx), density(levelsIdx)
    end do
    close(atmProfileUnit)
    ! -------------------------------------------------------- !

    ! ---------- reading Statistical Sums --------------------- !
    open(unit=stSumUnit, file='data/QofT_formatted.dat', status='old', action='read')
    read(stSumUnit, *) nIsotopes, nTemperatures
    allocate(QofT(nIsotopes, nTemperatures))
    do stSumIsoIdx = 1, nIsotopes
        read(stSumUnit, *) (QofT(stSumIsoIdx, stSumTIdx), stSumTIdx=1, nTemperatures)
    end do
    close(stSumUnit)
    ! --------------------------------------------------------- !

    reducedMoleculeName = inputMolecule
    
    open(atmControlUnit, file='control/PT-Protocol')

    do levelsIdx = 1, levels
        write(*,*) levelsIdx, levels ! for real time tracking how many levels has been processed
        P = pressure(levelsIdx)
        T = temperature(levelsIdx)
        RO = density(levelsIdx)
        write(atmControlUnit, *) levelsIdx, P, T, levels

        levelLabel = '___'
        if ( levelsIdx < 10 ) then
            write(levelLabel(1:1), '(I1)') levelsIdx
        else
            if ( levelsIdx < 100 ) then
                write(levelLabel(1:2), '(I2)') levelsIdx
            else
                write(levelLabel(1:3), '(I2)') levelsIdx
            end if
        end if
        
        open(outputUnit, access='DIRECT', form='UNFORMATTED', recl=NT*4, &
            file='output/PT_CALC/'//levelLabel//'.'//reducedMoleculeName)
        ! RECL = NT for Windows Fortrans !

        
        startDeltaWV = startWV 
        endDeltaWV = startDeltaWV + deltaWV
        do while (startDeltaWV < endWV)
            ! *** calculation inside 10.0 cm^-1 *** !

            outputRecNum = (startDeltaWV + 1.0) / 10.0 ! *** (0.0 -> 0 , 10.0 -> 1,..., 560.0 -> 56, etc.)
            ! write(*,*) startDeltaWV, V_END

            ! if ( extStartDeltaWV < startWV ) then
            !     extStartDeltaWV = startDeltaWV
            ! end if

            ! write(*,*) 'startDeltaWV before K_HITRAN call: ', startDeltaWV
            ! pause
            call K_HITRAN_3g(inputMolecule, levelsIdx)

            write(outputUnit, rec=outputRecNum) RK

            startDeltaWV = startDeltaWV + deltaWV
            endDeltaWV = startDeltaWV + deltaWV
            
            ! *** end of the calculation inside 10.0 cm^-1 *** !
        end do
        close(outputUnit)
        ! *** end of the calculation over Z *** !
    end do
    close(atmControlUnit)
    write(*,*) ' *** Congratulations! PT-table is READY! ***'
    deallocate(height)
    deallocate(pressure)
    deallocate(density)
    deallocate(temperature)
    deallocate(QofT)
end program main
