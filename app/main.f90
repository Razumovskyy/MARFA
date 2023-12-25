!*******************************************************************************
! main.f90
!*******************************************************************************
!
! Legacy Code Author: Boris Fomin
! Rework and adding new features: Mikhail Razumovskii
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
    use MESH1
    use K_HITRAN
    implicit none

    ! MAIN OUTPUT FILE !
    integer, parameter :: OUTPUT_UNIT = 47 ! IOUT ! ** ! output file unit where spectral PT-tables will be stored for each atmospheric level
    integer :: outputRecNum ! NW ! ** ! record number in the output file

    ! INPUT FILES !
    integer, parameter :: CONFIG_UNIT = 675 ! NEW ! ** !
    integer, parameter :: ATM_PROFILE_UNIT = 676 ! NEW ! ** ! 
    integer, parameter :: ST_SUM_UNIT = 677 ! NEW ! ** !

    ! CONTROL AND DEBUG FILES !
    integer, parameter :: ATM_CONTROL_UNIT = 5555 ! NEW ! ** !

    ! LIMITATIONS !
    integer, parameter :: levelsThreshold = 200 ! NEW ! ** !

    ! Atmospheric Profile Data
    integer :: numGasSpecies ! NGS ! ** ! number of species in the analysis
    integer :: levels ! levels ! ** ! number of levels in the header of the atmospheric profile file
    real :: ZZZ
    real(kind=DP) :: VSTARTT, V_END, DLT8

    ! Indices !
    integer :: levelsIdx ! levelsIdx ! ** ! 

    ! Arrays (pressure, temperature, density)
    real :: PPP(200), TTT(200), RORO(200)

    ! String Variables
    character(len=20) :: ATM
    character(len=5) :: N_HAUS
    character(len=5) :: MOLECULE
    character(len=3) :: MOL3
    character(len=3) :: JNAMB ! unique identifier for each Zj level: ('1__','2__',...,'100',...)
    ! character(len=20) :: LINE_PATH
    ! character(len=50) :: FI

    integer :: kk, ll, nMolecules, nTemperatures

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

    open(CONFIG_UNIT, file='simConfig.ini', status='old')
    read(CONFIG_UNIT, *) VSTARTT, V_END
    read(CONFIG_UNIT, '(A20)') ATM
    close(CONFIG_UNIT)

    open(ATM_PROFILE_UNIT, file='data/Atmospheres/'//ATM, status='old')
    read(ATM_PROFILE_UNIT, '(A5)') N_HAUS !!! change the length value (first line of atmospheric profile)
    read(ATM_PROFILE_UNIT, *) numGasSpecies, levels
    
    if (levels > 200) then
        write(*, '(A, I3, A)') 'WARNING: input number of atmospheric levels is // &
                                bigger than ', levelsThreshold, '. WARNING message here.'
    endif
    
    read(ATM_PROFILE_UNIT, '(A5)') MOLECULE
    do levelsIdx = 1, levels
        read(ATM_PROFILE_UNIT, *) ZZZ, PPP(levelsIdx), TTT(levelsIdx), RORO(levelsIdx)
    end do
    close(ATM_PROFILE_UNIT)

    MOL3 = MOLECULE

    DLT8 = 10.0

    open(unit=ST_SUM_UNIT, file='data/QofT_formatted.dat', status='old', action='read')
    read(ST_SUM_UNIT, *) nMolecules, nTemperatures
    allocate(QofT(nMolecules, nTemperatures))
    do kk = 1, nMolecules
        read(ST_SUM_UNIT, *) (QofT(kk, ll), ll=1, nTemperatures)
    end do
    close(ST_SUM_UNIT)
    
    open(5555, file='control/PT-Protocol')

    do levelsIdx = 1, levels
        write(*,*) levelsIdx, levels ! for real time tracking how many levels has been processed
        P = PPP(levelsIdx)
        T = TTT(levelsIdx)
        RO = RORO(levelsIdx)
        VSTART = VSTARTT - DLT8 ! ???
        write(5555, *) levelsIdx, P, T, levels
        
        !DEBUG SECTION
        !write(*,*) levelsIdx, P, T, levels

        JNAMB = '___'
        if ( levelsIdx < 10 ) then
            write(JNAMB(1:1), '(I1)') levelsIdx
        else
            if ( levelsIdx < 100 ) then
                write(JNAMB(1:2), '(I2)') levelsIdx
            else
                write(JNAMB(1:3), '(I2)') levelsIdx
            end if
        end if
        
        open(OUTPUT_UNIT, access='DIRECT', form='UNFORMATTED', recl=NT*4, &
            file='output/PT_CALC/'//JNAMB//'.'//MOL3)
        ! RECL = NT for Windows Fortrans !

        do while (VSTART < V_END)
            ! *** calculation inside 10.0 cm^-1 *** !
            ! write(*,*) 'VSTART: ', VSTART
            VSTART = VSTART + DLT8
            VFINISH = VSTART + DLT8
            ! write(*,*) 'V_END: ', V_END
            ! write(*,*) VSTART, VFINISH
            ! pause

            outputRecNum = (VSTART + 1.0) / 10.0 ! *** (0.0 -> 0 , 10.0 -> 1,..., 560.0 -> 56, etc.)
            ! write(*,*) VSTART, V_END

            WVA = VSTART - OBR25
            WVB = VSTART + DLT8 + OBR25

            if ( WVA < 10.0 ) then
                WVA = 10.0
            end if
            
            ! DEBUG SECTION
            ! write(*,*) '!---------input to subroutine K_HITRAN_3g-----------!'
            ! write(*,*) 'MOLECULE: ', MOLECULE
            ! write(*,*) 'levelsIdx: ', levelsIdx
            ! write(*,*) '---------end of input to subroutine K_HITRAN_3g------!'

            call K_HITRAN_3g(MOLECULE, levelsIdx)

            write(OUTPUT_UNIT, rec=outputRecNum) RK
            ! *** end of the calculation inside 10.0 cm^-1 *** !
        end do
        close(OUTPUT_UNIT)
        ! *** end of the calculation over Z *** !
    end do
    close(5555)
    write(*,*) ' *** Congratulations! PT-table is READY! ***'
    deallocate(QofT)
end program main
