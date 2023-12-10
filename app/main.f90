program main
    !*******************************************************************************
    !*  Principal feature: 9  internal grids                                       *
    !* --------------------------------------------------------------------------- * 
    !*  Calculates PT-TABLE for Venus Atmosphere (up to 1000 K!)                   *
    !* A.  HITRAN- ONLY 1 gas (H2O or CO2  or SO2);                                * 
    !* B. NO Continuum models                                                      *
    !!*******************************************************************************
    use globals
    use shared_vars_main
    use MESH1
    use K_HITRAN

    implicit none

    ! Parameters
    integer, parameter :: IOUT = 47
    real(kind=dp), parameter :: DIAP = 10.0

    ! Numeric Variables
    integer :: NGS, JMAX, J, JJJ
    integer :: NW ! record number in the derect access file
    real :: ZZZ
    real(kind=DP) :: VSTARTT, V_END, DLT8

    ! Arrays (pressure, temperature, density)
    real(kind=DP) :: PPP(200), TTT(200), RORO(200)

    ! String Variables
    character(len=20) :: ATM
    character(len=5) :: N_HAUS
    character(len=5) :: MOLECULE
    character(len=3) :: MOL3
    character(len=3) :: JNAMB ! unique identifier for each Zj level: ('1__','2__',...,'100',...)
    ! character(len=20) :: LINE_PATH
    ! character(len=50) :: FI

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

    open(675, file='data/input.txt', status='old')
    read(675, *) VSTARTT, V_END
    read(675, '(A20)') ATM
    close(675)

    open(676, file='data/Atmospheres/'//ATM, status='old')
    read(676, '(A5)') N_HAUS !!! change the length value (first line of atmospheric profile)
    read(676, *) NGS, JMAX
    if (JMAX > 200) then
        write(*,*) ' ***  JMAX>200  *** '
    endif
    read(676, '(A5)') MOLECULE
    do J = 1, JMAX
        read(676, *) ZZZ, PPP(J), TTT(J), RORO(J)
    end do
    close(676)

    MOL3 = MOLECULE
    ! <<<< End of the user's part >>>> !

    DLT8 = 10.0
    
    open(5555, file='control/PT-Protocol')

    do JJJ = 1, JMAX
        write(*,*) JJJ, JMAX ! for real time tracking how many levels has been processed
        P = PPP(JJJ)
        T = TTT(JJJ)
        RO = RORO(JJJ)
        VSTART = VSTARTT - DLT8
        write(5555, *) JJJ, P, T, JMAX

        JNAMB = '___'
        if ( JJJ < 10 ) then
            write(JNAMB(1:1), '(I1)') JJJ
        else
            if ( JJJ < 100 ) then
                write(JNAMB(1:2), '(I2)') JJJ
            else
                write(JNAMB(1:3), '(I2)') JJJ
            end if
        end if
        
        open(IOUT, access='DIRECT', form='UNFORMATTED', recl=NT*4, &
            file='output/PT_CALC/'//JNAMB//'.'//MOL3)
        ! RECL = NT for Windows Fortrans !

        do J = 0, 99999999
            ! *** calculation inside 10.0 cm^-1 *** !
            VSTART = VSTART + DLT8
            VFINISH = VSTART + DLT8
            
            if ( VSTART >= V_END ) exit

            NW = (VSTART + 1.0) / 10.0 ! *** (0.0 -> 0 , 10.0 -> 1,..., 560.0 -> 56, etc.)
            ! write(*,*) VSTART, V_END

            WVA = VSTART - OBR25
            WVB = VSTART + DLT8 + OBR25

            if ( WVA < 10.0 ) then
                WVA = 10.0
            end if
            
            call K_HITRAN_3g(MOLECULE, JJJ)

            write(IOUT, rec=NW) RK
            ! *** end of the calculation inside 10.0 cm^-1 *** !
        end do
        close(IOUT)
        ! *** end of the calculation over Z *** !
    end do
    close(5555)
    write(*,*) ' *** Congratulations! PT-table is READY! ***'
end program main
