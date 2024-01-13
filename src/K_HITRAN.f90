module K_HITRAN
    use kinds
    use shared_vars_main
    use mesh
    use molecule_vars
    use LBL
    use van_fleck_huber_vars
    implicit none

contains
    subroutine K_HITRAN_3g(molecule, loopLevel)
        character(len=5), intent(in) :: molecule
        integer, intent(in) :: loopLevel
        
        integer, parameter :: hitranFileUnit = 7777 ! Unit of file with HITRAN data

        logical, save :: isFirstCall = .true.
        real(kind=DP), save :: startingLineWV ! VAA0
        integer, save :: currentLevel ! JOLD
        integer, save :: startingLineIdx ! LINBEG0

        ! TODO understand better why `save` is needed for these guys:
        real(kind=DP), save :: lineWV ! VAA ! spectral position of the first line from HITRAN file which is in the loop interval
        integer, save :: lineIdx ! LINBEG
        
        ! real(kind=DP) :: VFISH ! *** ! extEndDeltaWV
        ! real(kind=DP) :: VA ! *** ! extStartDeltaWV
        ! real(kind=DP) :: VS ! the same as startDeltaWV
        ! real(kind=DP) :: VR4 ! the same as startDeltaWV

        integer :: I, J, M ! loop variables in lines 122 - 231 
        
        ! TODO: moving `ifFirstCall` part to the separate subroutine, because
        ! it is only for reading from the HITRAN file and must be done only at once.

        ! TODO: deal with `totalLines` (now it is from HITEMP, but it seems that any large number is sufficient)

        ! TODO: deal with `molType`
        
        ! write(*,*) 'loopLevel: ', loopLevel
        ! write(*,*) 'startDeltaWV: ', startDeltaWV
        ! pause

        ! The first call operations:
        ! 1. Fetch the molecule name and its molType
        ! 2. Fetch corresponding HITRAN file and a total number of lines in this file
        ! 3. Interrupt if unsuported molecule
        ! 4. 
        if (isFirstCall) then
            isFirstCall = .false.
            startingLineIdx = 1
            currentLevel = 0
            select case (molecule)
            case ('H2O')
                molType = 1
                hitranFile = 'data/HITRAN16/H16.01'
                totalLines = 3892633
            case ('CO2')
                molType = 2
                hitranFile = 'data/HITRAN16/H16.02'
                totalLines = 10868410
            case ('SO2')
                molType = 1
                hitranFile = 'data/HITRAN16/H16.09'
                totalLines = 95120
            case default
                write(*,*) 'Molecule = ', molecule
                write(*,*) 'Unsupported type of the molecule !'
                ! stopping the calculation because of the unknown type of the molecule
                stop
            end select

            open(hitranFileUnit, access='DIRECT', form='UNFORMATTED', recl=36, file=hitranFile)
            ! Uncomment if the direct access file was created under Windows
            ! open(7777, access='DIRECT', form='UNFORMATTED', recl=9, file=hitranFile)

            read(hitranFileUnit, rec=startingLineIdx) startingLineWV
            if (startDeltaWV > cutOff) then
                extStartDeltaWV = startDeltaWV - cutOff
                do while (startingLineWV <= extStartDeltaWV)
                    startingLineIdx = startingLineIdx + 1
                    read(hitranFileUnit, rec=startingLineIdx) startingLineWV
                end do
            end if
        end if

        if ( currentLevel /= loopLevel ) then
            currentLevel = loopLevel
            lineIdx = startingLineIdx
            lineWV = startingLineWV
        end if

        ! -----------------UNKNOWN FEATURE --------------------!
        lineWV = startDeltaWV + DELTA - cutOff
        if (lineWV <= startingLineWV) lineWV = startingLineWV
        ! ------------------------------------------------------!
        
        extEndDeltaWV = startDeltaWV + delta + cutOff

        ! VS = startDeltaWV
        !*-------------------------------------------------------------

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

        !! Note: Arrays are initialized to zero at the beginning of each subroutine call. 
        ! This approach ensures that the arrays have a clean state for each calculation. 
        ! If the subroutine is called multiple times and the arrays do not need to be 
        ! reset every time, consider optimizing by removing this initialization or 
        ! modifying it as per the specific requirements of the program.

        !*--------------------------------------------------------------
        
        ! TODO: remove lineWV from LBL2023 input parameters, because it is accessible through the `shared_vars_main` module
        
        ! write(*,*) 'lineIdx before LBL2023: ', lineIdx
        ! pause
        call LBL2023(molType, lineIdx, lineWV, totalLines)
        ! write(*,*) 'lineIdx after LBL2023: ', lineIdx
        ! pause
        ! ------------ END OF THE FIRST PART ---------------------------------- !
        
        ! ----------------------SECOND PART ----------------------------------- !
        ! TODO: CREATE SUBROUTINE TO REMOVE REPETITIONS !
        ! be careful that last loop differs from NT0-NT8!
        
        do J = 1, NT0
            I = J * 2 - 1
            RK1P(I) = RK1P(I) + RK0P(J)
            RK1(I) = RK1(I) + RK0P(J) * 0.375 + RK0(J) * 0.75 - RK0L(J) * 0.125
            RK1L(I) = RK1L(I) + RK0(J)
            M = I + 1
            RK1P(M) = RK1P(M) + RK0(J)
            RK1(M) = RK1(M) + RK0L(J) * 0.375 + RK0(J) * 0.75 - RK0P(J) * 0.125
            RK1L(M) = RK1L(M) + RK0L(J)
        end do

        do J = 1, NT1
            I = J * 2 - 1
            RK2P(I) = RK2P(I) + RK1P(J)
            RK2(I) = RK2(I) + RK1P(J) * 0.375 + RK1(J) * 0.75 - RK1L(J) * 0.125
            RK2L(I) = RK2L(I) + RK1(J)
            M = I + 1
            RK2P(M) = RK2P(M) + RK1(J)
            RK2(M) = RK2(M) + RK1L(J) * 0.375 + RK1(J) * 0.75 - RK1P(J) * 0.125
            RK2L(M) = RK2L(M) + RK1L(J)
        end do

        do J = 1, NT2
            I = J * 2 - 1
            RK3P(I) = RK3P(I) + RK2P(J)
            RK3(I) = RK3(I) + RK2P(J) * 0.375 + RK2(J) * 0.75 - RK2L(J) * 0.125
            RK3L(I) = RK3L(I) + RK2(J)
            M = I + 1
            RK3P(M) = RK3P(M) + RK2(J)
            RK3(M) = RK3(M) + RK2L(J) * 0.375 + RK2(J) * 0.75 - RK2P(J) * 0.125
            RK3L(M) = RK3L(M) + RK2L(J)
        end do

        do J = 1, NT3
            I = J * 2 - 1
            RK4P(I) = RK4P(I) + RK3P(J)
            RK4(I) = RK4(I) + RK3P(J) * 0.375 + RK3(J) * 0.75 - RK3L(J) * 0.125
            RK4L(I) = RK4L(I) + RK3(J)
            M = I + 1
            RK4P(M) = RK4P(M) + RK3(J)
            RK4(M) = RK4(M) + RK3L(J) * 0.375 + RK3(J) * 0.75 - RK3P(J) * 0.125
            RK4L(M) = RK4L(M) + RK3L(J)
        end do

        do J = 1, NT4
            I = J * 2 - 1
            RK5P(I) = RK5P(I) + RK4P(J)
            RK5(I) = RK5(I) + RK4P(J) * 0.375 + RK4(J) * 0.75 - RK4L(J) * 0.125
            RK5L(I) = RK5L(I) + RK4(J)
            M = I + 1
            RK5P(M) = RK5P(M) + RK4(J)
            RK5(M) = RK5(M) + RK4L(J) * 0.375 + RK4(J) * 0.75 - RK4P(J) * 0.125
            RK5L(M) = RK5L(M) + RK4L(J)
        end do

        do J = 1, NT5
            I = J * 2 - 1
            RK6P(I) = RK6P(I) + RK5P(J)
            RK6(I) = RK6(I) + RK5P(J) * 0.375 + RK5(J) * 0.75 - RK5L(J) * 0.125
            RK6L(I) = RK6L(I) + RK5(J)
            M = I + 1
            RK6P(M) = RK6P(M) + RK5(J)
            RK6(M) = RK6(M) + RK5L(J) * 0.375 + RK5(J) * 0.75 - RK5P(J) * 0.125
            RK6L(M) = RK6L(M) + RK5L(J)
        end do
        
        do J = 1, NT6
            I = J * 2 - 1
            RK7P(I) = RK7P(I) + RK6P(J)
            RK7(I) = RK7(I) + RK6P(J) * 0.375 + RK6(J) * 0.75 - RK6L(J) * 0.125
            RK7L(I) = RK7L(I) + RK6(J)
            M = I + 1
            RK7P(M) = RK7P(M) + RK6(J)
            RK7(M) = RK7(M) + RK6L(J) * 0.375 + RK6(J) * 0.75 - RK6P(J) * 0.125
            RK7L(M) = RK7L(M) + RK6L(J)
        end do

        do J = 1, NT7
            I = J * 2 - 1
            RK8P(I) = RK8P(I) + RK7P(J)
            RK8(I) = RK8(I) + RK7P(J) * 0.375 + RK7(J) * 0.75 - RK7L(J) * 0.125
            RK8L(I) = RK8L(I) + RK7(J)
            M = I + 1
            RK8P(M) = RK8P(M) + RK7(J)
            RK8(M) = RK8(M) + RK7L(J) * 0.375 + RK7(J) * 0.75 - RK7P(J) * 0.125
            RK8L(M) = RK8L(M) + RK7L(J)
        end do

        do J = 1, NT8
            I = J * 2 - 1
            RK9P(I) = RK9P(I) + RK8P(J)
            RK9(I) = RK9(I) + RK8P(J) * 0.375 + RK8(J) * 0.75 - RK8L(J) * 0.125
            RK9L(I)= RK9L(I) + RK8(J)
            M = I + 1
            RK9P(M) = RK9P(M) + RK8(J)
            RK9(M) = RK9(M) + RK8L(J) * 0.375 + RK8(J) * 0.75 - RK8P(J) * 0.125
            RK9L(M)= RK9L(M) + RK8L(J)
        end do

        I=1
        do J = 1, NT9
            I = I+1
            RK(I) = RK(I) + (RK9P(J) * 0.375 + RK9(J) * 0.75 - RK9L(J) * 0.125)
            I = I+1
            RK(I) = RK(I) + RK9(J)
            I = I+1
            RK(I) = RK(I) + (RK9L(J) * 0.375 + RK9(J) * 0.75 - RK9P(J) * 0.125)
            I = I + 1
            RK(I) = RK(I) + RK9L(J)
        end do

        ! ------------------------ END OF THE SECOND PART -------------------------- !

        ! ------------------------ THIRD PART ---------------------------------- !
        ! Conditional Applying the Van-Vleck-Weisskopf-Huber factor !
        
        ! VR4 = VS
        
        if ( startDeltaWV >= 2000. ) then
            JM1 = 0
            do J = 1, NT
                FACTV = startDeltaWV + H * JM1
                RK(J) = RK(J) * FACTV  ! <------ key line here
                if ( RK(J) < 0.) RK(J) = 0.
                JM1 = JM1 + 1
            end do
        else
            JM1 = 0
            do J = 1, NT
                VIVI = startDeltaWV + H * JM1
                EVV_ = VIVI * 1.438786 / T
                if (VIVI < 0.1) then
                    FACTV = VIVI * EVV_ / (2. - EVV_)
                else
                    EVV = exp(-EVV_)
                    FACTV = VIVI * (1. - EVV) / (1. + EVV)
                end if
                RK(J) = RK(J) * FACTV ! <------ key line here
                if (RK(J) < 0.) RK(J) = 0.
                JM1 = JM1 + 1
            end do
        end if

        ! ------------------------END OF THE THIRD PART ---------------------------- !
    end subroutine
end module
