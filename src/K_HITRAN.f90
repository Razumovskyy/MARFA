module K_HITRAN
    use globals
    use shared_vars_main
    use MESH1
    use molecule_vars
    implicit none

contains
    subroutine K_HITRAN_3g(MOLECULE, JO)
        character(len=5) :: MOLECULE
        integer :: JO
        
        ! LEGACY -> MODERN tips
        ! SAVE blocks -- consider to remove completely
        ! DATA blocks used for the initial declaration of the variable. 
        ! instead use inline initialization in type declaration.

        !SAVE ISTART,LINBEG,LINBEG0,MOTYPE,NLIN,LINE_PATH,JOLD,VS,VFISH,VA,VAA,VAA0    

        ! DATA ISTART/0/
        ! if (ISTART == 0) then
        !     JOLD = 0
        !     ISTART = 1
        ! end if

        select case (MOLECULE)
        case ('H2O')
            MOTYPE = 1
            LINE_PATH = './HITRAN16/H16.01'
            NLIN = 3892633
        case ('CO2')
            MOTYPE = 2
            LINE_PATH = './HITRAN16/H16.02'
            NLIN = 10868410
        case ('SO2')
            MOTYPE = 1
            LINE_PATH = './HITRAN16/H16.09'
            NLIN = 95120
        case default
            write(*,*) 'Molecule = ', MOLECULE
            ! stopping the calculation because of the unknown type of the molecule
            stop
        end select ! LINE 34 in h2.f90, legacy code
    end subroutine
end module
