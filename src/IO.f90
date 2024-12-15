module IO
    ! TODO:(!) get rid of this module, it might be redundant. It doesn't perform output operation.
    ! or on the contrary, move output operations here and make it more reasonable. TO REFACTOR.
    use Constants
    use Interfaces
    use Atmosphere
    use Shapes
    use Spectroscopy
    use Grids
    use ChiFactors, only: chiFactorFuncName
    implicit none
    
    ! pointer to the specific line shape function which will be used in grid calculations 
    procedure(shape), pointer :: shapeFuncPtr
    
    character(len=30) :: lineShapeFuncName ! name of the line shape function from the Shapes.f90 module

    character(len=300) :: databaseFile ! file name with line spectral data

    ! Variables and parameters related to TIPS:
    character(len=300), parameter :: TIPSFile = 'data/QofT_formatted.dat' ! path to the file with TIPS data
    integer, parameter :: TIPSUnit = 5467 ! unit for file with TIPS data
    integer :: nIsotopes, nTemperatures ! number of different isotopes and temperatures in the TIPS file
    integer :: stSumTIdx, stSumIsoIdx  ! loop indices for partition sums: temperatures and isotopes

    real(kind=DP) :: startWV, endWV ! [cm-1] -- boundaries of an initial spectral interval [startWV; EndWV]
    real(kind=DP) :: extStartWV, extEndWV ! boundaries of an extended initial interval: [startWV-cutOff; endWV+cutOff]
    real(kind=DP) :: startDeltaWV, endDeltaWV ! boundaries of a subinterval ! VSTART (also: VS,VR4), VFINISH (legacy) !
   
    ! boundaries of an extended subinterval (cutOffs included): [startDeltaWV-cutOff; endWV+cutOff] ! WVA, WVB (legacy) !
    real(kind=DP) :: extStartDeltaWV ! legacy: VA
    real(kind=DP) :: extEndDeltaWV ! legacy: VFISH

contains

    subroutine readTIPS()
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

end module IO
