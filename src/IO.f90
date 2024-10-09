module IO
    use Constants
    use ShapeFuncInterface
    use Atmosphere
    use Shapes
    use Spectroscopy
    use Mesh
    use ChiFactors, only: chiFactorFuncName
    implicit none
    integer, parameter :: TIPSUnit = 5467
    character(len=30), parameter :: TIPSFile = 'data/QofT_formatted.dat'
    integer, parameter :: configUnit = 555
    integer, parameter :: stSumUnit = 666
    real(kind=DP) :: startWV, endWV, calcPrecision ! [cm-1] -- boundaries and grid calcPrecision for the spectral interval for calculating spectra 
    character(len=30) :: lineShapeFuncName ! name of the line shape function (custom or standard)

    procedure(shape), pointer :: shapeFuncPtr ! pointer for implementing different line shapes functions

    integer :: totalLines ! number of lines in one HITRAN file (e.g. H16.01) 
    character(len=100) :: hitranFile ! file names with HITRAN raw data (one file per molecule)
    
    ! Input config file !
    !real(kind=DP) :: startWV, endWV
    !real(kind=DP) :: deltaWV ! NEW ! ** ! DLT8 ! also declared in LBL2023 !!! ! to be added to the input config file

    ! Statistical sums file !
    integer :: nIsotopes, nTemperatures
    integer :: stSumTIdx, stSumIsoIdx  ! indices for statistical sums

    ! Shape Vars
    real(kind=DP) :: VI
    real :: SLL
    real(kind=DP) :: LorHWHM ! AL ! Lorentz HWHM -- T and P dependent !
    real(kind=DP) :: DopHWHM ! ADD ! Doppler HWHM

    ! SharedVars Main
    real(kind=DP) :: extStartWV, extEndWV
    real(kind=DP) :: startDeltaWV, endDeltaWV ! starting and ending wavenumbers inside deltaWV interval *** ! VSTART, VFINISH !
    real(kind=DP) :: extStartDeltaWV, extEndDeltaWV ! extended deltaWV intervals including cut-offs ! WVA, WVB !

    real(kind=DP) :: EPS ! ??? epsilon ???
contains
    subroutine readInputParameters()
    ! --------- reading from the config file -------------- !
    open(configUnit, file='simConfig.ini', status='old')
    read(configUnit, '(A7)') inputMolecule
    read(configUnit, *) startWV, endWV
    read(configUnit, '(A20)') atmProfileFile
    read(configUnit, *) cutOff
    read(configUnit, '(A30)') chiFactorFuncName
    close(configUnit)
    ! ------------------------------------------------------ !
    end subroutine readInputParameters


    subroutine readTIPS()
    ! ---------- reading Statistical Sums --------------------- !
    open(unit=TIPSUnit, file='data/TIPS/TIPS.dat', status='old', action='read')
    read(TIPSUnit, *) nIsotopes, nTemperatures
    allocate(TIPS(nIsotopes, nTemperatures))
    do stSumIsoIdx = 1, nIsotopes
        read(TIPSUnit, *) (TIPS(stSumIsoIdx, stSumTIdx), stSumTIdx=1, nTemperatures)
    end do
    close(TIPSUnit)
    ! --------------------------------------------------------- !
    end subroutine readTIPS

    subroutine generateOutput()

    end subroutine generateOutput
end module IO
