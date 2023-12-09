program main
  use globals
  use shared_vars_main
  use MESH1

  implicit none

  integer, parameter :: IOUT = 47
  real(kind=dp), parameter :: DIAP = 10.0

  integer :: NGS, JMAX, J
  real :: ZZZ
  real(kind=DP) :: VSTARTT, V_END, DLT8

  
  real(kind=DP) :: PPP(200), TTT(200), RORO(200)
  
  character(len=20) :: ATM
  character(len=5) :: N_HAUS, MOLECULE
  character(len=3) :: MOL3, JNAMB
  ! character(len=20) :: LINE_PATH
  ! character(len=50) :: FI

  EPS = 0.0
  ! ******************** !
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
  ! ******************** !

  open(675, file='data/input.txt', status='old')
  read(675, *) VSTART, V_END
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
  WRITE(*,*) N_HAUS
  ! <<<<<<<<<<<<<<<<<< END OF THE USER's PART >>>>>>>>>>>>>>>>>>> !
end program main
