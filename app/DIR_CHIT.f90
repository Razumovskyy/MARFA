PROGRAM DIR_CHIT
implicit none
integer, PARAMETER :: NT=20481
integer :: I, NCASE
integer :: JN
integer :: M, N, NZ
real :: H
CHARACTER MET*3,DIR_NAME*4
real ::  RK(NT)
real :: V, VV

!!! INPUT EXTENSION !!!
MET='H2O'


H=10.0/(NT-1)
! -------------------------------!
DO NCASE=1,5
WRITE(*,*)'  V   J  (0.0 or 0 => STOP) '
READ(*,*)V,JN
DIR_NAME='___.'
IF(V==0.0.OR.JN==0)STOP 
IF(JN<10)THEN
WRITE(DIR_NAME(1:1),'(I1)')JN   
ELSE
     IF(JN<100)THEN
WRITE(DIR_NAME(1:2),'(I2)')JN
	 ELSE
WRITE(DIR_NAME(1:3),'(I3)')JN
	 END IF
END IF
NZ = (V+1.)/10.
OPEN(491,ACCESS='DIRECT', FORM='UNFORMATTED',	&
 RECL = NT*4, FILE='./output/PT_CALC/'//DIR_NAME//MET) ! NT*4 for UNIX
!WRITE(*,*) './PT_CALC/'//DIR_NAME//MET
!PAUSE 
READ(491,REC=NZ)RK
CLOSE(491)
OPEN(88,FILE='SPECTR')
	N=NT ! 10000
	M=1
	WRITE(88,*)M,N
             DO I=1,N 
VV=V+H*(I-1)
WRITE(88,*)VV,ALOG10(RK(I))
             END DO
CLOSE(88)
END DO
END PROGRAM DIR_CHIT