! legacy fortran script for plotting and data processing !
program chit
implicit none
integer, parameter :: NT=20481
integer :: I, NCASE
integer :: JN
integer :: M, N, NZ
real :: H
character MET*3, DIR_NAME*4
real ::  RK(NT)
real :: V, VV

MET='CO2' ! input extension

H = 10.0 / (NT-1)
! -------------------------------!
do NCASE = 1, 5
	write(*,*)'  V   J  (0.0 or 0 => STOP) '
	read(*,*) V, JN
	
	DIR_NAME='___.'
	
	if (V==0.0 .or. JN==0) STOP 
	
	if (JN < 10) then
	write(DIR_NAME(1:1),'(I1)') JN   
	else
		if (JN < 100) then
			write(DIR_NAME(1:2), '(I2)') JN
		else
			write(DIR_NAME(1:3), '(I3)') JN
		end if
	end if

	NZ = (V+1.)/10.
	
	open(491, access='DIRECT', form='UNFORMATTED',	&
 		recl = NT*4, file='./output/PT_CALC/'//DIR_NAME//MET) ! NT*4 for UNIX

	read(491, rec=NZ) RK ! reads only one record with deltaWV=10 cm-1 interval data
	close(491)
	
	open(88, file='SPECTR')
	N = NT ! 10000
	M = 1
	write(88,*) M, N
            do I = 1, N 		
				VV = V + H*(I-1)
				write(88,*) VV, ALOG10(RK(I))
			end do
			close(88)
		end do
end program chit
