PROGRAM MY_CHIT
    implicit none
    integer, PARAMETER :: NT=20481
    integer :: I, II, NCASE
    integer :: JN
    integer :: M, N, NZ
    integer :: NZ1, NZ2, count
    real :: H
    CHARACTER MET*3,DIR_NAME*4
    real ::  RK(NT)
    real :: V1, V2, V11, VV

    !!! INPUT EXTENSION !!!
    MET='CO2'

    ! -------------------------------!
    DO NCASE=1,5
        WRITE(*,*)'  V1  V2  J  (0.0 or 0 => STOP) '
        READ(*,*)V1,V2,JN
        DIR_NAME='___.'
        IF(V1==0.0.OR.V2==0.0.OR.JN==0) STOP 
        IF(JN<10)THEN
            WRITE(DIR_NAME(1:1),'(I1)')JN   
        ELSEIF(JN<100)THEN
            WRITE(DIR_NAME(1:2),'(I2)')JN
        ELSE
            WRITE(DIR_NAME(1:3),'(I3)')JN
        END IF
        NZ1 = (V1+1.)/10.
        NZ2 = (V2+1.)/10.
        count = NZ2 - NZ1
        OPEN(491, ACCESS='DIRECT', FORM='UNFORMATTED', &
             RECL = NT*4, FILE='./output/PT_CALC/'//DIR_NAME//MET) ! NT*4 for UNIX
        !WRITE(*,*) './PT_CALC/'//DIR_NAME//MET
        !PAUSE 
        
        OPEN(188, FILE='SPECTR')
            M = 1
        WRITE(188,*) M, NT*count
        
        DO II = 1, count
            READ(491, REC=NZ1+II-1) RK
            V11 = V1 + 10.0 * (II - 1)
            DO I=1,NT
                VV = V11 + (I-1) * (10.0) / (NT - 1)
                WRITE(188,*)VV,ALOG10(RK(I))
            END DO
        END DO
        CLOSE(491)
        CLOSE(188)
    END DO
END PROGRAM MY_CHIT
