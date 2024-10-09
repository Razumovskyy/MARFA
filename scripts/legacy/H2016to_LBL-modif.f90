! *** 29 June, 2018 ***
!   HITRAN-2016 => the first 7 gases, two digits for the isotop number
!  Numbers of isotops for the first 12 moleculecules to use TIPS-2017 
! Direct ACCESS files 'H16.01','H16.02',...'H16.07'
! from ASCII file 'H16par. 1','H16par. 2', ... ,'H16par. 7'
!  "Species by Species" 
! ************************************************************************

!   Each item is defined below, with its format shown in parenthesis.
! ------------------------------------------------------------------------------------------
! I-----------------------------------------------------------------------------------------
!      1        (I2)  ###        MO = molecule number {I use MO_ISO= MO*10+ISO (I3)}
!      1### now (A1) instead of  (I1)  ISO = isotope number {1 = most abundant, 2 = second, etc.}
!      2        (F12.6)       V = frequency of transition in wavenumbers (cm-1)
!      3        (E10.3)       S = intensity in cm-1/(molec * cm-2) at 296 Kelvin
!*  not used    (E10.3)       A = Einstein A coefficient
!      4        (F5.4)        AGAM = air-broadened halfwidth  in cm-1/atm at 296 Kelvin
!      5        (F5.4)        SGAM = self-broadened halfwidth in cm-1/atm at 296 Kelvin
!      6        (F10.4)       E = lower state energy in wavenumbers (cm-1)
!      7        (F4.2)        N = coeff. of temperat. dependence of air-broadened halfwidth
!      8        (F8.6)        d = shift of transition due to pressure (cm-1)
!*  not used    (A15)         V1 = upper state global quanta 
!*  not used    (A15)         V2 = lower state global quanta 
!*  not used    (A15)         Q1 = upper state local quanta
!*  not used    (A15)         Q2 = lower state local quanta
!*  not used    (6I1)         IERR = accuracy indexes
!*  not used    (6I2)         IREF = reference indexes
!*  not used    (A1)          FLAG = flag for line coupling {(*) if "Yes"}
!*  not used    (F7.1)        Gu = upper stat. weight
!*  not used    (F7.1)        Gl = lower stat. weight
! --------------------------------------------------------------------------------------------

!   The molecule numbers are encoded as shown in the table below:
! __________________________________________________________________
!  1=   H2O    2=  CO2    3=    O3    4=   N2O    5=  CO
!  6=   CH4    7=   O2    8=    NO    9=   SO2   10= NO2   11=    NH3
! 12=  HNO3   13=   OH   14=    HF   15=   HCl   16= HBr   17=     HI
! 18=   ClO   19=  OCS   20=  H2CO   21=  HOCl   22=  N2   23=    HCN
! 24= CH3Cl   25= H2O2   26=  C2H2   27=  C2H6   28= PH3   29=   COF2
! 30=   SF6   31=  H2S   32= HCOOH   33=   H02   34=   O   35= ClONO2
! 36=   NO+   37= HOBr   38=  C2H4   39= CH3OH 
! 2008 ->     40=CH3Br   41= CH3CN   42= CF4
! 2012 ->     43=C4H2    44=  HC3N   45=H2       46=CS     47= SO3 
! ___________________________________________________________________
PROGRAM legacyParConvertion       
       REAL*4 N
       REAL*8 V,S,E
       CHARACTER FNAMW*4,FNAMR*7, XXX*2, ISO*1, YYY*2
	   CHARACTER V1*15, V2*15, Q1*15, Q2*15, FLAG*1
       DIMENSION ISOTOP(18),IERR(6),IREF(6)  
    DIMENSION N_MOL_ISO(20,12)  ! <--- to use TIPS-2017 
 DATA N_MOL_ISO/1,2,3,4,5,6,7,8,9,11*0, & ! 1 H2O
 10,11,12,13,14,15,16,17,18,19,20,21,22,7*0, & ! 2 CO2
 23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,2*0,  & ! 3 O3
 41,42,43,44,45,15*0, & ! 4 N2O
 46,47,48,49,50,51,52,53,54,11*0, &  ! 5 CO
 55,56,57,58,16*0,               &  ! 6 CH4
 59,60,61,62,63,64,14*0, &  ! 7 O2
 65,66,67,17*0, & ! 8 NO
 68,69,18*0, & ! 9  SO2
 70,19*0,    & ! 10  NO2
 71,72,18*0, & ! 11  NH3
 73,74,18*0/ ! 12 HNO3


! ***  Specie by Specie *** !  <---- "by hand" 
       DATA FNAMW/'H16.'/
       DATA FNAMR/'H16par.'/ 

! ***********************************************************************
XXX=' 4'   ! Number of Specie  (H2O -> ' 1', CO2 -> ' 2', etc.)
YYY='04'   ! The same for DIRECT Files 

!XXX='12'   ! Number of Specie  (H2O -> ' 1', CO2 -> ' 2', etc.)
!YYY='12'   ! The same for DIRECT Files 


           OPEN(4,FILE='Report_H2016.'//XXX)
       WRITE(*,*)FNAMR//XXX

	 write(*,*)' RECL=9 - Attention!'
	OPEN(3,FILE=FNAMR//XXX)
OPEN(2, ACCESS='DIRECT', FORM='UNFORMATTED',RECL = 9,FILE=FNAMW//YYY) ! RECL=36 for other FORTRANS
 open(222,FILE=FNAMW//XXX//'ascii')
        ISOTOP=0
      DO NLINE =1, 10000000
!	write(*,*)nline
      READ(3,800,END=10) MO,ISO,V,S,A,AGAM,SGAM,E,N,D &
	  ,V1,V2,Q1,Q2,IERR,IREF,FLAG,Gu,Gl

                     SS=S  ! odinary precision is enough
                     EE=E  ! odinary precision is enough
IF(ISO=='1') IS=1 ; IF(ISO=='2') IS=2 ; IF(ISO=='3') IS=3 
IF(ISO=='4') IS=4 ; IF(ISO=='5') IS=5 ; IF(ISO=='6') IS=6 
IF(ISO=='7') IS=7 ; IF(ISO=='8') IS=8 ; IF(ISO=='9') IS=9 
IF(ISO=='0') IS=10 
IF(ISO=='A') IS=11 ; IF(ISO=='B') IS=12 ; IF(ISO=='C') IS=13 
IF(ISO=='D') IS=14 ; IF(ISO=='E') IS=15 ; IF(ISO=='F') IS=16 
IF(ISO=='G') IS=17 ; IF(ISO=='H') IS=18 ; IF(ISO=='I') IS=19

N_MOLIS=N_MOL_ISO(IS,MO) ! The isotop number in the internal database 
                         ! (to define mol. weight and TIPS-17).

MO_ISO=N_MOLIS*100+MO    ! [isotop number*100 + molecular number]
ISOTOP(IS)=ISOTOP(IS)+1
       WRITE(2,REC=NLINE)  V,SS,AGAM,SGAM,EE,N,MO_ISO,D
             WRITE(222,*)  V,SS,AGAM,SGAM,EE,N,MO_ISO,D
 END DO ! Line-loop
 800  FORMAT( I2,A1,F12.6,2E10.3,2F5.4,F10.4,F4.2,F8.6,4A15,6I1,6I2,A1,F7.1,F7.2)
 10      NLMAX=NLINE-1
        WRITE(4,*) ' Lines = ',nlmax,'   ',fnamw//XXX,fnamr//XXX
        WRITE(4,*)ISOTOP
                 CLOSE(3)
                 CLOSE(2)
      WRITE(*,*)'* This HITRAN-2016 file is ready for LBL calculations*'
             CLOSE(4)

END PROGRAM

