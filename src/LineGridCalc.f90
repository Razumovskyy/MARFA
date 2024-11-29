module LineGridCalc
    use Constants
    use Mesh
    use Interfaces
    implicit none
contains
    subroutine leftLBL(FREQ, UL, FSHAPE)
        ! calculation contributions to the subinterval on each grid
        ! FROM the LEFT PART of the extended subinterval: [startDeltaWV-cutOff; startDeltaWV]

        IMPLICIT INTEGER*4 (I-N)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU
        real(kind=DP) :: EPS
        real(kind=DP) :: XXX
        real :: FF
        procedure(shape), pointer :: FSHAPE

        UU=UL-FREQ
                            IF(UU >= 0.)RETURN
        IF(-UU>cutOff)RETURN
        FF=FSHAPE(UU)
            IF(FF < EPS)		RETURN
        RK(1)=FF+RK(1)
            IF(-UU < H0)		GOTO 20
        RK0P(1)=RK0P(1)+FF
        FF=FSHAPE(UU-H1)
        RK0(1)=RK0(1)+FF
        FF=FSHAPE(UU-H0)
        RK0L(1)=RK0L(1)+FF
                                    GOTO 10
     20	IF(-UU < H1)			GOTO 21
        RK1P(1)=RK1P(1)+FF
        FF=FSHAPE(UU-H2)
        RK1(1)=RK1(1)+FF
        FF=FSHAPE(UU-H1)
        RK1L(1)=RK1L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 121
    !* 21 IF(-UU < H2.and.a < h2)goto 22
     21	IF(-UU < H2)			GOTO 22
        RK2P(1)=RK2P(1)+FF
        FF=FSHAPE(UU-H3)
        RK2(1)=RK2(1)+FF
        FF=FSHAPE(UU-H2)
        RK2L(1)=RK2L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 122
     22	IF(-UU < H3)			GOTO 23
        RK3P(1)=RK3P(1)+FF
        FF=FSHAPE(UU-H4)
        RK3(1)=RK3(1)+FF
        FF=FSHAPE(UU-H3)
        RK3L(1)=RK3L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 123
     23	IF(-UU < H4)			GOTO 24
        RK4P(1)=RK4P(1)+FF
        FF=FSHAPE(UU-H5)
        RK4(1)=RK4(1)+FF
        FF=FSHAPE(UU-H4)
        RK4L(1)=RK4L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 124
     24	IF(-UU < H5)			GOTO 25
        RK5P(1)=RK5P(1)+FF
        FF=FSHAPE(UU-H6)
        RK5(1)=RK5(1)+FF
        FF=FSHAPE(UU-H5)
        RK5L(1)=RK5L(1)+FF
        IF(FF < EPS)			RETURN
                                    GOTO 125
     25	IF(-UU < H6)			GOTO 26
        RK6P(1)=RK6P(1)+FF
        FF=FSHAPE(UU-H7)
        RK6(1)=RK6(1)+FF
        FF=FSHAPE(UU-H6)
        RK6L(1)=RK6L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 126
     26	IF(-UU < H7)			GOTO 27
        RK7P(1)=RK7P(1)+FF
        FF=FSHAPE(UU-H8)
        RK7(1)=RK7(1)+FF
        FF=FSHAPE(UU-H7)
        RK7L(1)=RK7L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 127
     27	IF(-UU < H8)			GOTO 28
        RK8P(1)=RK8P(1)+FF
        FF=FSHAPE(UU-H9)
        RK8(1)=RK8(1)+FF
        FF=FSHAPE(UU-H8)
        RK8L(1)=RK8L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 128
     28	IF(-UU < H9)			GOTO 29
        RK9P(1)=RK9P(1)+FF
        FF=FSHAPE(UU-H-H)
        RK9(1)=RK9(1)+FF
        FF=FSHAPE(UU-H9)
        RK9L(1)=RK9L(1)+FF
            IF(FF < EPS)			RETURN
                                    GOTO 129
     29	RK(2)=RK(2)+FSHAPE(UU-H)
        RK(3)=RK(3)+FSHAPE(UU-H-H)
        RK(4)=RK(4)+FSHAPE(UU+H-H9)
        FF=FSHAPE(UU-H9)
    !* CASCADE
        RK(5)=RK(5)+FF
     129 RK9P(2)=RK9P(2)+FF
        FF=FSHAPE(UU-H9-H-H)
        RK9(2)=RK9(2)+FF
        FF=FSHAPE(UU-H8)
        RK9L(2)=RK9L(2)+FF
     128 RK8P(2)=RK8P(2)+FF
        FF=FSHAPE(UU-H8-H9)
        RK8(2)=RK8(2)+FF
        FF=FSHAPE(UU-H7)
        RK8L(2)=RK8L(2)+FF
        IF(FF < EPS)		RETURN
     127 RK7P(2)=RK7P(2)+FF
        FF=FSHAPE(UU-H7-H8)
        RK7(2)=RK7(2)+FF
        FF=FSHAPE(UU-H6)
        RK7L(2)=RK7L(2)+FF
        IF(FF < EPS)		RETURN
     126 RK6P(2)=RK6P(2)+FF
        FF=FSHAPE(UU-H6-H7)
        RK6(2)=RK6(2)+FF
        FF=FSHAPE(UU-H5)
        RK6L(2)=RK6L(2)+FF
        IF(FF < EPS)		RETURN
     125 RK5P(2)=RK5P(2)+FF
        FF=FSHAPE(UU-H5-H6)
        RK5(2)=RK5(2)+FF
        FF=FSHAPE(UU-H4)
        RK5L(2)=RK5L(2)+FF
        IF(FF < EPS)		RETURN
     124 RK4P(2)=RK4P(2)+FF
        FF=FSHAPE(UU-H4-H5)
        RK4(2)=RK4(2)+FF
        FF=FSHAPE(UU-H3)
        RK4L(2)=RK4L(2)+FF
        IF(FF < EPS)		RETURN
     123 RK3P(2)=RK3P(2)+FF
        FF=FSHAPE(UU-H3-H4)
        RK3(2)=RK3(2)+FF
        FF=FSHAPE(UU-H2)
        RK3L(2)=RK3L(2)+FF
        IF(FF < EPS)		RETURN
     122 RK2P(2)=RK2P(2)+FF
        FF=FSHAPE(UU-H2-H3)
        RK2(2)=RK2(2)+FF
        FF=FSHAPE(UU-H1)
        RK2L(2)=RK2L(2)+FF
        IF(FF < EPS)		RETURN
     121 RK1P(2)=RK1P(2)+FF
        FF=FSHAPE(UU-H1-H2)
        RK1(2)=RK1(2)+FF
        FF=FSHAPE(UU-H0)
        RK1L(2)=RK1L(2)+FF
    !*					*** FINISH LEFT SIDE ***
     10			CONTINUE
                XXX=H0
                DO I = 2,NT0
                RK0P(I)=RK0P(I)+FF
                FF=FSHAPE(UU-XXX-H1)
                RK0(I)=RK0(I)+FF
                XXX=XXX+H0
                FF=FSHAPE(UU-XXX)
                RK0L(I)=RK0L(I)+FF
                IF(FF < EPS) GOTO 500
                END DO
     30		CONTINUE
     500		RETURN
    end subroutine leftLBL
    
    subroutine centerLBL(FREQ, UL, FSHAPE)
        ! calculation contributions to the subinterval on each grid
        ! FROM the CENTRAL PART of the extended subinterval: [startDeltaWV; endDeltaWV]
        
        IMPLICIT INTEGER*4 (I-N)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU, CONSER, UUU
        real(kind=DP) :: EPS, EPS4
        real :: FF, FA
        procedure(shape), pointer :: FSHAPE
	UU=UL-FREQ
		IF(UU >= deltaWV)		GOTO 40
	FF=FSHAPE(0.D0)
		IF(FF < EPS)		RETURN
!*							* left-right side *
		NPOINT=1
		CONSER=UU-H
		FA=FSHAPE(UU)
		EPS4=EPS*0.25
		IF(FA > EPS4) RK(1)=RK(1)+FA
		IF(UU < H)			GOTO 211
		I=0
		UUU=UU
		IF(UUU < H0+H0)	GOTO 201
				DO I =1,NT0
				UUU=UUU-H0
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 190
				RK0P(I)=RK0P(I)+FA
				RK0(I)=RK0(I)+FSHAPE(UUU+H1)
				RK0L(I)=RK0L(I)+FF
 190			FA=FF
				IF(UUU-H0 < H0) GOTO 201
				END DO
 201			I=I*2
				IF(UUU < H0)	GOTO 202
				IB=I+1
				DO I =IB,NT1
				UUU=UUU-H1
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 191
				RK1P(I)=RK1P(I)+FA
				RK1(I)=RK1(I)+FSHAPE(UUU+H2)
				RK1L(I)=RK1L(I)+FF
 191			FA=FF
				IF(UUU-H1 < H1)GOTO 202
				END DO
 202			I=I*2
				IF(UUU < H1)		GOTO 203
				IB=I+1
				DO I =IB,NT2
				UUU=UUU-H2
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 192
				RK2P(I)=RK2P(I)+FA
				RK2(I)=RK2(I)+FSHAPE(UUU+H3)
				RK2L(I)=RK2L(I)+FF
 192			FA=FF
				IF(UUU-H2 < H2) GOTO 203
				END DO
 203			I=I*2
				IF(UUU < H2)	GOTO 204
				IB=I+1
				DO I =IB,NT3
				UUU=UUU-H3
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 193
				RK3P(I)=RK3P(I)+FA
				RK3(I)=RK3(I)+FSHAPE(UUU+H4)
				RK3L(I)=RK3L(I)+FF
 193			FA=FF
				IF(UUU-H3 < H3)GOTO 204
				END DO
 204			I=I*2
				IF(UUU < H3) GOTO 205
				IB=I+1
				DO I =IB,NT4
				UUU=UUU-H4
				FF=FSHAPE(UUU)
				IF(FF < EPS) GOTO 194
				RK4P(I)=RK4P(I)+FA
				RK4(I)=RK4(I)+FSHAPE(UUU+H5)
				RK4L(I)=RK4L(I)+FF
 194			FA=FF
				IF(UUU-H4 < H4) GOTO 205
				END DO
 205			I=I*2
				IF(UUU < H4)	GOTO 206
				IB=I+1
				DO I =IB,NT5
				UUU=UUU-H5
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 195
				RK5P(I)=RK5P(I)+FA
				RK5(I)=RK5(I)+FSHAPE(UUU+H6)
				RK5L(I)=RK5L(I)+FF
 195			FA=FF
				IF(UUU-H5 < H5) GOTO 206
				END DO
 206			I=I*2
				IF(UUU < H5)	GOTO 207
				IB=I+1
				DO I =IB,NT6
				UUU=UUU-H6
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 196
				RK6P(I)=RK6P(I)+FA
				RK6(I)=RK6(I)+FSHAPE(UUU+H7)
				RK6L(I)=RK6L(I)+FF
 196			FA=FF
				IF(UUU-H6 < H6) GOTO 207
				END DO
 207			I=I*2
				IF(UUU < H6)	GOTO 208
				IB=I+1
				DO I =IB,NT7
				UUU=UUU-H7
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 197
				RK7P(I)=RK7P(I)+FA
				RK7(I)=RK7(I)+FSHAPE(UUU+H8)
				RK7L(I)=RK7L(I)+FF
 197			FA=FF
				IF(UUU-H7 < H7) GOTO 208
				END DO
 208			I=I*2
				IF(UUU < H7)	GOTO 209
				IB=I+1
				DO I =IB,NT8
				UUU=UUU-H8
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 198
				RK8P(I)=RK8P(I)+FA
				RK8(I)=RK8(I)+FSHAPE(UUU+H9)
				RK8L(I)=RK8L(I)+FF
 198			FA=FF
				IF(UUU-H8 < H8) GOTO 209
				END DO
 209			I=I*2
				IF(UUU < H8)	GOTO 210
				IB=I+1
				DO I =IB,NT9
				UUU=UUU-H9
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 199
				RK9P(I)=RK9P(I)+FA
				RK9(I)=RK9(I)+FSHAPE(UUU+H+H)
				RK9L(I)=RK9L(I)+FF
 199			FA=FF
				IF(UUU-H9 < H9) GOTO 210
				END DO
 210			I=I*4
				IB=I+2
				CONSER=UU-(IB-1)*H
				DO ICON=IB,NT
				RK(ICON)=RK(ICON)+FSHAPE(CONSER)
				CONSER=CONSER-H
				IF(CONSER < 0.) GOTO 490
				END DO
									GOTO 11
 490				NPOINT=ICON
!*							* right-LEFT side *
 211	NPOINT=NPOINT+1
		UUU=deltaWV-UU
		FA=FSHAPE(UUU)
				III=0
			IF(UUU < H0+H0)		GOTO 301
				DO I =NT0,1,-1
				III=III+1
				UUU=UUU-H0
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 290
				RK0L(I)=RK0L(I)+FA
				RK0(I)=RK0(I)+FSHAPE(UUU+H1)
				RK0P(I)=RK0P(I)+FF
 290			FA=FF
				IF(UUU-H0 < H0) GOTO 1301
				END DO
 301			IF(UUU < H0)	GOTO 302
 1301			III=III*2
				IB=NT1-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H1
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 291
				RK1L(I)=RK1L(I)+FA
				RK1(I)=RK1(I)+FSHAPE(UUU+H2)
				RK1P(I)=RK1P(I)+FF
 291			FA=FF
				IF(UUU-H1 < H1) GOTO 1302
				END DO
 302	IF(UUU < H1)			GOTO 303
 1302			III=III*2
				IB=NT2-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H2
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 292
				RK2L(I)=RK2L(I)+FA
				RK2(I)=RK2(I)+FSHAPE(UUU+H3)
				RK2P(I)=RK2P(I)+FF
 292			FA=FF
				IF(UUU-H2 < H2) GOTO 1303
				END DO
 303			IF(UUU < H2)	GOTO 304
 1303	III=III*2
				IB=NT3-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H3
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 293
				RK3L(I)=RK3L(I)+FA
				RK3(I)=RK3(I)+FSHAPE(UUU+H4)
				RK3P(I)=RK3P(I)+FF
 293			FA=FF
				IF(UUU-H3 < H3) GOTO 1304
				END DO
 304			IF(UUU < H3)	GOTO 305
 1304			III=III*2
				IB=NT4-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H4
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 294
				RK4L(I)=RK4L(I)+FA
				RK4(I)=RK4(I)+FSHAPE(UUU+H5)
				RK4P(I)=RK4P(I)+FF
 294			FA=FF
					IF(UUU-H4 < H4) GOTO 1305
				END DO
 305			IF(UUU < H4)	GOTO 306
 1305			III=III*2
				IB=NT5-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H5
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 295
				RK5L(I)=RK5L(I)+FA
				RK5(I)=RK5(I)+FSHAPE(UUU+H6)
				RK5P(I)=RK5P(I)+FF
 295			FA=FF
				IF(UUU-H5 < H5) GOTO 1306
				END DO
 306			IF(UUU < H5)	GOTO 307
 1306			III=III*2
				IB=NT6-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H6
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 296
				RK6L(I)=RK6L(I)+FA
				RK6(I)=RK6(I)+FSHAPE(UUU+H7)
				RK6P(I)=RK6P(I)+FF
 296			FA=FF
				IF(UUU-H6 < H6) GOTO 1307
				END DO
 307			IF(UUU < H6)	GOTO 308
 1307			III=III*2
				IB=NT7-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H7
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 297
				RK7L(I)=RK7L(I)+FA
				RK7(I)=RK7(I)+FSHAPE(UUU+H8)
				RK7P(I)=RK7P(I)+FF
 297			FA=FF
				IF(UUU-H7 < H7) GOTO 1308
				END DO
 308			IF(UUU < H7)	GOTO 309
 1308			III=III*2
				IB=NT8-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H8
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 298
				RK8L(I)=RK8L(I)+FA
				RK8(I)=RK8(I)+FSHAPE(UUU+H9)
				RK8P(I)=RK8P(I)+FF
 298			FA=FF
				IF(UUU-H8 < H8) GOTO 1309
				END DO
 309			IF(UUU < H8)	GOTO 310
 1309			III=III*2
				IB=NT9-III
				DO I =IB,1,-1
				III=III+1
				UUU=UUU-H9
				FF=FSHAPE(UUU)
				IF(FF < EPS)	GOTO 299
				RK9L(I)=RK9L(I)+FA
				RK9(I)=RK9(I)+FSHAPE(UUU+H+H)
				RK9P(I)=RK9P(I)+FF
 299			FA=FF
				IF(UUU-H9 < H9) GOTO 310
				END DO
 310			III=III*4
				I=NT-III
				DO II=NPOINT,I
				RK(II)=RK(II)+FSHAPE(CONSER)
 		CONSER=CONSER-H
                                END DO
!*						*** FINISH CENTRAL PART ***
 40					CONTINUE
 11							RETURN
    end subroutine centerLBL

    subroutine rightLBL(FREQ, UL, FSHAPE)
        ! calculation contributions to the subinterval on each grid
        ! FROM the RIGHT PART of the extended subinterval: [endDeltaWV; endDeltaWV+cutOff]

        IMPLICIT INTEGER*4 (I-N)
        real(kind=DP) :: FREQ
        real(kind=DP) :: UL, UU
        real(kind=DP) :: EPS
        real(kind=DP) :: XXX
        real :: FF
        procedure(shape), pointer :: FSHAPE

	UU=UL-FREQ-deltaWV ! CHANGE DIAP
	IF(UU >= cutOff) RETURN
	FF=FSHAPE(UU)
		IF(FF < EPS)			RETURN
		IF(UU < H0)			GOTO 60
	RK0L(NT0)=RK0L(NT0)+FF
	FF=FSHAPE(UU+H1)
	RK0(NT0)=RK0(NT0)+FF
	FF=FSHAPE(UU+H0)
	RK0P(NT0)=RK0P(NT0)+FF
								GOTO 12
 60	IF(UU < H1)			GOTO 61
	RK1L(NT1)=RK1L(NT1)+FF
	FF=FSHAPE(UU+H2)
	RK1(NT1)=RK1(NT1)+FF
	FF=FSHAPE(UU+H1)
	RK1P(NT1)=RK1P(NT1)+FF
		IF(FF < EPS)			RETURN
								GOTO 131
!* 61	IF(-UU < H2.and.a < h2)goto 22
 61	IF(UU < H2)			GOTO 62
	RK2L(NT2)=RK2L(NT2)+FF
	FF=FSHAPE(UU+H3)
	RK2(NT2)=RK2(NT2)+FF
	FF=FSHAPE(UU+H2)
	RK2P(NT2)=RK2P(NT2)+FF
		IF(FF < EPS)			RETURN
								GOTO 132
 62	IF(UU < H3)			GOTO 63
	RK3L(NT3)=RK3L(NT3)+FF
	FF=FSHAPE(UU+H4)
	RK3(NT3)=RK3(NT3)+FF
	FF=FSHAPE(UU+H3)
	RK3P(NT3)=RK3P(NT3)+FF
		IF(FF < EPS)			RETURN
								GOTO 133
 63	IF(UU < H4)			GOTO 64
	RK4L(NT4)=RK4L(NT4)+FF
	FF=FSHAPE(UU+H5)
	RK4(NT4)=RK4(NT4)+FF
	FF=FSHAPE(UU+H4)
	RK4P(NT4)=RK4P(NT4)+FF
		IF(FF < EPS)			RETURN
								GOTO 134
 64	IF(UU < H5)			GOTO 65
	RK5L(NT5)=RK5L(NT5)+FF
	FF=FSHAPE(UU+H6)
	RK5(NT5)=RK5(NT5)+FF
	FF=FSHAPE(UU+H5)
	RK5P(NT5)=RK5P(NT5)+FF
		IF(FF < EPS)			RETURN
								GOTO 135
 65	IF(UU < H6)			GOTO 66
	RK6L(NT6)=RK6L(NT6)+FF
	FF=FSHAPE(UU+H7)
	RK6(NT6)=RK6(NT6)+FF
	FF=FSHAPE(UU+H6)
	RK6P(NT6)=RK6P(NT6)+FF
		IF(FF < EPS)			RETURN
								GOTO 136
 66	IF(UU < H7)			GOTO 67
	RK7L(NT7)=RK7L(NT7)+FF
	FF=FSHAPE(UU+H8)
	RK7(NT7)=RK7(NT7)+FF
	FF=FSHAPE(UU+H7)
	RK7P(NT7)=RK7P(NT7)+FF
		IF(FF < EPS)			RETURN
								GOTO 137
 67	IF(UU < H8)			GOTO 68
	RK8L(NT8)=RK8L(NT8)+FF
	FF=FSHAPE(UU+H9)
	RK8(NT8)=RK8(NT8)+FF
	FF=FSHAPE(UU+H8)
	RK8P(NT8)=RK8P(NT8)+FF
		IF(FF < EPS)			RETURN
								GOTO 138
 68	IF(UU < H9)			GOTO 69
	RK9L(NT9)=RK9L(NT9)+FF
	FF=FSHAPE(UU+H+H)
	RK9(NT9)=RK9(NT9)+FF
	FF=FSHAPE(UU+H9)
	RK9P(NT9)=RK9P(NT9)+FF
		IF(FF < EPS)			RETURN
								GOTO 139
 69	RK(NT)=RK(NT)+FSHAPE(UU)
	RK(NT-1)=RK(NT-1)+FSHAPE(UU+H)
	RK(NT-2)=RK(NT-2)+FSHAPE(UU+H+H)
	RK(NT-3)=RK(NT-3)+FSHAPE(UU+H9-H)
!* CASCADE
 139	N2=NT9-1
	FF=FSHAPE(UU+H9)
	RK9L(N2)=RK9L(N2)+FF
	FF=FSHAPE(UU+H9+H+H)
	RK9(N2)=RK9(N2)+FF
	FF=FSHAPE(UU+H8)
	RK9P(N2)=RK9P(N2)+FF
 138 N=NT8-1
	RK8L(N)=RK8L(N)+FF
	FF=FSHAPE(UU+H8+H9)
	RK8(N)=RK8(N)+FF
	FF=FSHAPE(UU+H7)
	RK8P(N)=RK8P(N)+FF
		IF(FF < EPS)			RETURN
 137 N=NT7-1
	RK7L(N)=RK7L(N)+FF
	FF=FSHAPE(UU+H7+H8)
	RK7(N)=RK7(N)+FF
	FF=FSHAPE(UU+H6)
	RK7P(N)=RK7P(N)+FF
		IF(FF < EPS)			RETURN
 136 N=NT6-1
	RK6L(N)=RK6L(N)+FF
	FF=FSHAPE(UU+H6+H7)
	RK6(N)=RK6(N)+FF
	FF=FSHAPE(UU+H5)
	RK6P(N)=RK6P(N)+FF
		IF(FF < EPS)			RETURN
 135 N=NT5-1
	RK5L(N)=RK5L(N)+FF
	FF=FSHAPE(UU+H5+H6)
	RK5(N)=RK5(N)+FF
	FF=FSHAPE(UU+H4)
	RK5P(N)=RK5P(N)+FF
		IF(FF < EPS)			RETURN
 134 N=NT4-1
	RK4L(N)=RK4L(N)+FF
	FF=FSHAPE(UU+H4+H5)
	RK4(N)=RK4(N)+FF
	FF=FSHAPE(UU+H3)
	RK4P(N)=RK4P(N)+FF
		IF(FF < EPS)			RETURN
 133 N=NT3-1
	RK3L(N)=RK3L(N)+FF
	FF=FSHAPE(UU+H3+H4)
	RK3(N)=RK3(N)+FF
	FF=FSHAPE(UU+H2)
	RK3P(N)=RK3P(N)+FF
		IF(FF < EPS)			RETURN
 132 N=NT2-1
	RK2L(N)=RK2L(N)+FF
	FF=FSHAPE(UU+H2+H3)
	RK2(N)=RK2(N)+FF
	FF=FSHAPE(UU+H1)
	RK2P(N)=RK2P(N)+FF
		IF(FF < EPS)			RETURN
 131 N=NT1-1
	RK1L(N)=RK1L(N)+FF
	FF=FSHAPE(UU+H1+H2)
	RK1(N)=RK1(N)+FF
	FF=FSHAPE(UU+H0)
	RK1P(N)=RK1P(N)+FF
!*
!*						*** FINISH RIGHT PART ***
 12		CONTINUE
		XXX=H0
		DO I = NT0-1,1,-1
		RK0L(I)=RK0L(I)+FF
		FF=FSHAPE(UU+XXX+H1)
		RK0(I)=RK0(I)+FF
		XXX=XXX+H0
		FF=FSHAPE(UU+XXX)
		RK0P(I)=RK0P(I)+FF
		IF(FF < EPS)GOTO 1000
		END DO
			RK(1)=RK(1)+FF
 1000							RETURN
    end subroutine rightLBL

end module LineGridCalc
