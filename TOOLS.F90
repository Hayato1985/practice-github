! ======================================================================
    SUBROUTINE INVERT(MN,N,A)
! ======================================================================
! É}ÉgÉäÉNÉXÇÃãtçsóÒÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    PARAMETER(MCHCK=0)
    PARAMETER(ERR=1.0D-20)
    DIMENSION A(MN,MN)
    DOUBLEPRECISION, ALLOCATABLE :: B(:,:)
    DOUBLEPRECISION, ALLOCATABLE :: C(:,:)
    DOUBLEPRECISION, ALLOCATABLE :: D(:,:)
!
    ALLOCATE(B(MN,MN),STAT=IERROR)
    B=A
    IF(MCHCK.GE.2)THEN
      ALLOCATE(C(MN,MN),STAT=IERROR)
      ALLOCATE(D(MN,MN),STAT=IERROR)
      C=A
    ENDIF
!
    CALL ZEROR2(A,MN,MN,N,N)
    DO IPVT=1,N
      A(IPVT,IPVT)=1.0D00
    ENDDO
!
    DO IPVT=1,N
      DIAG=B(IPVT,IPVT)
      IF(ABS(DIAG).LT.ERR)THEN
        WRITE(8,*)' DIAGONAL ELEMENT',IPVT,' IS NEARLY ZERO',DIAG
      ELSE
        IF(ABS(DIAG).LE.1.0D-4)THEN
          WRITE(8,*)' WARNING! DIAGONAL ELEMENT',IPVT,' IS SMALL',DIAG
        ENDIF
        DO ILINE=1,N
          B(IPVT,ILINE)=B(IPVT,ILINE)/DIAG                       ! çsIPVTÇÃëŒäpóvëfÇÇPÇ…Ç∑ÇÈ
          A(IPVT,ILINE)=A(IPVT,ILINE)/DIAG                       ! çsIPVTÇÃëŒäpóvëfÇÇPÇ…Ç∑ÇÈ
        ENDDO
        DO IROW=1,N
          IF(IROW.NE.IPVT)THEN
            COEF=B(IROW,IPVT)
            DO ILINE=1,N
              B(IROW,ILINE)=B(IROW,ILINE)-COEF*B(IPVT,ILINE)      ! çsIROWÇ©ÇÁçsIPVTÇÃCOEFî{Çà¯Ç≠
              A(IROW,ILINE)=A(IROW,ILINE)-COEF*A(IPVT,ILINE)      ! çsIROWÇ©ÇÁçsIPVTÇÃCOEFî{Çà¯Ç≠
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
!
!***** ãtçsóÒÇÃê∏ìxÇÃÉ`ÉFÉbÉN
 IF(MCHCK.GE.2)WRITE(8,*) ' '
 IF(MCHCK.GE.2)WRITE(8,*) 'ãtçsóÒÇÃê∏ìxÇÃÉ`ÉFÉbÉN'
 IF(MCHCK.GE.2)CALL MLTPLY(N,N,N,A,C,D)
 IF(MCHCK.GE.2)CALL WMTR2(2,C ,MN,MN,N,N,'A    ','     ',0)
 IF(MCHCK.GE.2)CALL WMTR2(2,A ,MN,MN,N,N,'Ainv ','     ',0)
 IF(MCHCK.GE.2)CALL WMTR2(2,D ,MN,MN,N,N,'A*Ain','     ',0)
 IF(MCHCK.GE.2)CALL WMTR2(2,B ,MN,MN,N,N,'B    ','     ',0)
!*****
!
    RETURN
    END            
!
!
!======================================================================
    SUBROUTINE MLTPLY(N1,N2,N3,A,B,C)
!======================================================================
! É}ÉgÉäÉNÉXÇÃêœÇbÅÅÇ`*ÇaÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    DIMENSION A(N1,N2),B(N2,N3),C(N1,N3)
    DOUBLEPRECISION, ALLOCATABLE :: AB(:,:)
    ALLOCATE(AB(N1,N3),STAT=IERROR)
!
    DO I=1,N1
      DO J=1,N3
        CC=0.0
        DO K=1,N2
          CC=CC+A(I,K)*B(K,J)
        END DO
        AB(I,J)=CC
      END DO
    END DO
    C=AB
    RETURN
    END
!
!
!======================================================================
    SUBROUTINE MLTPLY_ABC(N1,N2,N3,N4,A,B,C,D)
!======================================================================
! É}ÉgÉäÉNÉXÇÃÇRèdêœÇcÅÅÇ`*Ça*ÇbÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    DIMENSION A(N1,N2),B(N2,N3),C(N3,N4),D(N1,N4)
    DOUBLEPRECISION, ALLOCATABLE :: AB(:,:)
    ALLOCATE(AB(N1,N3),STAT=IERROR)
!
    CALL MLTPLY(N1,N2,N3,A,B,AB)
    CALL MLTPLY(N1,N3,N4,AB,C,D)
    RETURN
    END
!
!
!======================================================================
    SUBROUTINE MLTPLY_ATBC(N1,N2,N3,N4,A,B,C,D)
!======================================================================
! É}ÉgÉäÉNÉXÇÃÇRèdêœÇcÅÅÇ`(ì]íu)*Ça*ÇbÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    DIMENSION A(N2,N1),B(N2,N3),C(N3,N4),D(N1,N4)
    DOUBLEPRECISION, ALLOCATABLE :: AT(:,:)
    ALLOCATE(AT(N1,N2),STAT=IERROR)
!
    CALL TRANS(N1,N2,A,AT)
    CALL MLTPLY_ABC(N1,N2,N3,N4,AT,B,C,D)
    RETURN
    END
!
!
!
!======================================================================
    SUBROUTINE MLTPLY_ATB(N1,N2,N3,A,B,C)
!======================================================================
! É}ÉgÉäÉNÉXÇÃêœÇbÅÅÇ`Åiì]íuÅj*ÇaÇÃåvéZ
    IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    DIMENSION A(N2,N1),B(N2,N3),C(N1,N3)
    DOUBLEPRECISION, ALLOCATABLE :: AT(:,:)
    ALLOCATE (AT(N1,N2),STAT=IERROR)
!
    CALL TRANS(N2,N1,A,AT)
    CALL MLTPLY(N1,N2,N3,AT,B,C)
    RETURN
    END
!
!
!======================================================================
      SUBROUTINE MLTPLY_ABT(N1,N2,N3,A,B,C)
!======================================================================
! É}ÉgÉäÉNÉXÇÃêœÇbÅÅÇ`*ÇaÅiì]íuÅjÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    DIMENSION A(N1,N2),B(N3,N2),C(N1,N3)
    DOUBLEPRECISION, ALLOCATABLE :: BT(:,:)
    ALLOCATE (BT(N2,N3),STAT=IERROR)
!
    CALL TRANS(N3,N2,B,BT)
    CALL MLTPLY(N1,N2,N3,A,BT,C)
    RETURN
    END
!
!
!======================================================================
    SUBROUTINE MLTPLY_AS(N1,N2,A,SCALAR)
!======================================================================
! É}ÉgÉäÉNÉXÇ∆ÉXÉJÉâÅ[ÇÃêœÇaÅÅÇ`*sÇÃåvéZ
    IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    DIMENSION A(N1,N2)
!
    DO IROW=1,N1
      DO ILINE=1,N2
        A(IROW,ILINE)=A(IROW,ILINE)*SCALAR
      ENDDO
    ENDDO
    RETURN
      END
!
!
!======================================================================
    SUBROUTINE MLTPLY_AADDBS(N1,N2,A,B,SCALAR)
!======================================================================
! É}ÉgÉäÉNÉXÇ`Ç…É}ÉgÉäÉNÉXÇaÇ∆ÉXÉJÉâÅ[ÇÃêœÇâ¡Ç¶ÇÈÇ`ÅÅÇ`Å{Ça*sÇÃåvéZ
!
    DOUBLEPRECISION A(N1,N2)
    DOUBLEPRECISION B(N1,N2)
    DOUBLEPRECISION SCALAR  
!
    DO IROW=1,N1
      DO ILINE=1,N2
        A(IROW,ILINE)=A(IROW,ILINE)+B(IROW,ILINE)*SCALAR
      ENDDO
    ENDDO
    RETURN
    END
!
!
!======================================================================
    SUBROUTINE SUMMT(N1,N2,A,B,C)
!======================================================================
! É}ÉgÉäÉNÉXÇ∆É}ÉgÉäÉNÉXÇÃòaÇbÅÅÇ`Å{ÇaÇÃåvéZ
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N1,N2),B(N1,N2),C(N1,N2)
!
    DO IROW=1,N1
      DO ILINE=1,N2
        C(IROW,ILINE)=A(IROW,ILINE)+B(IROW,ILINE)
      ENDDO
    ENDDO
    RETURN
      END
!
!
!
!======================================================================
      SUBROUTINE TRANS(NL,NR,A,AT)
!======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION A(NL,NR)
	DIMENSION AT(NR,NL)
!
      DO IL=1,NL
	  DO IR=1,NR
	    AT(IR,IL)=A(IL,IR)
        END DO
      END DO
      RETURN 
	END
!
! ======================================================================
      SUBROUTINE ZEROR1(A,MR,NR)
! ======================================================================
      DOUBLE PRECISION A(MR)
      DO 100 I=1,NR
        A(I)=0.0
  100 CONTINUE
      RETURN
      END
!
      SUBROUTINE ZEROR2(A,MR,MC,NR,NC)
      DOUBLE PRECISION A(MR,MC)
      DO 100 I=1,NR
      DO 100 J=1,NC
        A(I,J)=0.0
  100 CONTINUE
      RETURN
      END
!
      SUBROUTINE ZEROR3(A,M1,M2,M3,N1,N2,N3)
      DOUBLE PRECISION A(M1,M2,M3)
      DO 100 I=1,N1
      DO 100 J=1,N2
      DO 100 K=1,N3
        A(I,J,K)=0.0
  100 CONTINUE
      RETURN
      END
!
      SUBROUTINE ZEROI1(NN,MR,NR)
      DIMENSION NN(MR)
      DO 100 I=1,NR
        NN(I)=0
  100 CONTINUE
      RETURN
      END
!
      SUBROUTINE ZEROI2(NN,MR,MC,NR,NC)
      DIMENSION NN(MR,MC)
      DO 100 I=1,NR
      DO 100 J=1,NC
        NN(I,J)=0
  100 CONTINUE
      RETURN
      END
!
! ======================================================================
      SUBROUTINE WMTR2(INDEX,X,ML,MR,NL,NR,NM1,NM2,NC)
! ======================================================================
!       X(ML,MR) 1-NL:LINE 1-NR:ROW 
!       NM1:NAME OF ARRAY NM2,NC:MEMO
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(ML,MR)
      CHARACTER*5 NM1,NM2
      WRITE(8,1)NM1,NL,NR,NM2,NC
      NLT=INT(NL/20)
      NRT=INT(NR/20)
      DO 200 II=1,NLT+1
        NL1=20*(II-1)+1
        NL2=20*II
        IF(NL2.GT.NL)NL2=NL
        IF(NL1.GT.NL)GOTO 200
        DO 100 JJ=1,NRT+1
          NR1=20*(JJ-1)+1
          NR2=20*JJ
          IF(NR2.GT.NR)NR2=NR
          IF(NR1.GT.NR)GOTO 200
          WRITE(8,2) (IR,IR=NR1,NR2)
         DO 100 I=NL1,NL2
          IF(INDEX.EQ.1)WRITE(8,3) I,(X(I,J),J=NR1,NR2)
          IF(INDEX.EQ.2)WRITE(8,4) I,(X(I,J),J=NR1,NR2)
  100  CONTINUE
  200  CONTINUE
    1 FORMAT(A5,'[',I3,','I3,' ]',A5,'=',I3)
    2 FORMAT(5X, 20(I7,4X))
    3 FORMAT(I5,20(1X,F10.4))
    4 FORMAT(I5,20(1X,E10.4))
      RETURN
      END
!
! ======================================================================
      SUBROUTINE WMTR1(INDEX,X,ML,MR,NL,NR,NM1,NM2,NC)
! ======================================================================
!       X(ML) 1-NL:LINE
!       NM1:NAME OF ARRAY NM2,NC:MEMO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(ML)
      CHARACTER*5 NM1,NM2
      IF(INDEX.LE.2)NDEVD=20
 	IF(INDEX.EQ.3)NDEVD=6
      WRITE(8,2)NM1,NL,NM2,NC,(I,I=1,NDEVD)
      NLT=NL/NDEVD
      DO 200 II=1,NLT+1
        NL1=NDEVD*(II-1)+1
        NL2=NDEVD*II
        IF(NL2.GT.NL)NL2=NL
        IF(NL1.GT.NL)GOTO 200
        IF(INDEX.EQ.1)WRITE(8,3)NL1-1,(X(I),I=NL1,NL2)
        IF(INDEX.EQ.2)WRITE(8,4)NL1-1,(X(I),I=NL1,NL2)
	  IF(INDEX.EQ.3)WRITE(8,5)II,(X(I),I=NL1,NL2)
  200  CONTINUE
    2 FORMAT(A5,'[',I3,']',5X,A5,'=',I3,/,5X,20(I7,4X))
    3 FORMAT('+',I3,1X,20(1X,F10.4))
    4 FORMAT('+',I3,1X,20(1X,E10.4))
    5 FORMAT(I4,12(E11.4))
      RETURN
      END
!
!
! ======================================================================
      SUBROUTINE WMTI2(INDEX,NN,ML,MR,NL,NR,NM1,NM2,NC)
! ======================================================================
!       NN(ML,MR) 1-NL:LINE 1-NR:ROW 
!       NM1:NAME OF ARRAY NM2,NC:MEMO
!
      DIMENSION NN(ML,MR)
      CHARACTER*5 NM1,NM2
      IF(INDEX.EQ.1)THEN
       NLT=NL/20
       NRT=NR/20
       DO 200 II=1,NLT+1
       NL1=20*(II-1)+1
       NL2=20*II
       IF(NL2.GT.NL)NL2=NL
       IF(NL1.GT.NL)GOTO 200
       DO 100 JJ=1,NRT+1
       NR1=20*(JJ-1)+1
       NR2=20*JJ
       IF(NR2.GT.NR)NR2=NR
       IF(NR1.GT.NR)GOTO 100
       WRITE(8,4) NM1,NL,NR,NM2,NC,NL1,NR1,NL2,NR2,(I,(NN(I,J),J=NR1,NR2),I=NL1,NL2)
  100  CONTINUE
  200  CONTINUE
      ELSEIF(INDEX.EQ.2)THEN
       IF(NR.EQ.3)WRITE(8,3)NM1,NL,NR,NM2,NC,(K,K=1,NR),(I,(NN(I,J),J=1,NR),I=1,NL)
      ENDIF
    4 FORMAT(A5,'[',I3,','I3,' ]',A5,'=',I3,5X,'[',2I3,']-[',2I3,']',/,(I5,20(1X,I5)))
    3 FORMAT(A5,'[',I3,','I3,' ]',A5,'=',I3,/,5X, 3I5,/,(4I5))
      RETURN
      END
!
! ======================================================================
      SUBROUTINE WMTI1(INDEX,NN,ML,MR,NL,NR,NM1,NM2,NC)
! ======================================================================
!       X(ML) 1-NL:LINE
!       NM1:NAME OF ARRAY NM2,NC:MEMO
!
      DIMENSION NN(ML)
      CHARACTER*5 NM1,NM2
       NROW=6
       NLT=NL/NROW
       DO 200 II=1,NLT+1
       NL1=NROW*(II-1)+1
       NL2=NROW*II
       IF(NL2.GT.NL)NL2=NL
       IF(NL1.GT.NL)GOTO 200
       WRITE(8,3)NM1,NL,NR,NM2,NC,NL1,NL2,(NN(I),I=NL1,NL2)
  100  CONTINUE
  200  CONTINUE
    3 FORMAT(A5,'[',I4,','I2,' ]',A5,'=',I3,7X,'[',I4,'-',I4,']',24(I5))
      RETURN
      END
!
! ======================================================================
      SUBROUTINE WMTR3(INDEX,X,M1,M2,M3,N1,N2,N3,NM1,NM2,NC)
! ======================================================================
!       X(M1,M2,M3) 1-N1, 1-N2, 1-N3 
!       NM1:NAME OF ARRAY NM2,NC:MEMO
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(M1,M2,M3)
      CHARACTER*5 NM1,NM2
      NLT=INT(N1/20)
      NRT=INT(N2/20)
      DO 200 KK=1,N3
      WRITE(8,1)NM1,N1,N2,KK,NM1,N1,N2,N3,NM2,NC
      DO 200 II=1,NLT+1
        NL1=20*(II-1)+1
        NL2=20*II
        IF(NL2.GT.N1)NL2=N1
        IF(NL1.GT.N1)GOTO 200
        DO 100 JJ=1,NRT+1
          NR1=20*(JJ-1)+1
          NR2=20*JJ
          IF(NR2.GT.N2)NR2=N2
          IF(NR1.GT.N2)GOTO 200
          WRITE(8,2) (IR,IR=NR1,NR2)
         DO 100 I=NL1,NL2
          IF(INDEX.EQ.1)WRITE(8,3) I,(X(I,J,KK),J=NR1,NR2)
          IF(INDEX.EQ.2)WRITE(8,4) I,(X(I,J,KK),J=NR1,NR2)
  100  CONTINUE
  200  CONTINUE
    1 FORMAT(A5,'[',3I3,' ] OF ',A5,'[',3I3,' ]',A5,'=',I3)
    2 FORMAT(5X, 20(I7,4X))
    3 FORMAT(I5,20(1X,F10.4))
    4 FORMAT(I5,20(1X,E10.4))
      RETURN
      END
!
! ======================================================================
      SUBROUTINE WRTLD(NM1,NM2,XLOAD,LNODS,NELEM,NPOIN,NNODE)
! ======================================================================
!     NM1:NAME OF ARRAY NM2,NC:MEMO
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*5 NM1,NM2
      DIMENSION XLOAD(NELEM,*),LNODS(NELEM,*),TOTAL(300,5),SUMLD(5)
!
! *** CALCULATE AND OUTPUT THE EACH NODAL LOAD DUE TO THE XLOAD
!
      CALL ZEROR2(TOTAL,300,5,NPOIN,5)
      DO 100 IELEM=1,NELEM
      DO 100 INODE=1,NNODE
        IPOIN=IABS(LNODS(IELEM,INODE))
        DO 100 IDOFN=1,5
          IEVAB=(INODE-1)*5+IDOFN
          TOTAL(IPOIN,IDOFN)=TOTAL(IPOIN,IDOFN)+XLOAD(IELEM,IEVAB)
  100 CONTINUE
      CALL ZEROR1(SUMLD,5,5)
      DO 200 IPOIN=1,NPOIN
      DO 200 IDOFN=1,5
        SUMLD(IDOFN)=SUMLD(IDOFN)+TOTAL(IPOIN,IDOFN)
  200 CONTINUE
      KPOIN=NPOIN
      CALL WMTR2(2,TOTAL,300,5,KPOIN,5,NM1,NM2,0)
      WRITE(8,1)(SUMLD(IDOFN),IDOFN=1,5)
    1 FORMAT(2X,60('-')/(5X,5E11.4))
      RETURN
      END
!
!
! ======================================================================
      SUBROUTINE WRTLD1(NM1,NM2,XLOAD,MTOTV,NTOTV,NPOIN)
! ======================================================================
!       NM1:NAME OF ARRAY NM2,NC:MEMO
      PARAMETER(NDOFN=5)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*5 NM1,NM2
      DIMENSION XLOAD(MTOTV),SUMLD(NDOFN)
!
      CALL ZEROR1(SUMLD,NDOFN,NDOFN)
      DO 200 IDOFN=1,NDOFN
      DO 200 IPOIN=1,NPOIN
	  ITOTV=(IPOIN-1)*NDOFN+IDOFN
        SUMLD(IDOFN)=SUMLD(IDOFN)+XLOAD(ITOTV)
  200 CONTINUE
      KPOIN=NPOIN
      CALL WMTR1(3,XLOAD,MTOTV,0,NTOTV,0,NM1,NM2,0)
      WRITE(8,1)(SUMLD(IDOFN),IDOFN=1,5)
    1 FORMAT(2X,60('-')/(4X,5E11.4))
      RETURN
      END

