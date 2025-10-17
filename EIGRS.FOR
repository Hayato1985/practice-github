C ======================================================================
      SUBROUTINE EIGEN(EIGN1,MSIZE,MTOTV,NFREE,NORDR,NPOIN,NWDTH,NWSUM,
     &                 RLOAD,TSTIF)
C ======================================================================
C *** 全体剛性マトリクスTSTIFの固有値を計算する
C     SUB EIGRS：実対称行列の固有値を求めるプログラム
C       A( , )：実対称行列
C       NE    ：求める固有値の個数、NE<0のとき小さいものから順番になる
C       NV    ：求める固有ベクトルの個数、NV>0のとき正規直交化を行う
C       E( )  ：求められた固有値
C       V( )  ：求められた固有ベクトル
      PARAMETER(MMTOTV=1000,MCHK=0)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NWDTH(MTOTV),NWSUM(0:MTOTV),TSTIF(MSIZE),RLOAD(MTOTV)
      DIMENSION A(MMTOTV,MMTOTV),E(MMTOTV),LW(MMTOTV),V(MMTOTV,MMTOTV),
     &          W(MMTOTV,7),NORDR(MTOTV)
      DIMENSION VECTOR(5)
C
      IF(MMTOTV.LT.NFREE)THEN
        WRITE(6,9003)MMTOTV,NFREE
        WRITE(8,9003)MMTOTV,NFREE
 9003   FORMAT(' (SUB EIGEN) MMTOTV:',I5,' IS LESS THAN NFREE:',I5)
        STOP
      ENDIF
C
      DO 100 ITOTV=1,NFREE
        IWDTH=NWDTH(ITOTV)
        DO 100 JTOTV=1,ITOTV
          IF(JTOTV.GT.(ITOTV-IWDTH))THEN
            IPOSI=NWSUM(ITOTV)-ITOTV+JTOTV
            A(ITOTV,JTOTV)=TSTIF(IPOSI)
          ELSE
            A(ITOTV,JTOTV)=0.0
          ENDIF
  100 CONTINUE
C
C 安定性判定のための最小固有値計算<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	DO 150 ITOTV=1,NFREE
	  W(ITOTV,1)=RLOAD(ITOTV)
  150 CONTINUE
C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      N=NFREE
      N1=MMTOTV
      NE=-1
      NV=1
      EPS=1.0E-16
      CALL EIGRS(A, N, N1, NE, NV, EPS, W, LW, E, V, IER)
      EIGN1=E(1)
      SPRD1=0.0
      DO 200 ITOTV=1,NFREE
        SPRD1=SPRD1+RLOAD(ITOTV)*V(ITOTV,1)
  200 CONTINUE
C
      WRITE(8,9001)IER,EIGN1,SPRD1
 9001 FORMAT(' ***** EIGRS ( IER=',I3,' EIGN1=',E11.4,
     # ' SPRD1=',E11.4,' ) *****')
C      CALL WMTR1(2,E,MMTOTV,0,-NE,0,'EIGEN','EIGRS',0)
C      CALL WMTR2(2,V,MMTOTV,MMTOTV,NTOTV,NV,'EIGVT','EIGRS',0)
C     WRITE(8,9002)(V(NORDR(ITOTV),1),ITOTV=3,NFREE,5)
C9002 FORMAT(11E11.4)
C
C *** 最小固有値に対応する固有ベクトルの出力
C
      IF(MCHK.EQ.1)THEN
        NDOFN=6
        DO 350 IPOIN=1,NPOIN
          DO 300 IDOFN=1,NDOFN
            ITOTV=(IPOIN-1)*NDOFN+IDOFN
            IFREE=NORDR(ITOTV)
            IF(IFREE.LE.NFREE)THEN
              VECTOR(IDOFN)=V(IFREE,1)
            ELSE
              VECTOR(IDOFN)=0.0
            ENDIF
  300     CONTINUE
          WRITE(8,9004)IPOIN,(VECTOR(IDOFN),IDOFN=1,NDOFN)
  350   CONTINUE
      ENDIF
 9004 FORMAT(I5,6(E15.5))
C
      RETURN
      END
C ======================================================================
      SUBROUTINE EIGRS(A, N, N1, NE, NV, EPS, W, LW, E, V, IER)
C ======================================================================
************************************************************************
*  EIGENVALUES AND EIGENVECTORS FOR A SYMMETRIC MATRIX                 *
*     BY HOUSEHOLDER-BISECTION-INVERSE ITERATION METHOD.               *
*  PARAMETERS                                                          *  
*    (1) A: 2-DIM. ARRAY CONTAINING THE SYMMETRIC MATRIX               * 
*    (2) N: ORDER OF THE MATRIX (A)                                    *
*    (3) N1: ROW SIZE OF THE 2-DIM. ARRAY (A)                          *
*    (4) NE: NUMBER OF NEEDED EIGENVALUES                              *
*    (5) NV: NUMBER OF NEEDED EIGENVECTORS                             *
*    (6) EPS: TOLERANCE FOR CONVERGENCE                                *
*    (7) W: 2-DIM. WORKING ARRAY                                       *
*    (8) LW: 1-DIM. WORKING ARRAY                                      *
*    (9) E: 1-DIM. ARRAY CONTAINING COMPUTED EIGENVALUES               *
*   (10) V: 2-DIM. ARRAY CONTAINING COMPUTED EIGENVECTORS              *
*   (11) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI   MARCH 25 1990    VERSION 1.1                 *
************************************************************************  
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 SW
      DIMENSION A(N1,N), W(N1,7), LW(N), E(N1), V(N1,1)
C
      NEA = IABS(NE)
      NVA = IABS(NV)
      IF (NEA .EQ. 0 .OR. N1 .LT. N .OR. N .LT. NEA .OR. NEA .LT. 
     * NVA .OR. N .LT. 2) THEN
       WRITE(8,*) '(SUBR. EIGRS) INVALID ARGUMENT.',NV,NE,N,N1
       IER = 2
       RETURN
      ENDIF
      IF (EPS .LT. 0.0) EPS = 1.0D-16
      IF (N .EQ. 2) THEN
       W(1,1) = A(2,1)
       T = 0.5D0 * (A(1,1) + A(2,2))
       R = A(1,1) * A(2,2) - A(2,1) * A(2,1)
       D = T * T - R
       Q = DABS(T) + DSQRT(D)
       IF (T .LT. 0.0) Q = - Q
       T = T * DFLOAT(NE)
       IF (T .GE.0.0) THEN
        E(1) = Q
        IF (NEA .EQ. 2) E(2) = R / Q
       ELSE
        E(1) = R / Q
        IF (NEA .EQ. 2) E(2) = Q
       ENDIF
C  REDUCE TO TRIDIAGONAL FORM BY HOUSEHOLDER'S METHOD.
      ELSE
       DO 130 K=1,N-2
        S = 0.0D0
        DO 60 I=K+1,N
   60    S = S + A(I,K) * A(I,K)
        W(K,1) = 0.0D0
       IF (S .NE. 0) THEN
        SR = DSQRT(S)
        A1 = A(K+1,K)
        IF (A1 .LT. 0.0) SR = - SR
        W(K,1) = - SR
        R = 1.0D0 / (S + A1 * SR)
        A(K+1,K) = A1 + SR
        DO 90 I=K+1,N
         S = 0.0D0
         DO 70 J=K+1,I
   70     S = S + A(I,J) * A(J,K)
         IF (I .NE. N) THEN
          DO 80 J=I+1,N
   80      S = S + A(J,I) * A(J,K)
         ENDIF
   90    W(I,1) = S * R
         S = 0.0D0
         DO 100 I=K+1,N
  100     S = S + A(I,K) * W(I,1)
         T = 0.5D0 * S * R
         DO 110 I=K+1,N
  110     W(I,1) = W(I,1) - T * A(I,K)
         DO 120 J=K+1,N
          WJ1 = W(J,1)
          AJK = A(J,K)
          DO 120 I=J,N
  120      A(I,J) = A(I,J) - A(I,K)*WJ1 - W(I,1)*AJK
        ENDIF
  130 CONTINUE
        W(N-1,1) = A(N,N-1)
C  COMPUTE EIGENVALUES BY BISECTION METHOD.
         DO 135 I=1,N
  135     W(I,6) = A(I,I)
      R=DMAX1((DABS(W(1,6))+DABS(W(1,1))),(DABS(W(N-1,1))
     *   +DABS(W(N,6))))
      DO 140 I=2,N-1
       T=DABS(W(I-1,1))+DABS(W(I,6))+DABS(W(I,1))
       IF (T .GT. R) R = T
  140 CONTINUE
      EPS1 = R * 1.0D-16
      EPS2 = R * EPS
      DO 150 I=1,N-1
  150  W(I,2) = W(I,1) * W(I,1)
      IF (NE .LT. 0) R = - R
      F = R
      DO 160 I=1,NEA
  160  E(I) = - R
      DO 240 K=1,NEA
       D = E(K)
  170  T = 0.5D0 * (D + F)
      IF (DABS(D-F).GT.EPS2 .AND. T.NE.D .AND. T.NE.F) THEN
       J = 0
       I = 1
  180 Q = W(I,6) - T
  190 IF (Q .GE. 0.0) J = J + 1
      IF (Q .EQ. 0.0) GO TO 200
      I = I + 1
      IF (I .GT. N) GO TO 210
      Q = W(I,6) - T - W(I-1,2) / Q
      CALL OVERFL(L)
      IF (L .NE. 1) GO TO 190
      J = J + 1
      I = I - 1
  200 I = I + 2
      IF (I .LE. N) GO TO 180
  210 IF (NE .LT. 0) J = N - J
      IF (J .LT. K) THEN
      F = T
      GO TO 170
      ENDIF
      D = T
      M = MIN0(J, NEA)
      DO 230 I=K,M
  230  E(I) = T
      GO TO 170
      ENDIF
  240 E(K) = T
      ENDIF
C  EIGENVECTORS BY INVERSE ITERATION.
        IF (NV .EQ. 0) RETURN
         IF (N .EQ. 2) THEN
          W(1,6) = A(1,1)
          W(2,6) = A(2,2)
         ENDIF
C
        W(N,1) = 0.0D0
        MM = 584287
        DO 410 I=1,NVA
         DO 260 J=1,N
          W(J,2) = W(J,6) - E(I)
          W(J,3) = W(J,1)
  260     V(J,I) = 1.0D0
         SW = 0
C  REDUCE TO TRIANGULAR FORM
         DO 280 J=1,N-1
          IF (DABS(W(J,2)) .GE. DABS(W(J,1))) THEN
           IF (W(J,2) .EQ. 0.0) W(J,2) = 1.0D-30
           W(J,5) = W(J,1) / W(J,2)
           LW(J) = 0
           W(J+1,2) = W(J+1,2) - W(J,5) * W(J,3)
           W(J,4) = 0.0D0
          ELSE
           W(J,5) = W(J,2) / W(J,1)
           LW(J) = 1
           W(J,2) = W(J,1)
           T = W(J,3)
           W(J,3) = W(J+1,2)
           W(J,4) = W(J+1,3)
           W(J+1,2) = T - W(J,5) * W(J,3)
           W(J+1,3) = - W(J,5) * W(J,4)
          ENDIF
  280 CONTINUE
          IF (W(N,2) .EQ. 0.0) W(N,2) = 1.0D-30
C  BEGIN BACK SUBSTITUTION
         IF (I .NE. 1) THEN
          IF (DABS(E(I) - E(I-1)) .LT. EPS1) THEN
C  GENERATE RANDOM NUMBERS
           DO 290 J=1,N
            MM = MM * 48828125
  290       V(J,I) = DFLOAT(MM) * 0.4656613D-9
          ENDIF
        ENDIF
  300 CONTINUE
         T = V(N,I)
         R = V(N-1,I)
  310    V(N,I) = T / W(N,2)
         V(N-1,I) = (R - W(N-1,3) * V(N,I)) / W(N-1,2)
         CALL OVERFL(L)
         IF (L .EQ. 1) THEN
          DO 320 J=1,N-2
  320      V(J,I) = V(J,I) * 1.0D-5
          T = T * 1.0D-5
          R = R * 1.0D-5
          GO TO 310
         ENDIF 
         IF (N .NE. 2) THEN
         K = N - 2
  340    T = V(K,I)
  350    V(K,I) = (T - W(K,3)*V(K+1,I) - W(K,4)*V(K+2,I)) / W(K,2)
        CALL OVERFL(L)
        IF (L .EQ. 1) THEN
         DO 360 J=1,N
  360     V(J,I) = V(J,I) * 1.0D-5
         T = T * 1.0D-5
         GO TO 350
        ENDIF
        K = K - 1
        IF (K .GT. 0) GO TO 340
       ENDIF
        IF (SW .NE. 0) THEN
         SW = 1
         DO 400 J=1,N-1
          IF (LW(J) .NE. 0) THEN
           V(J+1,I) = V(J+1,I) - W(J,5) * V(J,I)
          ELSE
           T = V(J,I)
           V(J,I) = V(J+1,I)
           V(J+1,I) = T - W(J,5) * V(J+1,I)
          ENDIF
  400    CONTINUE
        GO TO 300
       ENDIF
  410 CONTINUE
C  BEGIN BACK TRANSFORMATION
      IF (N .NE. 2) THEN
       DO 415 I=1,N-2
  415   W(I,1) = - W(I,1) * A(I+1,I)
       DO 460 I=1,NVA
        K = N - 2
  420   R = W(K,1)
        IF (R .NE. 0.0) THEN
         R = 1.0D0 / R
         S = 0.0D0
         DO 430 J=K+1,N
  430     S = S + A(J,K) * V(J,I)
         R = R * S
         DO 440 J=K+1,N
  440     V(J,I) = V(J,I) - R * A(J,K)
        ENDIF
        K = K - 1
        IF (K .GE. 1) GO TO 420
  460 CONTINUE
C  NORMALIZE EIGENVECTORS.  MAX ELEMENT IS 1.
      ENDIF
      DO 490 I=1,NVA
       T = DABS(V(1,I))
       K = 1
       DO 480 J=2,N
        R = DABS(V(J,I))
        IF (T .LT. R) THEN
         T = R
         K = J
        ENDIF
  480  CONTINUE
       T = 1.0D0 / V(K,I)
       DO 490 J=1,N
  490   V(J,I) = V(J,I) * T
       IF (NV .LT. 0) RETURN
C  ORTHONORMALIZE AS NORM IS 1.
      DO 550 I=1,NVA
       IF (I .NE. 1) THEN
        IF (DABS(E(I) - E(I-1)) .LT. EPS1) THEN
         DO 510 J=M,I-1
          S = 0.0D0
          DO 500 K=1,N
  500      S = S + V(K,J) * V(K,I)
          DO 510 K=1,N
  510      V(K,I) = V(K,I) - S * V(K,J)
        ELSE
         M = I
        ENDIF
       ENDIF
       S = 0.0D0
       DO 540 J=1,N
  540   S = S + V(J,I) * V(J,I)
       T = 0.0D0
       IF (S .NE. 0.0) T = DSQRT(1.0D0 / S)
       DO 550 J=1,N
  550   V(J,I) = V(J,I) * T
       RETURN
       END
C
       SUBROUTINE OVERFL(L)
       L = 2
       RETURN
       END
