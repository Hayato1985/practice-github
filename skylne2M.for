C ==================================================================================================
C プログラムの構成
C SKYLNE2M.FOR：スカイライン法による連立一次方程式の解法に関するプログラム、（CMEM7.FOR仕様）
C   FIXNO  NET系の拘束節点情報を、SKYLNE型に変更する
C   SIZING 全体剛性マトリクスの成分と１次元配列化されたスカイラインマトリクスの対応を示す指標配列を作成する
C   MATRIX 要素剛性マトリクスの成分を全体剛性マトリクス（１次元配列化されたスカイラインマトリクス）へ代入する
C   BCONUU 係数マトリックスおよび右辺ベクトルに対して、幾何学的な条件を処理する
C   RARNG  各種変位ベクトルの順序を元に戻す
C   SKYLNE スカイライン法により連立一次方程式の解を求める
C ======================================================================
      SUBROUTINE FIXNO(IFFIX,IFPRE,NODEM,NODEC,
     &                 MTOTV,NDOFN,NELEM,NELEC,NFPOIN,
     &                 NFREE,NOFIX,NORDR,NPOIN,NTOTV,NVFIX,
     &                 NELEF,NODEF)
C ======================================================================
C *** CMEM型の拘束節点情報を、SKYLNE型に変更する
C  * CMEM型 
C     NFPOIN            ：総自由節点数
C     NPOIN             ：総節点数
C  * SKYLNE型
C     LNODS(IELEM,INODE)：要素IELEMの要素内節点INODEの全体節点番号
C     NELEM             ：要素数
C     IFFIX(ITOTV)      ：拘束に関する情報（0:自由、1:拘束)
C     IFPRE(IVFIX)      ：拘束状態を表す変数（例　101:x,z方向拘束、y方向自由）
C     NDOFN             ：１節点当りの自由度数
C     NFREE             ：総自由度数（拘束を除く)
C     NOFIX(IVFIX)      ：拘束節点の番号
C     NTOTV             ：総自由度数（拘束を含む), ≦MTOTV
C     NVFIX             ：拘束節点の総数
C     
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IFFIX(MTOTV)
      DIMENSION NOFIX(NVFIX)
      DIMENSION NORDR(MTOTV)
	DIMENSION IFPRE(NVFIX)
      DIMENSION NODEM(NELEM,4)       ! 膜要素のデータ（1～3：節点番号、4：部材番号）
      DIMENSION NODEC(NELEC,3)     !ケーブル要素の節点番号(1～2)、部材番号(3)
      DIMENSION NODEF(NELEF,3)       ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** FIXNO **********'
C
C *** 拘束情報を示す配列IFFIXの設定、部材に接続されていない節点は拘束節点IIFIX=1として処理する
C
      DO 10 ITOTV=1,NTOTV
        IFFIX(ITOTV)=1
   10 CONTINUE
      DO 20 IELEM=1,NELEM
      DO 20 INODE=1,3
        IPOIN=NODEM(IELEM,INODE)
        DO 20 IDOFN=1,3
	    ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IFFIX(ITOTV)=0
   20 CONTINUE
      DO 30 IELEC=1,NELEC
      DO 30 INODE=1,2
        IPOIN=NODEC(IELEC,INODE)
        DO 30 IDOFN=1,3
	    ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IFFIX(ITOTV)=0
   30 CONTINUE
C
      DO 40 IELEF=1,NELEF
	DO 40 INODE=1,2
        IPOIN=NODEF(IELEF,INODE)
	  ITTV1=(IPOIN-1)*NDOFN+1
	  ITTV2=IPOIN*NDOFN
	  DO 40 ITOTV=ITTV1,ITTV2
	    IFFIX(ITOTV)=0
   40 CONTINUE
C
      DO 100 IVFIX=1,NVFIX
        NLOCA=(NOFIX(IVFIX)-1)*NDOFN
        IFDOF=10**(NDOFN-1)
        DO 100 IDOFN=1,NDOFN
          ITOTV=NLOCA+IDOFN
          IF(IFPRE(IVFIX).GE.IFDOF)THEN
            IFFIX(ITOTV)=1
            IFPRE(IVFIX)=IFPRE(IVFIX)-IFDOF
          ENDIF
          IFDOF=IFDOF/10
  100 CONTINUE
C
C *** 節点番号 NFPOIN+1 ～ NPOIN の自由度をすべて拘束する
C
      DO 400 IPOIN=NFPOIN+1,NPOIN
      DO 400 IDOFN=1,NDOFN
        ITOTV=(IPOIN-1)*NDOFN+IDOFN
        IFFIX(ITOTV)=1
  400 CONTINUE
C
C *** 総自由度数(拘束を除く)NFREEの計算
      NFREE=NTOTV
      DO 500 ITOTV=1,NTOTV
        NFREE=NFREE-IFFIX(ITOTV)
  500 CONTINUE
C
C *** NORDR(ITOTV):全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて
C                  並換える時の新しい全体自由度番号の配列の作成
      IFREE=0
      IFIX=NFREE
      DO 9 ITOTV=1,NTOTV
        IF(IFFIX(ITOTV).EQ.0)THEN
          IFREE=IFREE+1
          NORDR(ITOTV)=IFREE
        ELSE
          IFIX=IFIX+1
          NORDR(ITOTV)=IFIX
        ENDIF
    9 CONTINUE
C
      IF(MCHK.EQ.1)THEN
        CALL WMTI1(2,IFFIX,MTOTV,0,NTOTV,0,'IFFIX','FIXNO',0)
        CALL WMTI1(2,NOFIX,NVFIX,0,NVFIX,0,'NOFIX','NVFIX',NVFIX)
        CALL WMTI1(2,NORDR,MTOTV,0,NTOTV,0,'NORDR','FIXNO',0)
      ENDIF
      RETURN
	END
C
C
C
C ======================================================================
      SUBROUTINE SIZING(NODEM,NODEC,MSIZE,MTOTV,NELEM,NELEC,
     &                  NDOFN,NFREE,NORDR,NPOIN,NTOTV,NWDTH,NWSUM,
     &                  NELEF,NODEF)
C ======================================================================
C *** 全体剛性マトリクスの成分と１次元配列化されたスカイラインマトリクスの対応を示す指標配列を作成する。
C     NODEM(IELEM,INODE) :膜要素IELEMの要素内節点INODEの全体節点番号
C     NODEC(IELEC,INODE):ケーブル要素IELECの要素内節点INODEの全体節点番号
C     NTOTV       : 総自由度数(=NPOIN*NDOFN<MTOTV)
C     NELEM       : 膜要素数
C     NELEC       : ケーブル要素数
C     NDOFN       : １節点あたりの自由度数
C     NPOIN       : 節点数
C     NORDR(ITOTV): 全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて並換える時の新しい全体自由度番号
C     NWDTH(ITOTV): 全体剛性マトリクスの第ITOTV行の第１非ゼロ要素から対角要素までの個数
C     NWSUM(ITOTV): 全体剛性マトリクスの第１～ITOTV行の第１非ゼロ要素から対角要素までの要素数の累和
C     LEFTS(ITOTV): 全体剛性マトリクスの第IPOIN行目の第１非ゼロ要素の列番号
C
      PARAMETER(MCHK=0,MMTOTV=20000)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NODEM(NELEM,4)       ! 膜要素のデータ（1～3：節点番号、4：部材番号）
      DIMENSION NODEC(NELEC,3)       !ケーブル要素の節点番号(1～2)、部材番号(3)
      DIMENSION NORDR(MTOTV)
      DIMENSION NWDTH(MTOTV)
      DIMENSION NWSUM(0:MTOTV)
      DIMENSION LEFTS(MMTOTV)
      DIMENSION NODEF(NELEF,3)       ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
      IF(MMTOTV.LT.NTOTV)THEN
        WRITE(8,9002)MMTOTV,NTOTV
 9002   FORMAT(' (SUB SIZING) MMTOTV:',I6,' IS LESS THAN NTOTV:',I6)
        STOP
      ENDIF
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** SIZING **********'
C
      DO 100 ITOTV=1,NTOTV
        LEFTS(ITOTV)=ITOTV
  100 CONTINUE
C
      DO 200 IELEM=1,NELEM
        DO 200 INODE=1,3
          IPOIN=NODEM(IELEM,INODE)
          DO 200 JNODE=1,3
            JPOIN=NODEM(IELEM,JNODE)
            DO 200 IDOFN=1,3
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:並べ換え後の行番号
              DO 200 JDOFN=1,3                                          JTOTV:並べ換え後の列番号
                JTOTV=NORDR((JPOIN-1)*NDOFN+JDOFN)
                IF(LEFTS(ITOTV).GT.JTOTV)LEFTS(ITOTV)=JTOTV
  200 CONTINUE
C
      DO 250 IELEC=1,NELEC
        DO 250 INODE=1,2
          IPOIN=NODEC(IELEC,INODE)
          DO 250 JNODE=1,2
            JPOIN=NODEC(IELEC,JNODE)
            DO 250 IDOFN=1,3
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:並べ換え後の行番号
              DO 250 JDOFN=1,3                                          JTOTV:並べ換え後の列番号
                JTOTV=NORDR((JPOIN-1)*NDOFN+JDOFN)
                IF(LEFTS(ITOTV).GT.JTOTV)LEFTS(ITOTV)=JTOTV
  250 CONTINUE
C
      DO 270 IELEF=1,NELEF
        DO 270 INODE=1,2
          IPOIN=NODEF(IELEF,INODE)
          DO 270 JNODE=1,2
            JPOIN=NODEF(IELEF,JNODE)
            DO 270 IDOFN=1,6
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:並べ換え後の行番号
              DO 270 JDOFN=1,6                                          JTOTV:並べ換え後の列番号
                JTOTV=NORDR((JPOIN-1)*NDOFN+JDOFN)
                IF(LEFTS(ITOTV).GT.JTOTV)LEFTS(ITOTV)=JTOTV
  270 CONTINUE
C
      DO 300 IPOIN=1,NPOIN
      DO 300 IDOFN=1,NDOFN
        ITOTV=(IPOIN-1)*NDOFN+IDOFN                                     ITOTV:並べ換え前の行番号
        JTOTV=NORDR(ITOTV)                                              JTOTV:並べ換え後の行番号
        NWDTH(JTOTV)=JTOTV-LEFTS(JTOTV)+1
  300 CONTINUE
C
      NWSUM(0)=0
      DO 400 JTOTV=1,NTOTV
        NWSUM(JTOTV)=NWSUM(JTOTV-1)+NWDTH(JTOTV)
  400 CONTINUE
C
      WRITE(8,9001)NTOTV,NFREE,NWSUM(NFREE),NWSUM(NTOTV)
      IF(MCHK.EQ.1)THEN
        CALL WMTI1(2,LEFTS,MMTOTV,0,NTOTV,0,'LEFTS','SIZNG',0)
        CALL WMTI1(2,NWDTH,MTOTV,0,NTOTV,0,'NWDTH','SIZNG',0)
C        CALL WMTI1(2,NWSUM,MTOTV+1,0,NTOTV+1,0,'NWSUM','SIZNG',0)
 9001 FORMAT(/,' *** NTOTV=',I4,'  NFREE=',I4,
     $            '  NFSIZE=',I6,'  NSIZE=',I6)
      ENDIF
C
      NFSIZE=NWSUM(NFREE)
      MSIZE=NFSIZE    !動的配列のサイズ、TSTIF(MSIZE)
!
      IF(NFSIZE.GT.MSIZE)THEN
        WRITE(8,9003)NFSIZE,MSIZE
        STOP
 9003 FORMAT(/,'  *** NFSIZE(',I6,') EXCEEDS THE LIMIT MSIZE(',I6,')')
      ENDIF
C  
      RETURN
      END
C
C
C
C ======================================================================
      SUBROUTINE MATRIX(ESTIF,IELEM,LNODS,NELEM,MSIZE,NDOFN,NDOF1,NNODE,
     &                  NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
C ======================================================================
C *** 要素剛性マトリクスの成分を全体剛性マトリクス（１次元配列化されたスカイラインマトリクス）へ代入する。
C     ESTIF(IEVAB,JEVAB) : 要素剛性マトリクス(IEVAB,JEVAB=1～NDOF1*NNODE)
C     TSTIF(JPOSI)       : 全体剛性マトリクス（スカイラインマトリクス）
C     LNODS(IELEM,INODE) : 要素IELEMの要素内節点INODEの全体節点番号
C     IELEM   : 要素番号
C     NDOFN   : １節点あたりの自由度数 PARAMETER(NDOFN=6)
C     NDOF1   : 各要素剛性マトリクスにおいて定義されている１節点あたりの自由度数
C     NNODE   : 各要素剛性マトリクスにおいて定義されている１要素あたりの節点数
C     NPOIN   : 節点数
C     NSIZE   : スカイラインマトリクスのサイズ
C     NWSUM(ITOTV)：全体剛性マトリクス、第１～ITOTV行の第１非ゼロ要素から対角要素までの要素数の累和
C     IFFIX(ITOTV)：拘束に関する情報（0:自由、1:拘束)、ITOTV：並べ替え前の全体自由度番号
C
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ESTIF(NDOF1*NNODE,NDOF1*NNODE)
      DIMENSION LNODS(NELEM,NNODE)
      DIMENSION NWSUM(0:MTOTV)
      DIMENSION NORDR(MTOTV)
      DIMENSION TSTIF(MSIZE)
      DIMENSION IFFIX(MTOTV)
C
C *** IEVAB:要素剛性マトリクスの行番号
C     JEVAB:要素剛性マトリクスの列番号
C     ITOTV:並べ換え後の全体剛性マトリクスの行番号
C     JTOTV:並べ換え後の全体剛性マトリクスの列番号
C     IITTV:並べ換え前の全体剛性マトリクスの行番号
C     JJTTV:並べ換え前の全体剛性マトリクスの列番号
      DO 100 INODE=1,NNODE
      IPOIN=LNODS(IELEM,INODE)
      DO 100 JNODE=1,INODE
      JPOIN=LNODS(IELEM,JNODE)
         DO 100 IDOFN=1,NDOF1
         IEVAB=(INODE-1)*NDOF1+IDOFN
         IITTV=(IPOIN-1)*NDOFN+IDOFN
         ITOTV=NORDR(IITTV)                                     
         DO 100 JDOFN=1,NDOF1
         JEVAB=(JNODE-1)*NDOF1+JDOFN
         JJTTV=(JPOIN-1)*NDOFN+JDOFN
         JTOTV=NORDR(JJTTV)                                     
         IJTTV=IFFIX(IITTV)+IFFIX(JJTTV)
	     IF(JEVAB.LE.IEVAB)THEN
              IF(JTOTV.LE.ITOTV)THEN
                JPOSI=NWSUM(ITOTV)-ITOTV+JTOTV
	        ELSE
                JPOSI=NWSUM(JTOTV)-JTOTV+ITOTV
	        ENDIF
              IF(IJTTV.EQ.0)THEN
                TSTIF(JPOSI)=TSTIF(JPOSI)+ESTIF(IEVAB,JEVAB)
              ENDIF
	     ENDIF
  100 CONTINUE
      IF(MCHK.EQ.1)THEN
        MEVAB=NDOFN*NNODE
        NSIZE=1872
        CALL WMTR2(2,ESTIF,MEVAB,MEVAB,MEVAB,MEVAB,'ESTIF','MATRX',0)
        CALL WMTR1(2,TSTIF,MSIZE,0,NSIZE,0,'TSTIF','MATRX',0)
      ENDIF
      RETURN
      STOP
      END
C
C
C
C ======================================================================
      SUBROUTINE BCONUU(ELOAD,IFFIX,MSIZE,MTOTV,NDOFN,NFREE,
     &                  NOFIX,NORDR,NTOTV,NVFIX,NWSUM,PRESC,TSTIF)
C ======================================================================
C *** 係数マトリックスおよび右辺ベクトルに対して、幾何学的な条件を処理する
C     ELOAD(ITOTV)：残差力ベクトル
C     IFFIX(ITOTV)：拘束に関する情報（0:自由、1:拘束)
C     MSIZE       ：スカイライン･マトリクスのサイズ最大値
C     NDOFN       ：１節点当りの自由度数
C     NFREE       ：総自由度数
C     NOFIX(IVFIX)：拘束節点の番号
C     NTOTV       ：=NPOIN*NDOFN, ≦MTOTV
C     NVFIX       ：拘束節点の総数
C     NWSUM(ITOTV):全体剛性マトリクス、第１～ITOTV行の第１非ゼロ要素から対角要素までの要素数の累和
C     PRESC(IVFIX,IDOFN):IVFIX番目の拘束節点のIDOFN番目の自由度の拘束量
C     TSTIF(JPOSI): 全体剛性マトリクス（スカイラインマトリクス）
C *L  ET(IMTOTV)  ：作業用配列
C
      PARAMETER(MCHK=0,MMTOTV=20000)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IFFIX(MTOTV)
      DIMENSION NORDR(MTOTV)
      DIMENSION ELOAD(MTOTV)
      DIMENSION NOFIX(NVFIX)
      DIMENSION NWSUM(MTOTV)
      DIMENSION PRESC(NVFIX,NDOFN)
      DIMENSION TSTIF(MSIZE)
	DIMENSION ET(MMTOTV)
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** BCONUU **********'
	IF(MMTOTV.LT.NTOTV)THEN
        WRITE(8,*)' *** ERROR (SUB BCONUU) MMTOTV<NTOTV (',NTOTV,' )'
        STOP
	ENDIF
C
C *** ITOTV:並べ換え前の行番号、JTOTV:並べ換え後の行番号
C
	DO 100 ITOTV=1,NTOTV
	  JTOTV=NORDR(ITOTV)
        ET(JTOTV)=ELOAD(ITOTV)                                          残差力ベクトルの符号に注意！
  100 CONTINUE
C
C *** ELOAD(NFREE+1～NTOTVC)に強制変位量を入れておく。
C     SUBROUTINE SKYLINEにてELOAD(1～NFREE)に各自由度の変位が入る。
C
      DO 200 IVFIX=1,NVFIX
	  IPOIN=NOFIX(IVFIX)
	  DO 200 IDOFN=1,NDOFN
	    ITOTV=(IPOIN-1)*NDOFN+IDOFN
	    JTOTV=NORDR(ITOTV)
	    IF(IFFIX(ITOTV).NE.0)ET(JTOTV)=PRESC(IVFIX,IDOFN)
  200 CONTINUE
C
C *** ITOTV:並べ換え後の行番号順で増分荷重ELOADを作成する。
C
        DO 500 ITOTV=1,NTOTV
	    ELOAD(ITOTV)=ET(ITOTV)
  500   CONTINUE
      IF(MCHK.EQ.1)THEN
        NFSIZE=NWSUM(NFREE)
        CALL WMTR1(2,ELOAD,MTOTV,0,NTOTV,0,'ELOAD','BCONU',0)
      ENDIF
      RETURN
      END
C
C
C
C ======================================================================
      SUBROUTINE RARNG(ASDIS,ELOAD,MTOTV,NORDR,NTOTV)
C ======================================================================
C *** 各種変位ベクトルの順序を元に戻す
C     ASDIS(ITOTV)：増分変位
C     ELOAD(ITOTV)：残差力ベクトル
C     NORDR(ITOTV): 全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて並換える時の新しい全体自由度番号
C     NTOTV       ：=NPOIN*NDOFN, ≦MTOTV
C  *L ASD(IMTOTV) ：作業用配列
C
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASDIS(MTOTV)
      DIMENSION ELOAD(MTOTV)
      DIMENSION NORDR(MTOTV)
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** RARNG **********'
C
C *** ITOTV:並べ換え前の行番号、JTOTV:並べ換え後の行番号
C
      DO 100 ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        ASDIS(ITOTV)=ELOAD(JTOTV)
  100 CONTINUE
      RETURN
      END
C
C
C ======================================================================
      SUBROUTINE SKYLNE(INDEX,A, B, NSUM, LX, M, MX, IER, D, T, NWK)
C ======================================================================
C *** スカイライン法により連立一次方程式の解を求める
************************************************************************
C *** CALL SKYLNE(TSTIF,P,NWSUM,MSIZE,NFREE,MTOTV,IER)
*  SKYLINE METHOD FOR FINITE ELEMENT METHOD TO SOLVE AX = B.           *
*  PARAMETERS                                                          *
*    (1) A: 1-DIM. ARRAY CONTAINING THE SKYLINE MATRIX                 *
*    (2) B: 1-DIM. ARRAY CONTAINING THE RIGHT HAND VECTOR              *
*    (3) NSUM: 1-DIM. ARRAY CONTAINING ACCUMURATED SUM OF NON-ZERO     *
*              ELEMENTS ON EACH ROW                                    *
*    (4) D: 1-DIM. WORKING ARRAY                                       *
*    (5) L: SIZE OF THE ARRAY (A)                                      *
*    (6) M: ROW SIZE OF THE MATRIX (A)                                 *
*    (7) T: 1-DIM. WORKING ARRAY                                       *
*    (8) NWK: 1-DIM. WORKING ARRAY                                     *
*    (9) EPS: TOLERANCE FOR PIVOTAL ELEMENTS                           *
*   (10) IER: ERROR CODE                                               *
*  COPYRIGHT   T. OGUNI    JUNE 30 1989    VERSION 1.0                 *
************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER*4 FST, FSTI
       PARAMETER(EPS=1.0D-28)
       DIMENSION A(LX), B(MX), NSUM(0:MX)
       DIMENSION D(MX), T(MX), NWK(MX)
C
C      IF(MX.GT.MXX)THEN
C        WRITE(8,9001)MX,MXX
C        STOP
C 9001 FORMAT(/,'  *** MX(',I6,') EXCEEDS THE LIMIT MXX(',I6,')')
C      ENDIF
C        CALL WMTR1(2,A,LX,0,M,0,'A    ','TSTIF',0)
C        CALL WMTR1(2,B,MX,0,M,0,'B    ','P    ',0)
C
! *** INDEX=1：逆マトリクスを求める
      IF(INDEX.EQ.1)THEN
!
      D(1) = 1.0D0 / A(1)
      DO 5 I=1,M
    5  NWK(I) = I - (NSUM(I) - NSUM(I-1)) + 1
      DO 45 K=2,M
       FST = NWK(K)
       IF (K .GT. 2) THEN
        II = 2
        IF (FST .GT. 2) II = FST
        DO 30 I=II,K-1
         FSTI = NWK(I)
         LGTH = MAX(FST, FSTI)
         IF (LGTH .LE. I) THEN
          IK = NSUM(I) - I
          KK = NSUM(K) - K
          IF (FST .NE. I) THEN
           IF (LGTH .LT. I) THEN
            S = 0.0D0
            DO 20 J=LGTH,I-1
   20        S = S + A(IK+J) * A(KK+J)
            A(KK+I) = A(KK+I) - S
           ENDIF
          ENDIF
         ENDIF
   30   CONTINUE
       ENDIF
C
       S = 0.0D0
       II = NSUM(K-1) - FST + 1
       IF (FST .NE. K) THEN
        DO 35 I=FST,K-1
         T(I) = A(II+I)
   35    A(II+I) = D(I) * T(I)
        DO 40 I=FST,K-1
   40    S = S + A(II+I) * T(I)
        D(K) = 1.0D0 / (A(NSUM(K)) - S)
        IF (D(K) .LE. EPS) THEN
         WRITE(6,*) '(SUBR. SKYLNE) DIAGONAL ELEMENT IS SMALL.',K,D(K)
         WRITE(11,*) '(SUBR. SKYLNE) DIAGONAL ELEMENT IS SMALL.',K,D(K)
         IER = 2
C         RETURN
        ENDIF
        A(NSUM(K)) = 1.0D0 / D(K)
       ELSE
        D(K) = 1.0D0 / A(NSUM(K))
       ENDIF
   45 CONTINUE
!
      ENDIF
!
C  SKYLINE SUBSTITUTION
C     ENTRY SKYSUB(B)
C  FORWARD SUBSTITUTION
      DO 60 K=2,M
       FST = NWK(K)
       KK = NSUM(K-1) - FST + 1
       DO 50 J=FST,K-1
   50   B(K) = B(K) - A(KK+J) * B(J)
   60 CONTINUE
C  BACKWARD SUBSTITUTION
      DO 70 I=1,M
   70  B(I) = B(I) * D(I)
      DO 90 K=M,2,-1
       FST = NWK(K)
       KK = NSUM(K-1) - FST + 1
       IF (FST .NE. K) THEN
        DO 80 J=FST,K-1
   80    B(J) = B(J) - A(KK+J) * B(K)
       ENDIF
   90 CONTINUE
      IER = 0
      RETURN
      END
