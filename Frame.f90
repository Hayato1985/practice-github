! =====================================================================
      SUBROUTINE FRAME(INDEX,CALFR,COORD,VALK,&
                       NPOIN,NDOFN,NELEF,NODEF,NPIN,NPFRM,  &
                       NVARF,PROPF,STRSF,STRF1,                   &
                       BETA,IFFIX,MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
! =====================================================================
! *** ＜曲げ要素の入出力、剛性マトリクス作成に関するサブルーチン＞ ***
!     NELEF：曲げ要素数
!     NPFRM：初期張力（=0：なし、=1：あり）
!     NVARF：曲げ要素の種類数
!
      PARAMETER( MFORMAT=1 )    ! MFORMAT：固定フォーマット（=0）、配列データの一部はフリーフォーマット（=1）
!
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION CALFR(NELEF)              ! CALFR(IELEF)：部材の長さ
      DIMENSION COORD(NPOIN,3)            ! COORD(IPOIN,IDIR)：節点座標
      DIMENSION VALK(NELEF,2)             ! VALK(IELEF, )：材端1,2の接合条件（＝1.0E+20：剛接、＝0.0：ピン接、その他：半剛接）
      DIMENSION NODEF(NELEF,3)            ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
      DIMENSION NPIN(NELEF)               ! NPIN(IELEF)：部材端部の接合条件（=0：剛剛、=1：ピン剛、=2：剛ピン、=3：ピンピン）
      DIMENSION STRSF(NELEF,NDOFN*2)      ! STRSF(IELEF,IDOFN)：曲げ要素の応力（基準座標系における部材端荷重)
      DIMENSION STRF1(NELEF,NDOFN*2)      ! STRSF(IELEF,IDOFN)：曲げ要素の応力（部材座標系における部材端荷重)
      DIMENSION BETA(NELEF)               ! BETA(IELEF)：曲げ部材の回転角
      DIMENSION IFFIX(MTOTV)              ! IFFIX(ITOTV)：拘束に関する情報（0:自由、1:拘束)、ITOTV：並べ替え前の全体自由度番号
      DIMENSION NORDR(MTOTV)              ! NORDR(ITOTV): 全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて並換える時の新しい全体自由度番号
      DIMENSION NWSUM(0:MTOTV)            ! NWSUM(ITOTV):全体剛性マトリクス、第１～ITOTV行の第１非ゼロ要素から対角要素までの要素数の累和
      DIMENSION PROPF(NVARF,10)           ! 曲げ要素の材料特性（種類番号、特性項目(1～10)）
      DIMENSION TSTIF(MSIZE)              ! TSTIF(JPOSI): 全体剛性マトリクス（スカイラインマトリクス）
      DIMENSION DIS(NPOIN*NDOFN)          ! DIS(IPOIN*IDOFN)：変位増分
      DIMENSION Q(NPOIN*NDOFN)            ! Q(ITOTV)：残差力ベクトル
      DIMENSION C(12,6),CSK(12,6),CT(6,12),SK(6,6),SKP(12,12)
      DIMENSION C0(12,6),CSK0(12,6),SKP0(12,12)
      DIMENSION SKG(12,12),SKGT(12,12),TSKGT(12,12)
      DIMENSION TR(12,12),TRT(12,12),SKK(12,12)
      DIMENSION DM(12,1),PM(12,1)
      DIMENSION VECLY(3)                                                ! VECLY(3)：局所座標系のｙ軸方向を表わすベクトル
character*170 text

      NEVAB=NDOFN*2
!
! ---------------------------------------------------------------------
      IF(INDEX.EQ.1)THEN ! *** ＜曲げ要素に関するデータをファイル＃７から読込む＞ ***
! ---------------------------------------------------------------------

! *** ＜曲げ要素の重量、剛性、せん断補正係数、部材主軸の法線ベクトル＞

do i=1,6
  read(7,'(A)')text
  write(8,'(A)')text
enddo
do i=1,NVARF
  read(7,*,err=999)num,(PROPF(i,j),j=1,10)
enddo
write(8,8002)(i,(PROPF(i,j),j=1,10),i=1,NVARF)
!read(7,*,err=999)(num,PROPF(i,7),(PROPF(i,j),j=1,6),(PROPF(i,j),j=8,10),i=1,NVARF)
!write(8,8002)(i,PROPF(i,7),(PROPF(i,j),j=1,6),(PROPF(i,j),j=8,10),i=1,NVARF)
8002 format(I15,F15.4,5G15.4,4F15.4)

! *** ＜曲げ要素の節点番号、材端バネの係数＞

do i=1,6
  read(7,'(A)')text
  write(8,'(A)')text
enddo
read(7,*,err=999)(num,NODEF(i,1),NODEF(i,2),VALK(i,1),VALK(i,2),NODEF(i,3),i=1,NELEF)
WSUM=0.0
do IE=1,NELEF
  CAL=0.0
  NOD1=NODEF(IE,1)
  NOD2=NODEF(IE,2)
  IV=NODEF(IE,3)
  do IDIR=1,3
    CAL=CAL+(COORD(NOD2,IDIR)-COORD(NOD1,IDIR))**2
  enddo
  CALFR(IE)=SQRT(CAL)
!!!!!!単位系変更  WEIGHT=PROPF(IV,7)*CALFR(IE)  !重量[kgf]＝単位重量[kg/m]*長さ[m]
  WEIGHT=PROPF(IV,7)*CALFR(IE)*9.8  !重量[N]＝単位重量[kg/m]*長さ[m]*9.8
  WSUM=WSUM+WEIGHT
  write(8,8004)IE,NODEF(IE,1),NODEF(IE,2),VALK(IE,1),VALK(IE,2),NODEF(IE,3),CALFR(IE),WEIGHT
enddo
8003 format(//&
'曲げ要素',/&
'       要素番号       節点番号                材端のバネ係数                  部材種類番号       部材長さ           重量',/&
'                           I端            J端            I端            J端                           [m]          [kgf]',/&
'     I=1～NELEF     NODEF(I,1)     NODEF(I,2)      VALK(I,1)      VALK(I,2)     NODEF(I,3)')
8004 format(3I15,2G15.4,I15,2F15.4)
WRITE(8,8005) WSUM
8005 FORMAT(4X,120('-'),/,85X,'  曲げ部材の総重量 (',F15.4,' [N])')

        CALL ZEROR2(STRSF,NELEF,NEVAB,NELEF,NEVAB)
        IF(NPFRM.EQ.1)THEN
          WRITE(8,8600)
          DO IE=1,NELEF
!           !READ(7,7600,err=999)KK,(STRSF(IE,IEVAB),IEVAB=1,NEVAB)
            !READ(7,*,err=999)KK,(STRSF(IE,IEVAB),IEVAB=1,NEVAB)
            WRITE(8,7600)IE,(STRSF(IE,IEVAB),IEVAB=1,NEVAB)
          ENDDO
        ENDIF
 7600   FORMAT(I5,12F15.4)
 8600   FORMAT(/,' (PRETENSION OF FRAME ELEMENTS)',/, &
               '  ELEM      PX1       PY1       PZ1      Myz1      Mxz1      Mxy1',&
                    '       PX2       PY2       PZ2      Myz2      Mxz2      Mxy2')
!
        CALL ZEROR1(BETA,NELEF,NELEF)
!
! ---------------------------------------------------------------------
      ELSEIF(INDEX.EQ.2)THEN ! *** ＜曲げ要素の剛性マトリクスを作成する＞ ***
! ---------------------------------------------------------------------
        DO 800 NE=1,NELEF
          IPOIN=NODEF(NE,1)
          JPOIN=NODEF(NE,2)
          IV=NODEF(NE,3)
          ITOTV=(IPOIN-1)*NDOFN
          JTOTV=(JPOIN-1)*NDOFN
          DELX=COORD(JPOIN,1)-COORD(IPOIN,1)
          DELY=COORD(JPOIN,2)-COORD(IPOIN,2)
          DELZ=COORD(JPOIN,3)-COORD(IPOIN,3)
          CAL1=DSQRT(DELX**2+DELY**2+DELZ**2)
          CAL2=DSQRT(DELX**2+DELY**2) 
          DBETAX=(DIS(ITOTV+1)+DIS(JTOTV+1))*DELX/CAL1
          DBETAY=(DIS(ITOTV+2)+DIS(JTOTV+2))*DELY/CAL1
          DBETAZ=(DIS(ITOTV+3)+DIS(JTOTV+3))*DELZ/CAL1
          DBETA=(DBETAX+DBETAY+DBETAZ)/2
!          BETA(NE)=BETA(NE)+DBETA   

EA=PROPF(IV,1)
GA=PROPF(IV,2)
EIy=PROPF(IV,3)
EIz=PROPF(IV,4)
GIp=PROPF(IV,5)
VALKAPPA=PROPF(IV,6)
!
!  *** ＜部材座標系における微小変位剛性マトリクスの作成＞ ***
!
          CALL STIFF(SK,EA,EIz,EIy,GIp,GA,CAL1,VALKAPPA,VALK(NE,1),VALK(NE,2))   
!         CALL WMTR2(2,SK,6,6,6,6,'SK   ','NE   ',NE)
!
!  *** ＜マトリクスＣ0，Ｃt，座標変換マトリクスＴＲの作成＞ ***
!
          DO IDIR=1,3
            VECLY(IDIR)=PROPF(IV,IDIR+7)    ! 主軸法線ベクトル
          ENDDO
          CALL TRANST(TRT,CAL1,CAL2,DELX,DELY,DELZ,BETA(NE),VECLY)
          CALL TRANS(12,12,TRT,TR)
          CALL CNECT1(C0,CT,CAL1,CAL2,DELX,DELY,DELZ,BETA(NE),TRT)
          CALL TRANS(6,12,CT,C)
!
!  *** ＜弾性剛性マトリクスＳＫＰの作成＞ ***
!
          CALL MLTPLY(12,6,6,C,SK,CSK)
          CALL MLTPLY(12,6,12,CSK,CT,SKP)
!
!  *** ＜部材座標系における部材端荷重算定のためのマトリクスＳＫＰ０の作成＞ ***
!
          CALL MLTPLY(12,6,6,C0,SK,CSK0)
          CALL MLTPLY(12,6,12,CSK0,CT,SKP0)
!
! *** ＜曲げ要素の部材応力を計算する＞ ***
!
          DO INODE=1,2
            IPOIN=NODEF(NE,INODE)
            DO IDOFN=1,6
              ITOTV=(IPOIN-1)*NDOFN+IDOFN
              IEVAB=(INODE-1)*NDOFN+IDOFN
              DM(IEVAB,1)=DIS(ITOTV)
            ENDDO
          ENDDO
!
          CALL MLTPLY(12,12,1,SKP0,DM,PM)
!
           DO IEVAB=1,NEVAB
            STRSF(NE,IEVAB)=STRSF(NE,IEVAB)+PM(IEVAB,1)
          ENDDO
!
!  *** ＜部材座標系における幾何剛性マトリクスの作成＞ ***
!
          CALL GEOMAT(NELEF,NDOFN,NE,SKG,CAL1,STRSF,EIy,EIz,EA)
          CALL MLTPLY(12,12,12,SKG,TRT,SKGT)
          CALL MLTPLY(12,12,12,TR,SKGT,TSKGT)
!
          DO I=1,12
            DO J=1,12
              SKK(I,J)=SKP(I,J)+TSKGT(I,J)
            ENDDO
          ENDDO

!         CALL WMTR2(2,C0,12,6,12,6,'C0   ','CNCT1',1)
!         CALL WMTR2(2,C,12,6,12,6,'C    ','NE   ',NE)
!         CALL WMTR2(2,CT,6,12,6,12,'CT   ','NE   ',NE)
!         CALL WMTR2(2,SK,6,6,6,6,'SK   ','NE   ',NE)
!         CALL WMTR2(2,SKP,12,12,12,12,'SKP  ','NE   ',NE)
!         CALL WMTR2(2,TSKGT,12,12,12,12,'TSKGT','NE   ',NE)
!         CALL WMTR2(2,SKK,12,12,12,12,'SKK  ','NE   ',NE)
!
! *** ＜要素剛性マトリックス→全体剛性マトリックス＞ ***
!
          CALL MATRIX(SKK,NE,NODEF,NELEF,MSIZE,NDOFN,6,2,NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
!
! *** ＜曲げ部材の部材端荷重ベクトルを部材座標系ＳＴＲＦＲから基準座標系ＳＴＲＦ１に変換する＞ ***
!
          DO IDOFN=1,NDOFN
            DO INODE=1,2
              IPOIN=NODEF(NE,INODE)
              IEVAB=(INODE-1)*NDOFN+IDOFN
              PM(IEVAB,1)=STRSF(NE,IEVAB)
            ENDDO
          ENDDO
          CALL MLTPLY(12,12,1,TR,PM,DM)
          DO IEVAB=1,NEVAB
            STRF1(NE,IEVAB)=DM(IEVAB,1)
          ENDDO
!
! *** ＜曲げ部材応力の等価節点力ベクトルDMを全体の節点力ベクトルＱに加える＞ ***
!
          DO INODE=1,2
            IPOIN=NODEF(NE,INODE)
            DO IDOFN=1,6
              ITOTV=(IPOIN-1)*NDOFN+IDOFN
              IEVAB=(INODE-1)*NDOFN+IDOFN
              Q(ITOTV)=Q(ITOTV)+STRF1(NE,IEVAB)
            ENDDO
  IEVAB1=(INODE-1)*NDOFN+1
  IEVAB2=(INODE-1)*NDOFN+NDOFN
!  IF(NE.EQ.1.AND.INODE.EQ.1)WRITE(8,*)' (STRESS OF FRAME ELEMENTS IN LOCAL & GROBAL COORDINATES)'
!  WRITE(8,8009) NE,IPOIN,(STRSF(NE,IEVAB),IEVAB=IEVAB1,IEVAB2),(STRF1(NE,IEVAB),IEVAB=IEVAB1,IEVAB2)
  8009 FORMAT(2I3,6E12.4,3X,6E12.4)
          ENDDO
!
  800   CONTINUE
!
! ---------------------------------------------------------------------
      ELSEIF(INDEX.EQ.3)THEN ! *** ＜曲げ要素の部材応力をファイル＃８に出力する＞ ***
! ---------------------------------------------------------------------
        WRITE(8,8010)
 8010   FORMAT(//,' (MEMBER FORCES OF FRAME ELEMENTS IN LOCAL COORDINATES)',/,' ELEM',&
               12X,'PX1',12X,'PY1',12X,'PZ1',11X,'Myz1',11X,'Mxz1',11X,'Mxy1',&
               17X,'PX2',12X,'PY2',12X,'PZ2',11X,'Myz2',11X,'Mxz2',11X,'Mxy2')
        DO NE=1,NELEF
          WRITE(8,40) NE,(STRSF(NE,IEVAB),IEVAB=1,NDOFN*2)
        ENDDO
   40   FORMAT(I5,2(6F15.5,5X))
!
      ENDIF
goto 1000
  999 call echo
 1000 RETURN
      END
!
!
!
! ======================================================================  
        SUBROUTINE STIFF(SK,EA,EIz,EIy,GIp,GA,CAL,VALKAPPA,VALK1,VALK2)
!=======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SK(6,6)
      DIMENSION H(6,6)  ! マトリクスＨ
      DIMENSION F(6,6)  ! マトリクスｆ
!
      CALL ZEROR2(SK,6,6,6,6)
!
      VALKAPPAy=VALKAPPA
      VALKAPPAz=VALKAPPA
      GAy=GA
      GAz=GA
      GAMMAy=EIz*VALKAPPAy/GAy/CAL**2
      GAMMAz=EIy*VALKAPPAz/GAz/CAL**2
      VALK12=VALK1*VALK2
      VALKy=4.0*(VALK12+VALK1+VALK2)+3.0+12.0*GAMMAy*(4.0*VALK12+VALK1+VALK2)
      VALKz=4.0*(VALK12+VALK1+VALK2)+3.0+12.0*GAMMAz*(4.0*VALK12+VALK1+VALK2)
      SK(1,1)=EA/CAL
      SK(4,4)=GIp/CAL
      SK(2,2)=12.0*EIz/CAL**3*(4.0*VALK12+VALK1+VALK2)/VALKy
      SK(2,6)=-12.0*EIz/CAL**2*VALK2*(2.0*VALK1+1.0)/VALKy
      SK(3,3)=12.0*EIy/CAL**3*(4.0*VALK12+VALK1+VALK2)/VALKz
      SK(3,5)=12.0*EIy/CAL**2*VALK2*(2.0*VALK1+1.0)/VALKz
      SK(5,5)=4.0*EIy/CAL*VALK2*(4.0*VALK1+12.0*GAMMAz*VALK1+3.0)/VALKz
      SK(6,6)=4.0*EIz/CAL*VALK2*(4.0*VALK1+12.0*GAMMAy*VALK1+3.0)/VALKy
      SK(6,2)=SK(2,6)
      SK(5,3)=SK(3,5)
!
      RETURN
      END   
!
!
!
!         
!======================================================================
      SUBROUTINE GEOMAT(NELEM,NDOFN,NE,SKG,CAL1,STRSM,EIy,EIz,EA)
!======================================================================
! *** 曲げ要素の幾何剛性マトリクスを作成するサブルーチン。
!     各成分は「座屈問題解析」（川井忠彦著・培風館）PP.91～96を参考にしている。
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SKG(12,12)                            ! SKG(12,12)：幾何剛性マトリクス
      DIMENSION SKGF(4,4)                             ! SKMNF(4,4)：ｋGf→(Ａ-1)t・ｋGf→(Ａ-1)t・ｋGf・Ａ-1
      DIMENSION AINV(4,4)                             ! AINV(4,4) ：Ａ-1
      DIMENSION FK(0:6)                               ! FK(6)     ：ｆk＝∫f(z)z^kdz
      DIMENSION CL(0:6)                               ! CL(6)     ：CL(n)=CAL1**n
      DIMENSION STRSM(NELEM,NDOFN*2)
      INTEGER NORDR(12)                               ! NORDR(12) ：SKGF(4,4)をSKG(12,12) に代入する時の代入先
      DATA NORDR /2,6,8,12,3,5,9,11,4,0,10,0/
!
!     NDOFN ：１節点当りの自由度数(6)
!     CAL1  ：部材の長さ
!     FZ    ：部材端荷重（Px,Py,Pz,Myz,Mzx,Mxy）
!      PI=ATAN(1.0)*4.0
!
     CALL ZEROR2(SKG,12,12,12,12)
! 
      DO 1000 IDOFN=1,NDOFN
        FZ=STRSM(NE,IDOFN)
!
        DO K=0,6
          FK(K)=CAL1**(K+1)/FLOAT(K+1)
          CL(K)=CAL1**K
          FK(K)=FK(K)*FZ
        END DO
!
!  WRITE(8,8009) NE,FZ,(STRSM(NE,IEVAB),IEVAB=1,12)
!  8009 FORMAT(I3,3X,E12.4,6X,6E12.4,3X,6E12.4)
!  WRITE(8,8010) (FK(K),K=0,6)
!  8010 FORMAT(24X,7E12.4)
!
          IF(IDOFN.EQ.5.OR.IDOFN.EQ.6)THEN
          DO K=0,6
              WOO1=STRSM(NE,IDOFN)
              WOO2=STRSM(NE,IDOFN+NDOFN)
                FK(K)=CAL1**(K+1)/FLOAT(K+2)*(WOO1/FLOAT(K+1)+WOO2)
          END DO
!  WRITE(8,8010) (FK(K),K=0,6)
           ENDIF  
!
!
        CALL ZEROR2(SKGF,4,4,4,4)
        IF(IDOFN.EQ.1.OR.IDOFN.EQ.5.OR.IDOFN.EQ.6)THEN
          SKGF(1,1)=36.0*FK(2)/CL(4)-72.0*FK(3)/CL(5)+36.0*FK(4)/CL(6)
          SKGF(1,2)=-6.0*FK(1)/CL(2)+30.0*FK(2)/CL(3)-42.0*FK(3)/CL(4)+18.0*FK(4)/CL(5)
          SKGF(1,3)=-SKGF(1,1)
          SKGF(1,4)=12.0*FK(2)/CL(3)-30.0*FK(3)/CL(4)+18.0*FK(4)/CL(5)
          SKGF(2,1)=SKGF(1,2)
          SKGF(2,2)=FK(0)-8.0*FK(1)/CL(1)+22.0*FK(2)/CL(2)-24.0*FK(3)/CL(3)+9.0*FK(4)/CL(4)
          SKGF(2,3)=-SKGF(2,1)
          SKGF(2,4)=-2.0*FK(1)/CL(1)+11.0*FK(2)/CL(2)-18.0*FK(3)/CL(3)+9.0*FK(4)/CL(4)
          SKGF(3,1)=SKGF(1,3)
          SKGF(3,2)=SKGF(2,3)
          SKGF(3,3)= SKGF(1,1)
          SKGF(3,4)=-SKGF(1,4)
          SKGF(4,1)=SKGF(1,4)
          SKGF(4,2)=SKGF(2,4)
          SKGF(4,3)=SKGF(3,4)
          SKGF(4,4)= 4.0*FK(2)/CL(2)-12.0*FK(3)/CL(3)+9.0*FK(4)/CL(4)
        ELSEIF(IDOFN.EQ.2.OR.IDOFN.EQ.3)THEN
          SKGF(1,1)=-6.0*FK(1)/CL(2)+6.0*FK(2)/CL(3)+18.0*FK(3)/CL(4)-30.0*FK(4)/CL(5)+12.0*FK(5)/CL(6)
          SKGF(1,2)=FK(0)-4.0*FK(1)/CL(1)+14.0*FK(3)/CL(3)-17.0*FK(4)/CL(4)+6.0*FK(5)/CL(5)
          SKGF(1,3)=-SKGF(1,1)
          SKGF(1,4)=-2.0*FK(1)/CL(1)+3.0*FK(2)/CL(2)+6.0*FK(3)/CL(3)-13.0*FK(4)/CL(4)+6.0*FK(5)/CL(5)
          SKGF(2,1)=-6.0*FK(2)/CL(2)+18.0*FK(3)/CL(3)-18.0*FK(4)/CL(4)+6.0*FK(5)/CL(5)
          SKGF(2,2)=FK(1)-6.0*FK(2)/CL(1)+12.0*FK(3)/CL(2)-10.0*FK(4)/CL(3)+3.0*FK(5)/CL(4)
          SKGF(2,3)=-SKGF(2,1)
          SKGF(2,4)=-2.0*FK(2)/CL(1)+7.0*FK(3)/CL(2)-8.0*FK(4)/CL(3)+3.0*FK(5)/CL(4)
          SKGF(3,1)=-18.0*FK(3)/CL(4)+30.0*FK(4)/CL(5)-12.0*FK(5)/CL(6)
          SKGF(3,2)=3.0*FK(2)/CL(2)-14.0*FK(3)/CL(3)+17.0*FK(4)/CL(4)-6.0*FK(5)/CL(5)
          SKGF(3,3)=-SKGF(3,1)
          SKGF(3,4)=-6.0*FK(3)/CL(3)+13.0*FK(4)/CL(4)-6.0*FK(5)/CL(5)
          SKGF(4,1)=6.0*FK(3)/CL(3)-12.0*FK(4)/CL(4)+6.0*FK(5)/CL(5)
          SKGF(4,2)=-1.0*FK(2)/CL(1)+5.0*FK(3)/CL(2)-7.0*FK(4)/CL(3)+3.0*FK(5)/CL(4)
          SKGF(4,3)=-SKGF(4,1)
          SKGF(4,4)=2.0*FK(3)/CL(2)-5.0*FK(4)/CL(3)+3.0*FK(5)/CL(4)
        ELSEIF(IDOFN.EQ.4)THEN
          SKGF(1,1)=36.0*FK(1)/CL(4)-108.0*FK(2)/CL(5)+72.0*FK(3)/CL(6)
          SKGF(1,2)=24.0*FK(1)/CL(3)-60.0*FK(2)/CL(4)+36.0*FK(3)/CL(5)
          SKGF(1,3)=-SKGF(1,1)
          SKGF(1,4)=12.0*FK(1)/CL(3)-48.0*FK(2)/CL(4)+36.0*FK(3)/CL(5)
          SKGF(2,1)=-6.0*FK(0)/CL(2)+36.0*FK(1)/CL(3)-66.0*FK(2)/CL(4)+36.0*FK(3)/CL(5)
          SKGF(2,2)=-4.0*FK(0)/CL(1)+22.0*FK(1)/CL(2)-36.0*FK(2)/CL(3)+18.0*FK(3)/CL(4)
          SKGF(2,3)=-SKGF(2,1)
          SKGF(2,4)=-2.0*FK(0)/CL(1)+14.0*FK(1)/CL(2)-30.0*FK(2)/CL(3)+18.0*FK(3)/CL(4)
          SKGF(3,1)=SKGF(1,3)
          SKGF(3,2)=-24.0*FK(1)/CL(3)+60.0*FK(2)/CL(4)-36.0*FK(3)/CL(5)
          SKGF(3,3)=-SKGF(3,1)
          SKGF(3,4)=-12.0*FK(1)/CL(3)+48.0*FK(2)/CL(4)-36.0*FK(3)/CL(5)
          SKGF(4,1)=12.0*FK(1)/CL(3)-42.0*FK(2)/CL(4)+36.0*FK(3)/CL(5)
          SKGF(4,2)= 8.0*FK(1)/CL(2)-24.0*FK(2)/CL(3)+18.0*FK(3)/CL(4)
          SKGF(4,3)=-SKGF(4,1)
          SKGF(4,4)= 4.0*FK(1)/CL(2)-18.0*FK(2)/CL(3)+18.0*FK(3)/CL(4)
        ENDIF
!
        DO IL=1,4
          DO IR=1,4
            IF(IDOFN.EQ.1)THEN
              JL1=NORDR(IL)
              JR1=NORDR(IR)
              JL2=NORDR(IL+4)
              JR2=NORDR(IR+4)
                    JL3=NORDR(IL+8)
              JR3=NORDR(IR+8)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)+SKGF(IL,IR)
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)+SKGF(IL,IR)
              IF(JL3*JR3.NE.0)SKG(JL3,JR3)=SKG(JL3,JR3)+(EIy+EIz)/EA*SKGF(IL,IR)
                  ELSEIF(IDOFN.EQ.2)THEN
              JL1=NORDR(IL+4)
              JR1=NORDR(IR+8)
              JL2=NORDR(IL+8)
              JR2=NORDR(IR+4)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)+SKGF(IR,IL)
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)+SKGF(IL,IR)
            ELSEIF(IDOFN.EQ.3)THEN
              JL1=NORDR(IL)
              JR1=NORDR(IR+8)
              JL2=NORDR(IL+8)
              JR2=NORDR(IR)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)-SKGF(IR,IL)
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)-SKGF(IL,IR)
            ELSEIF(IDOFN.EQ.4)THEN
              JL1=NORDR(IL)
              JR1=NORDR(IR+4)
              JL2=NORDR(IL+4)
              JR2=NORDR(IR)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)-0.5*(SKGF(IL,IR)-SKGF(IR,IL))
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)+0.5*(SKGF(IL,IR)-SKGF(IR,IL))
            ELSEIF(IDOFN.EQ.5)THEN
              JL1=NORDR(IL)
              JR1=NORDR(IR+8)
              JL2=NORDR(IL+8)
              JR2=NORDR(IR)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)-SKGF(IL,IR)
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)-SKGF(IL,IR)
            ELSEIF(IDOFN.EQ.6)THEN
              JL1=NORDR(IL+4)
              JR1=NORDR(IR+8)
              JL2=NORDR(IL+8)
              JR2=NORDR(IR+4)
              IF(JL1*JR1.NE.0)SKG(JL1,JR1)=SKG(JL1,JR1)-SKGF(IL,IR)
              IF(JL2*JR2.NE.0)SKG(JL2,JR2)=SKG(JL2,JR2)-SKGF(IL,IR)
            ENDIF
          END DO
        END DO
 1000 END DO
!
        DO JR=1,12
               SKG(5,JR)=-SKG(5,JR)
               SKG(11,JR)=-SKG(11,JR)
            END DO
          DO JL=1,12
               SKG(JL,5)=-SKG(JL,5)
               SKG(JL,11)=-SKG(JL,11)
            END DO
!
      RETURN 
      END
!
!
!
! ======================================================================
      SUBROUTINE CNECT1(C0,CT,CAL1,CAL2,DX,DY,DZ,BETA,TRT)
!=======================================================================
! *** Ｃ0，ＣＴマトリクスの作成
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION C0(12,6)
      DIMENSION CT(6,12)
      DIMENSION TRT(12,12)    ! TRT(12,12)：座標変換マトリクスＴＲの転置
!
      CALL ZEROR2(C0,12,6,12,6)
      CALL ZEROR2(CT,6,12,6,12)
!
      DO I=1,6
        C0(I,I)=-1
        C0(I+6,I)=1
      ENDDO
      C0(5,3)=CAL1
      C0(6,2)=-CAL1
!     CALL WMTR2(2,C0,12,6,12,6,'C0   ','CNCT1',1)
!
      DO I=1,3
        DO J=1,3
          CT(I,J+6)=TRT(I,J)
        ENDDO
      ENDDO
!
      DO I=1,3
        DO J=7,9
          CT(I+3,J+3)=CT(I,J)
        ENDDO
      ENDDO
!
      DO J=1,6
        DO I=1,6
          CT(I,J)=-CT(I,J+6)
        ENDDO
      ENDDO
!
      DO J=1,3
        CT(2,J+3)=-CAL1*TRT(3,J)
        CT(3,J+3)=CAL1*TRT(2,J)
      ENDDO
!     CALL WMTR2(2,CT,6,12,6,12,'CT   ','CNCT1',1)
!
      RETURN
      END
!
!
!
!======================================================================
      SUBROUTINE TRANST(TRT,CAL1,CAL2,DX,DY,DZ,BETA,VECLY)
!======================================================================
! *** 座標変換マトリクスの作成
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TRT(12,12)    ! TRT(12,12)：座標変換マトリクスＴの転置
      DIMENSION VECLX(3)      ! VECLX(3)：局所座標系のｘ軸方向を表わすベクトル
      DIMENSION VECLY(3)      ! VECLY(3)：局所座標系のｙ軸方向を表わすベクトル
      DIMENSION VECLZ(3)      ! VECLZ(3)：局所座標系のｚ軸方向を表わすベクトル
!
      CALL ZEROR2(TRT,12,12,12,12)
!
      DLENG=SQRT(DX**2+DY**2+DZ**2)
      VECLX(1)=DX/DLENG
      VECLX(2)=DY/DLENG
      VECLX(3)=DZ/DLENG
!      IF(VECLX(1).EQ.0.0.AND.VECLX(2).EQ.0.0)THEN     !      IF(VECLX(1).EQ.0.0.AND.VECLX(3).EQ.0.0)THEN    !
!        VECLY(1)=0.0                                  !        VECLZ(1)=0.0                                 !
!        VECLY(2)=1.0                                  !        VECLZ(2)=0.0                                 !
!        VECLY(3)=0.0                                  !        VECLZ(3)=1.0                                 !
!      ELSE                                            !      ELSE                                           !
!        DLENG=SQRT(VECLX(1)**2+VECLX(2)**2)           !        DLENG=SQRT(VECLX(1)**2+VECLX(3)**2)          !
!        VECLY(1)=-VECLX(2)/DLENG                      !        VECLZ(1)=-VECLX(3)/DLENG                     !
!        VECLY(2)=VECLX(1)/DLENG                       !        VECLZ(2)=0.0                                 !
!        VECLY(3)=0.0                                  !        VECLZ(3)=VECLX(1)/DLENG                      !
!      ENDIF                                           !      ENDIF                                          !
!        VECLY(1)=-1.0/SQRT(2.0)
!        VECLY(2)= 1.0/SQRT(2.0)
!        VECLY(3)=0.0
      CALL VECTPRD(VECLX,VECLY,VECLZ)                 !      CALL VECTPRD(VECLZ,VECLX,VECLY)                !
      TRT(1,1)=VECLX(1)
      TRT(1,2)=VECLX(2)
      TRT(1,3)=VECLX(3)
      TRT(2,1)=VECLY(1)
      TRT(2,2)=VECLY(2)
      TRT(2,3)=VECLY(3)
      TRT(3,1)=VECLZ(1)
      TRT(3,2)=VECLZ(2)
      TRT(3,3)=VECLZ(3)
!       CALL WMTR2(2,TRT,12,12,6,6,'TRT  ','TRNST',2)
!
      DO N=3,9,3
        DO I=1,3
          DO J=1,3
            TRT(I+N,J+N)=TRT(I,J)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN 
      END
!
!
!
!======================================================================
      SUBROUTINE VECTPRD(A,B,C)
!======================================================================
! *** ベクトルＡ、Ｂの外積Ｃを求める。
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3),B(3),C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
