! =====================================================================
SUBROUTINE INPUT(ANG,COORD,DDV,DISMAX,ERR,FCOE,IFPRE,LODFIX, &
                 MAXCYL,MEIGN,MSTRE,NANG,NANGLE,NCONC,NDOFN,NELEM,NFPOIN,&
                 NLOAD,NODEM,NOFIX,NONISO,NPOIN,NPRET,NPRINT,NPROB,NREPT,  &
                 NSTEP,NTIMES,NVFIX,NVARM,PMGNF,PPRES,PRESC,   &
                 PRETEN,PROPM,SCF,SSNOW,STRSM,WCF,WIND,WWIND,MINCR)
! =====================================================================
      PARAMETER( MFORMAT=1 )    ! MFORMAT：固定フォーマット（=0）、配列データの一部はフリーフォーマット（=1）
! 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ANG(NELEM)
      DIMENSION COORD(NPOIN,3)
      DIMENSION FCOE(NPOIN*NDOFN)
      DIMENSION IFPRE(NVFIX)
      DIMENSION NANG(NELEM)
      DIMENSION NODEM(NELEM,4)          ! 膜要素のデータ（1～3：節点番号、4：部材番号）
      DIMENSION NOFIX(NVFIX)
      DIMENSION PMGNF(NDOFN,NTIMES)
      DIMENSION PPRES(NTIMES)
      DIMENSION PRESC(NVFIX,NDOFN)
      DIMENSION SCF(NELEM)
      DIMENSION SSNOW(NTIMES)
      DIMENSION STRSM(NELEM,MSTRE)
      DIMENSION WWIND(NTIMES)
      DIMENSION WCF(NELEM)

DIMENSION PROPM(NVARM,7,3)        !膜要素の材料特性（種類番号、特性項目(1～7)、段階Tri-Linear(1～3)）
character*150 text
!
! *** ＜荷重倍率、内圧、速度圧、雪荷重の設定＞ ***
!
do i=1,7
  read(7,'(A)')text
  write(8,'(A)')text
enddo
do i=1,NTIMES
  read(7,*,err=999) NUM,(PMGNF(j,i),j=1,NDOFN),PPres(i),WWind(i),SSnow(i)
enddo
write(8,117)(i,(PMGNF(j,i),j=1,NDOFN),PPres(i),WWind(i),SSnow(i),i=1,NTIMES)
117 format(I15,9G15.4)

! *** ＜節点座標、直接節点荷重の設定＞ ***

do i=1,6
  read(7,'(A)')text
  write(8,'(A)')text
enddo
if(NCONC.eq.0)then
  do i=1,NPOIN
    read(7,*,err=999) k,(COORD(i,j),j=1,3)
  enddo
  write(8,130)(i,(COORD(i,j),j=1,3),i=1,NPOIN)
else
  do i=1,NPOIN
    read(7,*,err=999) k,(COORD(i,j),j=1,3),(FCOE((i-1)*NDOFN+j),j=1,NDOFN)
  enddo
  write(8,131)(i,(COORD(i,j),j=1,3),(FCOE((i-1)*NDOFN+j),j=1,NDOFN),i=1,NPOIN)
endif
130 format(I15,3F15.4)
131 format(I15,9F15.4)

! *** ＜拘束条件の設定＞ ***

call ZEROR2(PRESC,NVFIX,NDOFN,NVFIX,NDOFN)
do i=1,7
  read(7,'(A)')text
  write(8,'(A)')text
enddo
do i=1,NVFIX
  read(7,*,err=999) NOFIX(i),IFPRE(i),(PRESC(i,j),j=1,NDOFN)
enddo
write(8,908)(NOFIX(i),IFPRE(i),(PRESC(i,j),j=1,NDOFN),i=1,NVFIX)
908 format(2I15,6F15.4)

! *** ＜三角形膜要素の部材特性の設定＞ ***

if(NVARM.NE.0)then
  do i=1,6
    read(7,'(A)')text
    write(8,'(A)')text
  enddo
  do i = 1,NVARM
!PROPM(NVARM,7,3)   !膜要素の材料特性（種類番号IVARM、特性項目(1～7)、段階Tri-Linear(1～3)）
    read(7,*,err=999)text,(PROPM(i,j,1),j=1,7)
    do k=2,3
      read(7,*,err=999)text,(PROPM(i,j,k),j=1,5)
    enddo
    write(8,168)i,(PROPM(i,j,1),j=1,7),((PROPM(i,j,k),j=1,5),k=2,3)
  enddo
endif
168 format(&
I15,'     第１勾配  ',7F15.4,/&
15X,'     第２勾配  ',5F15.4,/&
15X,'     第３勾配  ',5F15.4)

! *** ＜三角形膜要素の節点番号、初期張力、荷重係数の設定＞ ***

!NSURF=0   !NSURF=0：風雪荷重なし、=1:風、=2：雪、=3:エラー
do i=1,NTIMES
  if(WWIND(i).eq.0.0.and.SSNOW(i).eq.0.0)then
    NSURF=0
  elseif(WWIND(i).ne.0.0.and.SSNOW(i).eq.0.0)then
    NSURF=1
  elseif(WWIND(i).eq.0.0.and.SSNOW(i).ne.0.0)then
    NSURF=2
  else
    NSURF=3
    write(8,901)
    901 format(///,'error　風荷重と雪荷重を同時に計算することはできません',/&
                  '       WWINDとSSNOWのどちらか一方は全てゼロにしてください')
    stop
  endif
  if(NSURF.NE.0)exit
enddo

if(NELEM.NE.0)then
  do i=1,6
    read(7,'(A)')text
    write(8,'(A)')text
  enddo
  !read(7,'(A)')(text,i=1,4)
  !write(8,135)
  if(NSURF.EQ.0)then
    do i=1,NELEM
      read(7,*,err=999) num,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3)
    enddo
    write(8,136)(i,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3),i=1,NELEM)
  elseif(NSURF.EQ.1)then
    do i=1,NELEM
      read(7,*,err=999) num,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3),WCF(i)
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@@@@ GRANPA DOME の風力係数の計算 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(NSURF.eq.1)then
        NWCOE=0  ! =1 以下のプログラムで風力係数を設定する =0 何もしない（入力データの風力係数を使って計算する）
        if(NWCOE.eq.1)then
          center=0.0
          span=13.6*2.0
          do j=1,3
            center=center+COORD(NODEM(i,j),1)/3.0
          enddo
          if(center.LE.-span/8.0)then  !風上 端部～1/8
            WCF(i)=-0.25
          elseif(center.LE.-span/4.0)then  !風上1/8～1/4
            WCF(i)=-0.25
          elseif(center.LE.0.0)then   !風上1/4～中央
            WCF(i)=-0.61
          elseif(center.LE.span/4.0)then  !中央～風下1/4
            WCF(i)=-0.61
          elseif(center.LE.span/8.0)then  !風下1/4～1/8
            WCF(i)=-0.01
          else  !風下 1/8～端部
          WCF(i)=+0.36
          endif
        endif
      endif
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    enddo
    write(8,137)(i,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3),WCF(i),i=1,NELEM)
  elseif(NSURF.EQ.2)then
    do i=1,NELEM
      read(7,*,err=999) num,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3),SCF(i)
    enddo
    write(8,137)(i,(NODEM(i,j),j=1,4),(STRSM(i,j),j=1,3),SCF(i),i=1,NELEM)
  endif
endif

135 format(//&
'膜要素                                                                            初期張力                                     荷重係数',/&
'          番号               I              J              K       種類番号     IJ方向張力 IJ直交方向張力       せん断力     （風or雪）',/&
'                                                                                     [N/m]          [N/m]          [N/m]',/&
'     I=1～NELEM     NODEM(I,1)     NODEM(I,2)     NODEM(I,3)                    STRSM(I,1)     STRSM(I,2)     STRSM(I,3)')
136 format(5I15,3F15.4)
137 format(5I15,4F15.4)

!!
!! *** ＜膜要素のIJ方向と材料の主軸方向のなす角度の設定＞ ***
!!
!            IF(NONISO.EQ.2)THEN
!              WRITE(8,182)
!    IF(MFORMAT.EQ.0)THEN
!      !READ(7,197)(NANG(IANGLE),ANG(IANGLE),IANGLE=1,NANGLE)
!    ELSE
!      !READ(7,*)(NANG(IANGLE),ANG(IANGLE),IANGLE=1,NANGLE)
!    ENDIF
!              WRITE(8,183)(NANG(IANGLE),ANG(IANGLE),IANGLE=1,NANGLE)
!            ENDIF
!  183       FORMAT(5(I5,F8.2,7X))
!  182       FORMAT(//,' (ORTHOTROPIC PROPERTY --ANGLE BETWEEN THE LOCAL COORDINATE AND THE DIRECTION OF E1 IN DEGREE)')
!!
goto 1000
999 call echo
1000 RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE CABLE1(BL,CAL,COORD,DV,NODEC,NPOIN,NELEC,&
                        NPRSTR,NVARC,STRSC,STRC1,PROPC,&
                        MCTEMP,NTIMES,CTEMP1,CTEMP2)
! =====================================================================
! *** ＜ケーブルまたはトラスに関するデータをファイル＃７から読込む＞ ***
!     NELEC：ケーブルまたはトラスの要素数
!     NVARC：部材数
!     NPRSTR：ケーブルまたはトラスの初期応力（=0：なし、=1：あり）
!     DV：ケーブルが緩んだ場合に、剛性をEA/DVとする。通常は1.0とする。
!
      PARAMETER( MFORMAT=1 )    ! MFORMAT：固定フォーマット（=0）、配列データの一部はフリーフォーマット（=1）
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION NODEC(NELEC,3)   !ケーブル要素の節点番号(1～2)、部材番号(3)
      DIMENSION BL(NELEC)        !   BL(IELEC)：現在の長さ
      DIMENSION CAL(NELEC)       !   CAL(IELEC)：無ひずみ時の長さ
      DIMENSION STRSC(NELEC)     !   STRSC( )：初期応力
      DIMENSION STRC1(NELEC)
      dimension PROPC(NVARC,3)   !ケーブル要素の材料データ（ヤング率、断面積、単位重量）
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
      DIMENSION CTEMP1(NTIMES)   ! CTEMP1(LOADCREM)：ケーブルに与える温度増分[°]、MCTEMP=1（ケーブルの温度応力解析を行う）の場合
      DIMENSION CTEMP2(2,NELEC)  ! CTEMP2(1,IELEC)：温度応力を与える(=1)、与えない(=0)、CTEMP2(2,IELEC)：温度ひずみ
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
character*150 text

! *** ＜ケーブルまたはトラスの部材番号、断面積、密度、ヤング率＞

do i=1,6
  read(7,'(A)')text
  write(8,'(A)')text
enddo
read(7,*,err=999) (num,PROPC(i,1),PROPC(i,2),PROPC(i,3),i=1,NVARC)
write(8,101) (i,PROPC(i,1),PROPC(i,2),PROPC(i,3),i=1,NVARC)
101 format(I15,2F15.4,G15.4)

! *** ＜ケーブルまたはトラス要素の節点番号、部材番号、初期張力＞

read(7,'(A)')(text,i=1,6)
read(7,*,err=999)(num,(NODEC(i,j),j=1,3),STRSC(i),i=1,NELEC)

! *** ＜ケーブルまたはトラスの要素諸元（長さ、重量、剛性）の計算と出力＞

call ZEROR1(STRC1,NELEC,NELEC)
WSUM=0.0
write(8,102)
DO NE=1,NELEC
  STRC1(NE)=STRSC(NE)
  I=NODEC(NE,1)
  J=NODEC(NE,2)
  IV=NODEC(NE,3)
  BL(NE)=SQRT((COORD(I,1)-COORD(J,1))**2+(COORD(I,2)-COORD(J,2))**2+(COORD(I,3)-COORD(J,3))**2)
  EA=PROPC(IV,1)*DABS(PROPC(IV,2))   !EA[N]=ヤング率[N/mm2]*断面積[mm2]
  IF((STRSC(NE).EQ.0.0).OR.(EA.EQ.0.0))THEN
    CAL(NE)=BL(NE)
  ELSE
    CAL(NE)=BL(NE)/(STRSC(NE)/EA+1.0)
  ENDIF
!!!!!単位系変更  W=PROPC(IV,3)*CAL(NE)   !重量[kgf]＝単位重量[kg/m]*長さ[m]
  W=PROPC(IV,3)*CAL(NE)*9.8   !重量[N]＝単位重量[kg/mm2/m]*長さ[m]
  WSUM=WSUM+W
  write(8,103)NE,(NODEC(NE,i),i=1,3),STRSC(NE),PROPC(IV,2),BL(NE),W
ENDDO
102 format(//&
'ケーブル要素',/&
'       要素番号       節点番号                      種類番号       初期張力         断面積       部材長さ           重量',/&
'                           I端            J端                           [N]          [cm2]            [m]            [N]',/&
'      I=1～NELEC     NODEC(I,1)     NODEC(I,2)     NODEC(I,3)         STRSC(I)     PROPC(I,2)                              ')
103 format(4I15,4F15.4)
WRITE(8,337) WSUM
337 FORMAT(4X,120('-'),/,85X,'  ケーブルの総重量 (',F15.4,' [N]')
!
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
      IF(MCTEMP.EQ.1)THEN
        READ(7,*,err=999) (I,CTEMP1(ITIMES),ITIMES=1,NTIMES)    ! CTEMP1(LOADCREM)：ケーブルに与える温度増分[°]、MCTEMP=1（ケーブルの温度応力解析を行う）の場合
        WRITE(8,410) (ITIMES,CTEMP1(ITIMES),ITIMES=1,NTIMES)
  410   FORMAT(/,5X,' TEMPERATURE INCREMENT OF CABLE ELEMENTS',/,'  LOADCREM TEMPERATURE [DEGREE]',/,(I10,F10.4))
        READ(7,*,err=999) NCTEMP    ! NCTEMP：ケーブルの温度応力計算の対象とする要素の数
        IF(NCTEMP.LT.NELEC)THEN
          READ(7,*,err=999) (IELEC,CTEMP2(1,IELEC),ICTEMP=1,NCTEMP)   ! CTEMP2(1,IELEC)：温度応力を与える(=1)、与えない(=0)、CTEMP2(2,IELEC)：温度ひずみ
        ELSE
          DO IELEC=1,NELEC    ! NCTEMP=NELECとした場合には、CTEMP2(1,IELEC)を入力しなくても自動的に1.0に設定する。
            CTEMP2(1,IELEC)=1.0
          ENDDO
        ENDIF
        WRITE(8,420)
  420   FORMAT(/,5X,' CABLE ELEMENT NUMBER THAT TEMPERATURE STRAIN IS CONSIDERED')
        DO IELEC=1,NELEC
          IF(CTEMP2(1,IELEC).NE.0.0)WRITE(8,421)IELEC,CTEMP2(1,IELEC)
  421     FORMAT(I10,F10.4)
        ENDDO
      ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
goto 1000
  999 call echo
1000  RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE PROC(NFLAG,DIS,EIGMN,IFFIX,INCREM,LOADCREM,  &
                      NPOIN,NCYCLE,NDOFN,NFPOIN,NNN,NPRINT,NPROB,   &
                      NTIMES,NTOTV,PRES,Q,COORD,QMAX,MINCR)
! =====================================================================
! *** 計算過程を画面、ファイルに出力する
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NP(6),NQ(6),QQ(6),DD(6)
      DIMENSION DIS(NPOIN*NDOFN)
      DIMENSION IFFIX(NPOIN*NDOFN)
      DIMENSION Q(NPOIN*NDOFN)
      DIMENSION COORD(NPOIN,3)
!
      CALL ZEROR1(QQ,6,6)
      CALL ZEROR1(DD,6,6)
!
! *** 最大変位DDと最大残差力QQの計算
      DO IPOIN=1,NFPOIN
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IF(IFFIX(ITOTV).EQ.0)THEN
            IF(ABS(DIS(ITOTV)).GT.ABS(DD(IDOFN)))THEN
              DD(IDOFN)=DIS(ITOTV)
              NP(IDOFN)=IPOIN
            ENDIF
            IF(ABS(Q(ITOTV)).GT.ABS(QQ(IDOFN)))THEN
              QQ(IDOFN)=Q(ITOTV)
              NQ(IDOFN)=IPOIN
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
      IF(NFLAG.EQ.0)THEN
        IF(NPRINT.EQ.0)THEN
          WRITE(8,570) INCREM,NCYCLE,NNN
          WRITE(8,580) ( K,(DIS((K-1)*NDOFN+IDOFN),IDOFN=1,NDOFN),K,(Q((K-1)*NDOFN+IDOFN),IDOFN=1,NDOFN),K=1,NFPOIN)
        END IF
        IF(NCYCLE.EQ.1)WRITE(8,180)
        WRITE(8,170) NCYCLE,(NP(I),DD(I),I=1,6),(NQ(I),QQ(I),I=1,6)
      ELSE IF(NFLAG.EQ.1)THEN
        Qmax = QQ(1)
        IF(ABS(QQ(2)).GE.ABS(Qmax))Qmax=QQ(2)
        IF(ABS(QQ(3)).GE.ABS(Qmax))Qmax=QQ(3)
        WRITE(8,710)
        IF(NNN.EQ.1)THEN
          WRITE(6,722)
          IF(MINCR.EQ.1)WRITE(11,722)
        ENDIF
        NPROB1=ABS(NPROB)
        WRITE(6,723)  loadcrem, ntimes, nnn, increm, ncycle, nprob1,(COORD(NPROB1,IDOFN),IDOFN=1,3),PRES,Qmax,EIGMN
        IF(NPRINT.LE.2)WRITE(8,723) loadcrem, ntimes, nnn, increm, ncycle, nprob1,(COORD(NPROB1,IDOFN),IDOFN=1,3),PRES,Qmax,EIGMN
        IF(MINCR.EQ.1)WRITE(11,723) loadcrem, ntimes, nnn, increm, ncycle, nprob1,(COORD(NPROB1,IDOFN),IDOFN=1,3),PRES,Qmax,EIGMN
      ENDIF
!
  570 FORMAT(//,' INCREM=',I3,5X,'NCYCLE=',I3,5X,' NNN=',I3)
  580 FORMAT(' NODE',5X,'X-DISP.',8X,'Y-DISP.',8X,'Z-DISP.',7X,'XX-DISP.',7X,'YY-DISP.',7X,'ZZ-DISP.',  &
                    13X,'X-UNBAL',8X,'Y-UNBAL',8X,'Z-UNBAL',7X,'XX-UNBAL',7X,'YY-UNBAL',7X,'ZZ-UNBAL',/,    &
             (I5,6E15.5,I10,6E15.5))
  170 FORMAT(I5,6(1X,'(',I5,')',E11.3),5X,6(1X,'(',I5,')',E11.3))
  180 FORMAT(/,' NEWTON       MAX.               MAX.               MAX.               MAX.               MAX.               MAX.                  MAX. X-            MAX. Y-            MAX. Z-           MAX. XX-           MAX. YY-           MAX. ZZ-',/,&
               ' ITERAT       X-DISP.            Y-DISP.            Z-DISP.           XX-DISP.           YY-DISP.           ZZ-DISP.               UNBALANCED         UNBALANCED         UNBALANCED        UNBALANCED         UNBALANCED         UNBALANCED')
  710 FORMAT(//,' (CONVERGED VALUES)')
  722 FORMAT(8X,' NNN INC CYC  No( X-COOD, Y-COOD, Z-COOD)    PRES  MAXUNBAL  MINEIGEN')
  723 FORMAT(' (',I2,'/',I2,')',3I4,I6,'(',2(F7.4,','),F7.4,') ',F7.2,1X,2E10.3)
!
      IF(Qmax.GT.1.0E+08)THEN
        WRITE(6,900)
        WRITE(8,900)
        WRITE(11,900)
        STOP
      ENDIF
  900 FORMAT(' This calculation tends to diverge.')
!
      RETURN
      END
!
!
!
! =====================================================================
SUBROUTINE RESLT(BL,CAL,COOR0,COORD,FF,FFF,IFFIX,MSTRE,NDOFN,&
                 NELEM,NELEC,NFINAL,NFPOIN,NPOIN,NVARC,NVARM,&
                 NODEC,NODEM,Q,STRSM,STRSC,TTDIS,STRC1,PROPC,PROPM,NELEF,STRSF,&
                 MCTEMP,CTEMP2,NPROB,NVARF,NODEF)
! =====================================================================
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION IFFIX(NPOIN*NDOFN)
DIMENSION TTDIS(NPOIN*NDOFN)
DIMENSION STRSM(NELEM,MSTRE)
DIMENSION STRSF(NELEF,NDOFN*2)                                    ! STRSF(IELEF,IEVAB)：曲げ要素の応力（部材座標系における部材端荷重)
DIMENSION NODEF(NELEF,3)            ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
DIMENSION FFF(NPOIN*NDOFN)
DIMENSION BL(NELEC)
DIMENSION CAL(NELEC)
DIMENSION STRSC(NELEC)
DIMENSION STRC1(NELEC)
DIMENSION NODEC(NELEC,3)   !ケーブル要素の節点番号(1～2)、部材番号(3)
dimension PROPC(NVARC,3)    !ケーブル要素の材料データ（ヤング率、断面積、単位重量）
DIMENSION COORD(NPOIN,3)
DIMENSION COOR0(NPOIN,3)
DIMENSION FF(NPOIN*NDOFN)
DIMENSION Q(NPOIN*NDOFN)

DIMENSION NODEM(NELEM,4)          ! 膜要素のデータ（1～3：節点番号、4：部材番号）
DIMENSION PROPM(NVARM,7,3)        !膜要素の材料特性（種類番号、特性項目(1～7)、段階Tri-Linear(1～3)）

! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
!     MCTEMP                        ! MCTEMP：ケーブルの温度応力解析を行わない（=0）、行う（=1）
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)：温度応力を与える(=1.0)、与えない(=0.0)、CTEMP2(2,IELEC)：温度ひずみ
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
DIMENSION FALL(6),QALL(6),UNB(6)
DIMENSION PST(3),VALMM(7,2),NVALMM(7,2)
!
! *** ＜膜要素の張力をファイル＃８に出力する＞ ***
!
IF(NELEM.GT.0)THEN
  WRITE(8,951)
  WRITE(8,975)
  951 FORMAT(///,20X'**********************************************',/,30X'RESULTS OF EXTERNAL FORCES',/,20X,'**********************************************')
  975 FORMAT(//,' (MEMBRANE STRESSES AND PRINCIPAL STRESSES)[N/m]',/,49X,'WARP-STRESS : I-J DIRECTION OF AN ELEMENT',  &
             /,49X,'PRINCIPAL ANGLE : DEGREE',/,2X,'ELEM',3X,'WARP-STRESS',4X,'FILL-STRESS',4X,'SHEAR-STRES',9X,'STRESS-1',7X,'STRESS-2',7X,'PRINCPL ANG')
  DO IELEM=1,NELEM
    PX=STRSM(IELEM,1)
    PY=STRSM(IELEM,2)
    PXY=STRSM(IELEM,3)
    PXPY=(PX+PY)*0.5
    SQ=SQRT(((PX-PY)*(PX-PY))*0.25+PXY*PXY)
    PST(1)=PXPY+SQ
    PST(2)=PXPY-SQ
    IF(PX.LT.0.0.OR.PY.LT.0.0.OR.PX.EQ.PY)THEN
      PST(3)=0.0
    ELSE
      PST(3)=28.6*ATAN((2.0*PXY)/(PX-PY))
    ENDIF
    IF(PX.LT.PY)THEN
      IF(PXY.GE.0.0)PST(3)=PST(3)+90.0
      IF(PXY.LT.0.0)PST(3)=PST(3)-90.0
    ENDIF
    WRITE(8,960)IELEM,(STRSM(IELEM,ISTRE),ISTRE=1,3),IELEM,(PST(ISTRE),ISTRE=1,3)
    960 FORMAT(I5,3F15.5,I5,6F15.5)
  ENDDO
ENDIF
!
! *** ＜膜張力の最大値VALMM(ISTRE,1)と最小値VALMM(ISTRE,2)をファイル＃８に保存する＞
!
DO IELEM=1,NELEM
  DO ISTRE=1,3
    IF(IELEM.EQ.1)THEN
      VALMM(ISTRE,1)=STRSM(IELEM,ISTRE)
      NVALMM(ISTRE,1)=IELEM
      VALMM(ISTRE,2)=STRSM(IELEM,ISTRE)
      NVALMM(ISTRE,2)=IELEM
    ELSE
      IF(STRSM(IELEM,ISTRE).GT.VALMM(ISTRE,1))THEN
        VALMM(ISTRE,1)=STRSM(IELEM,ISTRE)
        NVALMM(ISTRE,1)=IELEM
      ELSEIF(STRSM(IELEM,ISTRE).LT.VALMM(ISTRE,2))THEN
        VALMM(ISTRE,2)=STRSM(IELEM,ISTRE)
        NVALMM(ISTRE,2)=IELEM
      ENDIF
    ENDIF
  ENDDO
ENDDO
if(NELEM.gt.0)WRITE(8,816)((NVALMM(ISTRE,I),VALMM(ISTRE,I),ISTRE=1,3),I=1,2)
816 FORMAT(/,' (MAXIMUM AND MINIMUM OF MEMBRANE STRESS) [N/m]',/,'           ELEM    WARP-STRESS      ELEM    FILL-STRESS      ELEM    SHEAR-STRES',/,'  MAX',3(I10,F15.5),/,'  MIN',3(I10,F15.5))
!
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ETFEフィルムの相当応力を計算し，その最大値VALMM(1,1)と最小値VALMM(1,2)をファイル＃８に保存する
CNV=1.0E+00               ! CNV：単位の換算係数（長さがmm単位の場合1.0、長さがm単位の場合1000.0）
DO IELEM=1,NELEM
  IV=NODEM(IELEM,4)          ! NODEM(IELEM,4)：部材番号
  THICK=PROPM(IV,6,1)        ! PROPM(IV,6,1)：フィルムの厚さ[mm]
  VALUE=1.0/CNV/THICK        ! VALUE：単位を変換するための係数
  PX=STRSM(IELEM,1)*VALUE
  PY=STRSM(IELEM,2)*VALUE
  PXY=STRSM(IELEM,3)*VALUE
  F=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)    ! F：相当応力[N/mm2]
  IF(IELEM.EQ.1)THEN
    VALMM(1,1)=F
    NVALMM(1,1)=1
    VALMM(1,2)=F
    NVALMM(1,2)=1
  ELSE
    IF(F.GT.VALMM(1,1))THEN
      VALMM(1,1)=F
      NVALMM(1,1)=IELEM
    ELSEIF(F.LT.VALMM(1,2))THEN
      VALMM(1,2)=F
      NVALMM(1,2)=IELEM
    ENDIF
  ENDIF
  WRITE(10,128) IELEM,PX,PY,PXY     !相当応力をファイル＃１０に出力する
ENDDO
if(NELEM.gt.0)WRITE(8,817)(NVALMM(1,I),VALMM(1,I),I=1,2)
817   FORMAT(/,' (MAXIMUM AND MINIMUM OF EQUIVARENT MEMBRANE STRESS) [N/mm2]',/,'           ELEM         STRESS',/,'  MAX',I10,F15.5,/,'  MIN',I10,F15.5)
!
! *** ＜変形後の応力をファイル＃１０に保存する＞ ***
!
!IF(NELEM.GE.1)WRITE(10,128) (K,(STRSM(K,ISTRE),ISTRE=1,3),K=1,NELEM)
128 FORMAT(I5,3F15.4)
!
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  個別の計算結果出力
!      IP=111
!      WRITE(11,8170)IP,COORD(IP,3),NVALMM(1,1),VALMM(1,1),NVALMM(1,2),VALMM(1,2)
! 8170 FORMAT('  Z(',I4,')=',G12.5,'   ST_MAX=(',I4,',',G12.5,')','   ST_MIN=(',I4,',',G12.5,')')
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! *** ＜ケーブル要素の部材応力をファイル＃８に出力する＞ ***
!
      IF(NELEC.NE.0)CALL CABLE3(1,BL,CAL,NELEC,NVARC,NODEC,STRSC,STRC1,PROPC,MCTEMP,CTEMP2)
!
! *** ＜曲げ要素の部材応力をファイル＃８に出力する＞ ***
!
      IF(NELEF.GT.0)WRITE(8,8040) (NE,(STRSF(NE,IEVAB),IEVAB=1,NDOFN*2),NE=1,NELEF)
 8040 FORMAT(//,' (MEMBER FORCES OF FRAME ELEMENTS IN LOCAL COORDINATES)',/,' ELEM',    &
            12X,'PX1',12X,'PY1',12X,'PZ1',11X,'MYZ1',11X,'MZX1',11X,'MXY1',/,           &
            17X,'PX2',12X,'PY2',12X,'PZ2',11X,'MYZ2',11X,'MZX2',11X,'MXY2',/,           &
            (I5,6E15.4,/,5X,6E15.4))
!
! *** ＜曲げ要素の部材種類IVごとに応力の絶対値が最大となる部材の番号NVALMM(ISTRE,1)と応力の値VALMM(ISTRE,1)をファイル＃８に保存する＞
!
    IF(NELEF.GT.0)THEN
      WRITE(8,8051)
	  WRITE(11,8051)
      8051 FORMAT(/,' (MAXIMUM AND MINIMUM OF MEMBER FORCES OF FRANE ELEMENTS) ',/,&
                  17X,'PX',14X,'PY',14X,'PZ',13X,'MYZ',13X,'MZX',13X,'MXY')
	  DO IV=1,NVARF
		CALL  ZEROR2(VALMM,7,2,7,2)
		CALL  ZEROI2(NVALMM,7,2,7,2)
        DO IE=1,NELEF
          IF(NODEF(IE,3).EQ.IV)THEN
            DO ISTRE=1,6
              DO INODE=1,2
                IEVAB=(INODE-1)*6+ISTRE
                IF(ABS(STRSF(IE,IEVAB)).GT.ABS(VALMM(ISTRE,1)))THEN
                  VALMM(ISTRE,1)=STRSF(IE,IEVAB)
                  NVALMM(ISTRE,1)=IE
                ENDIF
              ENDDO
            ENDDO
	      ENDIF
        ENDDO
        WRITE(8,8052)IV,(NVALMM(ISTRE,1),VALMM(ISTRE,1),ISTRE=1,6)
        8052 FORMAT('IV=',I4,6(I4,E12.4))
!       WRITE(11,8052)IV,(NVALMM(ISTRE,1),VALMM(ISTRE,1),ISTRE=1,6)
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  個別の計算結果出力
!		  JE=NVALMM(5,1)
!          WRITE(11,8053)IV,JE,(STRSF(JE,ISTRE),ISTRE=1,6)
! 8053 FORMAT(17X,10X,'PX',10X,'PY',10X,'PZ',9X,'MYZ',9X,'MZX',9X,'MXY',/,&
!             'MAX-MYY(IV',I1,':',I4,'):',6E12.4)
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	  ENDDO
    ENDIF
!
! *** ＜変形後の座標をファイル＃９に保存する＞ ***
!
      CALL WRTCOD(2,COORD,NODEC,NELEM,NELEC,NODEM,NPOIN,NELEF,NODEF)
!
! *** ＜変形後のケーブル張力、曲げ部材力をファイル＃１０に保存する＞ ***
!
      IF(NELEC.GE.1)WRITE(10,990) (K,STRSC(K),K=1,NELEC)
  990 FORMAT(4(I5,F15.4,5X))
!
!     WRITE(10,8041) (NE,(STRSF(NE,IEVAB),IEVAB=1,NDOFN*2),NE=1,NELEF)
      WRITE(10,8042) (NE,(STRSF(NE,IEVAB),IEVAB=7,7),NE=1,NELEF)
 8041   FORMAT(I5,12E15.4)
 8042   FORMAT(4(I5,E15.4))
!
! *** ＜変形後の座標と初期形状に対する変位をファイル＃８に保存する＞ ***
!
      WRITE(8,940)
      DO IPOIN=1,NPOIN
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          Q(ITOTV)=Q(ITOTV)-FFF(ITOTV)
          FF(ITOTV)=FFF(ITOTV)
          IF(IPOIN.GT.NFPOIN)TTDIS(ITOTV)=0.0
        ENDDO
        ITOTV1=(IPOIN-1)*NDOFN+1
        ITOTV2=IPOIN*NDOFN
        WRITE(8,980) IPOIN,(COORD(IPOIN,IDOFN),IDOFN=1,3),IPOIN,(TTDIS(ITOTV),ITOTV=ITOTV1,ITOTV2)
      ENDDO
  940 FORMAT(///,' (DEFORMED CO-ORDINATES AND TOTAL DISPLACEMENTS) ',' [m]',/,  &
            '  NODE',4X,'X',14X,'Y',14X,'Z',19X,'X-DISP.',8X,'Y-DISP.',8X,'Z-DISP.',7X,'XX-DISP.',7X,'YY-DISP.',7X,'ZZ-DISP.')
  980 FORMAT(I5,3E15.5,I5,6E15.5)
!
! *** ＜節点変位の最大値VALMM(IDOFN,1)と最小値VALMM(IDOFN,2)をファイル＃８に保存する＞
!
      DO IPOIN=1,NPOIN
        VLENG=0.0   ! VLENG：変位ベクトルの長さ
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IF(IDOFN.LE.3)VLENG=VLENG+TTDIS(ITOTV)**2
          IF(IPOIN.EQ.1)THEN
            VALMM(IDOFN,1)=TTDIS(ITOTV)
            NVALMM(IDOFN,1)=IPOIN
            VALMM(IDOFN,2)=TTDIS(ITOTV)
            NVALMM(IDOFN,2)=IPOIN
          ELSE
            IF(TTDIS(ITOTV).GT.VALMM(IDOFN,1))THEN
              VALMM(IDOFN,1)=TTDIS(ITOTV)
              NVALMM(IDOFN,1)=IPOIN
            ELSEIF(TTDIS(ITOTV).LT.VALMM(IDOFN,2))THEN
              VALMM(IDOFN,2)=TTDIS(ITOTV)
              NVALMM(IDOFN,2)=IPOIN
            ENDIF
          ENDIF
        ENDDO
        VLENG=SQRT(VLENG)
        IF(IPOIN.EQ.1)THEN
          VALMM(NDOFN+1,1)=VLENG
          NVALMM(NDOFN+1,1)=1
        ELSEIF(VLENG.GT.VALMM(NDOFN+1,1))THEN
          VALMM(NDOFN+1,1)=VLENG
          NVALMM(NDOFN+1,1)=IPOIN
        ENDIF
      ENDDO
      WRITE(8,871)((NVALMM(IDOFN,I),VALMM(IDOFN,I),IDOFN=1,NDOFN),I=1,2)
      WRITE(8,872)NVALMM(NDOFN+1,1),VALMM(NDOFN+1,1)
  871 FORMAT(/,' (MAXIMUM AND MINIMUM OF NODAL DISPLACEMENT) [m] OR [radian]',/,    &
      '           NODE         X-DISP      NODE         Y-DISP      NODE         Z-DISP      NODE        XX-DISP      NODE        YY-DISP      NODE        ZZ-DISP',/,  &
      '  MAX',6(I10,E15.5),/,'  MIN',6(I10,E15.5))
  872 FORMAT(' (MAXIMUM OF VECTOR LENGTH OF NODAL DISPLACEMENT) [m]',/, &
      '           NODE         LENGTH',/,I15,E15.5)
!
! *** ＜残差力（不釣合力）、節点外力、ノルムをファイル＃８に保存する＞
!
!      CALL UNBPRI(FF,IFFIX,NDOFN,NPOIN,NFPOIN,Q)
      CALL UNBPRI(FF,FFF,IFFIX,NDOFN,NPOIN,NFPOIN,Q)
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE UNBPRI(FF,FFF,IFFIX,NDOFN,NPOIN,NFPOIN,Q)
! =====================================================================
! *** 残差力（不釣合力）、節点外力、ノルムをファイル＃８に保存する
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FF(NPOIN*NDOFN)
      DIMENSION FFF(NPOIN*NDOFN)
      DIMENSION IFFIX(NPOIN*NDOFN)
      DIMENSION Q(NPOIN*NDOFN)
      DIMENSION FALL(6,2),RALL(6,2)
      CHARACTER*2 NB(24)
!
! *** 残差力と固定点反力
      WRITE(8,400)
      DO IPOIN=1,NPOIN
        DO K=1,12
          NB(K)=' '
        ENDDO
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IF(IFFIX(ITOTV).EQ.1)THEN
            IPOSI=(IDOFN-1)*2
            NB(IPOSI+1)='('
            NB(IPOSI+2)=')'
          ENDIF
        ENDDO
        WRITE(8,410)IPOIN,(NB(2*IDOFN-1),Q((IPOIN-1)*NDOFN+IDOFN),NB(2*IDOFN),IDOFN=1,NDOFN)
      ENDDO
  400 FORMAT(///,' (FINAL OUT OF BALANCE FORCES AT EACH NODAL POINT) [N]   (   ):BOUNDARY REACTION FORCE',//,  &
            ' NODE       X-UNBAL        Y-UNBAL        Z-UNBAL       XX-UNBAL       YY-UNBAL       ZZ-UNBAL')
  410 FORMAT(I5,6(1X,A1,E12.5,A1))
!
! *** 変形後の節点外力
      WRITE(8,420)(IPOIN,(FF((IPOIN-1)*NDOFN+IDOFN),IDOFN=1,NDOFN),IPOIN=1,NPOIN)
  420 FORMAT(///,' (TOTAL LOAD AT EACH NODAL POINT AFTER DEFORMATION) [N]',//, &
            '  NODE',4X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD',8X,'XX-LOAD',8X,'YY-LOAD',8X,'ZZ-LOAD',/,  &
            (I5,6(F14.4,1X)))
!
! *** 収束状況の検討
      CALL ZEROR2(FALL,NDOFN,2,NDOFN,2) ! FALL：外力の総和
      CALL ZEROR2(RALL,NDOFN,2,NDOFN,2) ! RALL：残差力の総和
      FNORM=0.0
      RNORM=0.0
	  DO IDOFN=1,NDOFN
        DO IPOIN=1,NPOIN
	      ITOTV=(IPOIN-1)*NDOFN+IDOFN
	      IF(IFFIX(ITOTV).EQ.0)THEN
            FALL(IDOFN,1)=FALL(IDOFN,1)+FFF(ITOTV)  ! FALL(IDOFN,1)：外力(自由節点）の総和
            RALL(IDOFN,1)=RALL(IDOFN,1)+Q(ITOTV)    ! RALL(IDOFN,1)：残差力(自由節点）の総和
          ELSE
            FALL(IDOFN,2)=FALL(IDOFN,2)+FFF(ITOTV)  ! FALL(IDOFN,2)：外力(拘束節点）の総和
            RALL(IDOFN,2)=RALL(IDOFN,2)+Q(ITOTV)    ! RALL(IDOFN,2)：残差力(節点）の総和
          ENDIF
        ENDDO
        FNORM=FNORM+FALL(IDOFN,1)**2
        RNORM=RNORM+RALL(IDOFN,1)**2
      ENDDO
!
      FNORM=SQRT(FNORM)                 ! FNORM：外力のノルム
      RNORM=SQRT(RNORM)                 ! RNORM：残差力のノルム
!
      WRITE(8,935) (FALL(IDOFN,1),IDOFN=1,NDOFN),(FALL(IDOFN,2),IDOFN=1,NDOFN),(FALL(IDOFN,1)+FALL(IDOFN,2),IDOFN=1,NDOFN)
      WRITE(8,936) (RALL(IDOFN,1),IDOFN=1,NDOFN)
      WRITE(8,937) (RALL(IDOFN,2),IDOFN=1,NDOFN),(FALL(IDOFN,1)+FALL(IDOFN,2)+RALL(IDOFN,2),IDOFN=1,NDOFN)
      IF(FNORM.NE.0.0)WRITE(8,980) RNORM,FNORM,RNORM/FNORM
  935 FORMAT(//,' (TOTAL EXTERNAL FORCES OF NODAL POINTS) [N] OR [Nm]',/,   &
            ' CONDITION         X-LOAD         Y-LOAD         Z-LOAD      XX-MOMENT      YY-MOMENT      ZZ-MOMENT',/,   &
            '   FREE   ',6E15.5,/,'   FIXED  ',6E15.5,/,100('-'),/,'   TOTAL  ',6E15.5)
  936 FORMAT(/,' (TOTAL RESIDUAL FORCES OF NODAL POINTS) [N] OR [Nm]',/,   &
            ' CONDITION         X-LOAD         Y-LOAD         Z-LOAD      XX-MOMENT      YY-MOMENT      ZZ-MOMENT',/,   &
            '   FREE   ',6E15.5)
  937 FORMAT(/,' (TOTAL REACTION FORCES OF NODAL POINTS) [N] OR [Nm]',/,   &
            ' CONDITION         X-LOAD         Y-LOAD         Z-LOAD      XX-MOMENT      YY-MOMENT      ZZ-MOMENT',/,   &
            '   FIXED  ',6E15.5,/,100('-'),/,'DIFFERENCE',6E15.5)
  980 FORMAT(//,' (NORM OF RESIDUAL FORCES AND EXTERNAL FORCES AT FREE NODAL POINTS) ',/,' RNORM',E15.5,' / FNORM',E15.5,' = ',E15.5,////)
!
      RETURN
      END
!
!
!
! ======================================================================
      SUBROUTINE WRTCOD(INDEX,COORD,NODEC,NELEM,NELEC,NODEM,NPOIN,NELEF,NODEF)
! ======================================================================
! *** ＜形状のデータをファイル＃９に保存する＞ ***
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION NODEM(NELEM,4)   ! 膜要素のデータ（1～3：節点番号、4：部材番号）
      DIMENSION NODEC(NELEC,3)   !ケーブル要素の節点番号(1～2)、部材番号(3)
      DIMENSION NODEF(NELEF,3)   ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
!
      IF(INDEX.EQ.1)THEN    ! INDEX=1（SUBROUTINE INPUT からのCALL）の場合には要素情報を記録する
!
        WRITE(9,*) -3,NPOIN,NPOIN,NELEM,NELEC+NELEF
        DO IELEM=1,NELEM,3
          IIELEM=IELEM+2
          IF(IIELEM.GT.NELEM)IIELEM=NELEM
          WRITE(9,9100)(JELEM,(NODEM(JELEM,INODE),INODE=1,3),JELEM=IELEM,IIELEM)
        ENDDO
 9100   FORMAT(3(4I5,5X))
!
        DO IELEC=1,NELEC,3
          IIELEC=IELEC+2
          IF(IIELEC.GT.NELEC)IIELEC=NELEC
          WRITE(9,9200)(JELT,(NODEC(JELT,INODE),INODE=1,2),JELT=IELEC,IIELEC)
        ENDDO
 9200   FORMAT(3(5X,3I5))
!
        DO IELEF=1,NELEF,3
          IIFRM=IELEF+2
          IF(IIFRM.GT.NELEF)IIFRM=NELEF
          WRITE(9,9200)(IE+NELEC,(NODEF(IE,INODE),INODE=1,2),IE=IELEF,IIFRM)
        ENDDO
      ENDIF
      WRITE(9,9300)(JPOIN,(COORD(JPOIN,IDIR),IDIR=1,3),JPOIN=1,NPOIN)
 9300 FORMAT(I5,3F15.4)
!
      RETURN
      END
!
!
!
! ======================================================================
subroutine echo
! ======================================================================
! *** インプットデータにエラーがあった場合は、残りのデータを出力して終了する。
character*150 text
write(6,*)'***** インプットデータにエラーがあるため残りのインプットデータを出力して終了します *****'
write(8,100)
100 format(/////,'***** 残りのインプットデータを以下に出力します *****',//)
do i=1,10000
  read(7,'(A)',end=999)text
  write(6,'(A)')text
  write(8,'(A)')text
enddo
999 stop
end