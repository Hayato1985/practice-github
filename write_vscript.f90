MODULE WRITE_VSCRIPT_MOD
    ! ======================================================================
    SUBROUTINE WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                             NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,NPROB,NTIMES,STRSF,STRSM,STRSC,TTDIS)
! ======================================================================
! 要素分割図と応力分布図を VECTOR SCRIPT 形式のファイルで出力するためのサブルーチン
!
! サブルーチンの処理を制御するためのパラメータ
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    PARAMETER(M_DIMSWITCH=0)        ! =1のときはＹ座標とＺ座標を入換えて表示する。
!    PARAMETER(M_PRNSTRESS=1)        ! 主応力の表示形式（=0：膜材料Ａ種想定許容応力度区分、=1：ETFEフィルム想定の降伏応力による区分）
    PARAMETER(M_ANGLE=1)            ! 応力の矢印表示方向（=0：IJ方向に平行、=1：主応力方向）
    PARAMETER(M_PALETTE=0)          ! 膜応力の表示色パレットの設定（=0：均等に区分、=1：初期張力状態を想定して区分、=2：安全率4.0を想定して区分）
    PARAMETER(M_ARROWDIM=3)         ! 要素分割図での主軸方向ijを２次元図形で表す（=2）、３次元（=3）
!    PARAMETER(VAL_THICK=100.0)      ! ETFEフィルム想定（M_PRNSTRESS=1）の場合の膜厚[μm]
    PARAMETER(VALST_MAX_FIX=0.0)    ! 複数のデータについて張力矢印の長さを統一したい場合にはVALST_MAX_FIXに0.0以外の数字を設定する
    PARAMETER(TMGNFY=2.0)           ! TMGNFY：張力矢印の長さを調節する係数、2.0の場合には最小内接円の直径に一致する。
    PARAMETER(FMGNFY=1.0000)        ! FMGNFY：図形の拡大倍率、例えばｍ単位をｍｍ単位に変換する場合には1000.0とすればよい。
    PARAMETER(CMGNFY=0.8E+00)       ! CMGNFY：膜要素のIJ方向を表す矢印のIJ長さに対する倍率
!   PARAMETER(VALMARGIN=2.0)        ! VALMARGIN：図形相互の間隔
    PARAMETER( NNODE = 3      )     ! NNODE ：膜要素を構成する節点の数
    DIMENSION COORD(NPOIN,3)        ! COORD(IPOIN,IDIM)：節点座標
    DIMENSION COOR0(NPOIN,3)        ! COOR0(IPOIN,IDIM)：初期形状時の節点座標
    DIMENSION NODEM(NELEM,4)        ! 膜要素のデータ（1～3：節点番号、4：部材番号）
    DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)：膜要素の応力（非抗圧縮性評価）
    DIMENSION NODEC(NELEC,3)        !ケーブル要素の節点番号(1～2)、部材番号(3)
    DIMENSION STRSC(NELEC)          ! STRSC(IELEC)：ケーブルの張力
    DIMENSION NODEF(NELEF,3)        ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
    DIMENSION STRSF(NELEF,NDOFN*2)  ! STRSF(IELEF,IEVAB)：曲げ要素IELEFの部材力（部材座標系における部材端荷重)
    DIMENSION TTDIS(NPOIN*NDOFN)    ! TTDIS(ITOTV)：初期形状に対する変位
    DIMENSION VAL(3),VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(2,3)
    DIMENSION COD_MM(3,2)           ! COD_MM(1, )：座標の最小値、COD_MM(2, )：最大値、、COD_MM(3, )：幅、IDIM=1,2（Ｘ、Ｙ座標）
    DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)：応力の区分、STDIV(0,ISTRE)：最小値、STDIV(9,ISTRE)：最大値、ISTRE=4：主応力の場合
    CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)：カラーパレット
    INTEGER,ALLOCATABLE :: NFILL(:,:)   ! NFILL(IELEM,ISTRE)=IDIV：パレットの番号０～９（０はリンクリングを表す）
    DOUBLE PRECISION,ALLOCATABLE :: COOR1(:,:)  ! COOR1(IPOIN,IDIM)：図形を表示するときに使用する座標値
    DOUBLE PRECISION,ALLOCATABLE :: STRSM1(:,:)    ! STRSM1(IELEM,ISTRE)：膜要素の主応力（非抗圧縮性評価）
!
IF(NPROB.GE.50000.AND.NPROB.LE.50999)THEN   ! ETFEフィルムの材料非線形性に対応する場合（NPROB=50000～50999番の場合）
  M_PRNSTRESS=1        ! 主応力の表示形式（=1：ETFEフィルム想定の降伏応力による区分）
  VAL_THICK=FLOAT(NPROB-50000)
ELSE
  M_PRNSTRESS=0        ! 主応力の表示形式（=0：膜材料Ａ種想定許容応力度区分）
ENDIF
!    OPEN(7,FILE='sample.txt')
!    READ(7,*)NPOIN,NELEM,NELEC,NELEF
!    LOADCREM=1  ! LOADCREM
!    NTIMES=1    ! NTIMES
!==========================================================
ALLOCATE(NFILL(NELEM,4),STAT=IERROR)        ! NFILL(IELEM,ISTRE)=IDIV：パレットの番号０～９（０はリンクリングを表す）、ISTRE=4：主応力の場合
ALLOCATE(COOR1(NPOIN,3),STAT=IERROR)        ! COOR1(IPOIN,IDIM)：図形を表示するときに使用する座標値
ALLOCATE(STRSM1(NELEM,MSTRE),STAT=IERROR)   ! STRSM1(IELEM,ISTRE)：膜要素の主応力（非抗圧縮性評価）
!
!READ(7,*)(I,(COOR0(IPOIN,IDIM),I D IM=1,3),IPOIN=1,NPOIN)
!READ(7,*)(I,(COORD(IPOIN,IDIM),IDIM=1,3),IPOIN=1,NPOIN)
!READ(7,*)(I,(NODEM(IELEM,INODE),INODE=1,3),IELEM=1,NELEM)
!READ(7,*)(I,(STRSM(IELEM,ISTRE),ISTRE=1,3),IELEM=1,NELEM)
!READ(7,*)(I,(NODEC(IELEC,INODE),INODE=1,2),IELEC=1,NELEC)
!READ(7,*)(I,STRSC(IELEC),IELEC=1,NELEC)
!READ(7,*)(I,(NODEF(IELEF,INODE),INODE=1,2),IELEF=1,NELEF)
!READ(7,*)(I,(STRSF(IELEF,IEVAB),IEVAB=1,NDOFN*2),IELEF=1,NELEF)
!==========================================================
! *** 図形を表示するときに使用する座標値COOR1を設定する
  DO IPOIN=1,NPOIN
    DO IDIM=1,3
!     COOR1(IPOIN,IDIM)=COOR0(IPOIN,IDIM)*FMGNFY
      COOR1(IPOIN,IDIM)=COORD(IPOIN,IDIM)*FMGNFY
    ENDDO
    ! M_DIMSWITCH=1のときはＹ座標とＺ座標を入換える。
    IF(M_DIMSWITCH.EQ.1)THEN
      COOR_TEMP=COOR1(IPOIN,2)
      COOR1(IPOIN,2)=COOR1(IPOIN,3)     ! Ｚ←Ｙ
      COOR1(IPOIN,3)=-COOR_TEMP         ! Ｙ←－Ｚ（右手系を保つためにマイナスとする）
    ENDIF
  ENDDO
!==========================================================
! *** Ｘ、Ｙ座標の最小値、最大値、幅の計算
  DO IDIM=1,2
    COD_MM(1,IDIM)=COOR1(1,IDIM)    ! COD_MM(1, )：座標の最小値
    COD_MM(2,IDIM)=COOR1(1,IDIM)    ! COD_MM(2, )：座標の最大値
    DO IPOIN=2,NPOIN
      IF(COOR1(IPOIN,IDIM).LT.COD_MM(1,IDIM))COD_MM(1,IDIM)=COOR1(IPOIN,IDIM)
      IF(COOR1(IPOIN,IDIM).GT.COD_MM(2,IDIM))COD_MM(2,IDIM)=COOR1(IPOIN,IDIM)
    ENDDO
    COD_MM(3,IDIM)=COD_MM(2,IDIM)-COD_MM(1,IDIM)    ! COD_MM(3, )：幅
  ENDDO
  DO IDIM=1,3
    DO JDIM=1,2
      COD_MM(IDIM,JDIM)=COD_MM(IDIM,JDIM)*FMGNFY
    ENDDO
  ENDDO
!==========================================================
! *** 図形相互の間隔、文字サイズの設定
  VALMARGIN=COD_MM(3,1)/6.0
  TEXTSIZE=6.0
!==========================================================
! *** LOADCREM=0 or 1 の場合は、ヘッダーと要素分割図を出力する
IF(LOADCREM.EQ.0.OR.LOADCREM.EQ.1)THEN
!==========================================================
! *** データを出力するファイルの設定
! OPEN(12,FILE='vector_script.txt')
!==========================================================
! *** ヘッダーの出力
  CALL vscript_header()
!==========================================================
! *** LAYER（要素分割図）の出力
  WRITE(12,*) '{Layer Characteristics}'
  WRITE(12,*) ''
  WRITE(12,*) 'Layer(''要素分割図'');'
  WRITE(12,*) 'SetScale(100);'
  WRITE(12,*) 'ShowLayer;'
  WRITE(12,*) 'CopyMode(8);'
  WRITE(12,*) 'LFillFore(0,0,0);'
  WRITE(12,*) 'LFillBack(65535,65535,65535);'
  WRITE(12,*) 'LPenFore(0,0,0);'
  WRITE(12,*) 'LPenBack(65535,65535,65535);'
  WRITE(12,*) ''
  WRITE(12,*) 'Projection(6,0,24.7904,-12.3952,12.3952,12.3952,-12.3952);'
  WRITE(12,*) ''
  WRITE(12,*) '{End of Layer Characteristics}'
  WRITE(12,*) ''
  WRITE(12,*) '{Object Creation Code}'
!==========================================================
! *** 膜要素の出力
  CALL vscript_element_division(1,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** 主軸方向ijを表す線分の出力
  CALL vscript_element_division(2,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** 線材要素の出力
  CALL vscript_element_division(3,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** 節点番号の出力
  CALL vscript_element_division(4,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** 膜要素番号の出力
  CALL vscript_element_division(5,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!========================================================== 
! *** 線材要素番号の出力
  CALL vscript_element_division(6,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
ENDIF
!
!
!==========================================================
! *** LOADCREM≧1 の場合は、応力分布図を出力する
IF(LOADCREM.GE.1) THEN
  DO IPOIN=1,NPOIN
    DO IDIM=1,3
      COOR1(IPOIN,IDIM)=COORD(IPOIN,IDIM)*FMGNFY
    ENDDO
    ! M_DIMSWITCH=1のときはＹ座標とＺ座標を入換える。
    IF(M_DIMSWITCH.EQ.1)THEN
      COOR_TEMP=COOR1(IPOIN,2)
      COOR1(IPOIN,2)=COOR1(IPOIN,3)     ! Ｚ←Ｙ
      COOR1(IPOIN,3)=-COOR_TEMP         ! Ｙ←－Ｚ（右手系を保つためにマイナスとする）
    ENDIF
  ENDDO
!==========================================================
! *** LAYER（case-xx）の出力
  WRITE(12,*) 'SetZVals(0,0);'
  WRITE(12,1007)
  1007 FORMAT(///////////////)
  WRITE(12,*) '{Layer Characteristics}'
  WRITE(12,*) ''
  IF(LOADCREM.LT.10)THEN
    WRITE(12,1003)LOADCREM
    1003 FORMAT('Layer(''case-0',I1,''');')
  ELSE
    WRITE(12,1004)LOADCREM
    1004 FORMAT('Layer(''case-',I2,''');')
  ENDIF
  WRITE(12,*) 'SetScale(100);'
  WRITE(12,*) 'ShowLayer;'
  WRITE(12,*) 'LFillFore(0,0,0);'
  WRITE(12,*) 'LFillBack(65535,65535,65535);'
  WRITE(12,*) 'LPenFore(0,0,0);'
  WRITE(12,*) 'LPenBack(65535,65535,65535);'
  WRITE(12,*) ''
  WRITE(12,*) 'Projection(6,0,24.7904,-12.3952,12.3952,12.3952,-12.3952);'
! WRITE(12,*) 'Projection(0,2,49.5808,-24.7904,24.7904,24.7904,-24.7904);'
  WRITE(12,*) 'SetView(#0｡ 0'' 0" ,#0｡ 0'' 0" ,#0｡ 0'' 0" ,0,0,0);'
  WRITE(12,*) ''
  WRITE(12,*) '{End of Layer Characteristics}'
  WRITE(12,*) ''
  WRITE(12,*) '{Object Creation Code}'
  WRITE(12,*) ''
!==========================================================
IF(NELEM.GT.0)THEN
! *** 膜張力の表示色パレットの設定
  CALL vscript_palette(M_PALETTE,MSTRE,NFILL,FILCOLOR,NELEM,STDIV,STRSM)
!==========================================================
! *** 膜張力（ｘ方向）の出力
  CALL vscript_stress(1,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** 膜張力（ｙ方向）の出力
  CALL vscript_stress(2,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** 膜のせん断力（ｘｙ）の出力
  CALL vscript_stress(3,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** 区分＿膜張力（ｘ方向）の出力
  CALL vscript_stress_div(1,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** 区分＿膜張力（ｙ方向）の出力
  CALL vscript_stress_div(2,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** 区分＿膜のせん断力（ｘｙ）の出力
  CALL vscript_stress_div(3,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** 膜の主応力および張力矢印の出力
  CALL vscript_stress_principal(M_ANGLE,M_PRNSTRESS,VAL_THICK,COD_MM,COOR1,FILCOLOR,MSTRE,NELEM,NFILL,NNODE,NODEM,NPOIN,STDIV,STRSM,STRSM1,TMGNFY,VALMARGIN,VALST_MAX_FIX)
!==========================================================
! *** 区分＿膜の主応力の出力
  CALL vscript_stress_div(4,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
ENDIF
! *** 線材応力（軸力)の出力
IF(NELEC.GT.0)&
  CALL vscript_stress_axialforce(COD_MM,COOR1,NODEC,NDOFN,NELEF,NELEC,NODEF,NPOIN,STRSF,STRSC,TEXTSIZE,VALMARGIN)
!==========================================================
IF(NELEF.GT.0)THEN
! *** 線材応力（１：x軸力、２：yせん断力、３：zせん断力、４：xxねじり、５：yy曲げ、６：zz曲げ)の出力
  !線材の長さの最大値VALLNG_MAXの計算
  VALLNG_MAX=0.0
  DO IELEF=1,NELEF
    IPOIN=NODEF(IELEF,1)
    JPOIN=NODEF(IELEF,2)
    VALLNG=SQRT((COOR1(JPOIN,1)-COOR1(IPOIN,1))**2+(COOR1(JPOIN,2)-COOR1(IPOIN,2))**2+(COOR1(JPOIN,3)-COOR1(IPOIN,3))**2)
    IF(VALLNG_MAX.LT.VALLNG)VALLNG_MAX=VALLNG
  ENDDO
!部材力の最大値FORCE_MEX,MOMENT_MAXの計算
  FORCE_MAX=0.0
  MOMENT_MAX=0.0
  DO IELEF=1,NELEF
    DO ISTRE=1,6
      DO INODE=1,2
        IEVAB=(INODE-1)*6+ISTRE
        IF(ISTRE.LE.3.AND.FORCE_MAX.LT.ABS(STRSF(IELEF,IEVAB)))FORCE_MAX=ABS(STRSF(IELEF,IEVAB))
        IF(ISTRE.GE.4.AND.MOMENT_MAX.LT.ABS(STRSF(IELEF,IEVAB)))MOMENT_MAX=ABS(STRSF(IELEF,IEVAB))
      ENDDO
    ENDDO
  ENDDO
!
  DO ISTRE=1,6
    IF(ISTRE.LE.3)STRSF_MAX=FORCE_MAX
    IF(ISTRE.GE.4)STRSF_MAX=MOMENT_MAX
    IF(STRSF_MAX.EQ.0.0)STRSF_MAX=1.0   ! STRSF_MAX=0の場合のエラー（Zero Divide）対策
IF(ISTRE.EQ.1.OR.ISTRE.EQ.5.OR.ISTRE.EQ.6)THEN  !軸力、強軸曲げ、弱軸曲げのみ出力する
    CALL vscript_stress_frame(ISTRE,COD_MM,COOR0,COOR1,NVARF,NDOFN,NELEF,NODEF,NPOIN,PROPF,STRSF,STRSF_MAX,VALLNG_MAX,TEXTSIZE,TTDIS,VALMARGIN)
ENDIF
  ENDDO
ENDIF
!==========================================================
ENDIF
!==========================================================
! *** LOADCASE=0 or NTIMES(最後)の場合は、クラスを出力してファイルを閉じる
IF(LOADCREM.EQ.0.OR.LOADCREM.EQ.NTIMES)THEN
  CALL vscript_class()    ! クラスの出力
  CLOSE(12)
ENDIF
!
RETURN
END
!
! ======================================================================
    SUBROUTINE vscript_header()    ! ヘッダーの出力
! ======================================================================
WRITE(12,*) 'Procedure LoadFile;'
WRITE(12,*) 'VAR'
WRITE(12,*) '   hatchName:STRING;'
WRITE(12,*) '	result, index:INTEGER;'
WRITE(12,*) '	boolResult:BOOLEAN;'
WRITE(12,*) '    tempHandle, tempHandle1, tempHandle2:HANDLE;'
WRITE(12,*) ''
WRITE(12,*) 'BEGIN'
WRITE(12,*) '{VectorWorks Version 8.0.1}'
WRITE(12,*) ''
WRITE(12,*) '{Global Characteristics}'
WRITE(12,*) ''
WRITE(12,*) 'DrwSize(1,1);'
WRITE(12,*) 'SetUnits(1000000,100,0,0.0254,''m'','' sq m'');'
WRITE(12,*) 'PrimaryUnits(0,100,100,3,1,FALSE,FALSE);'
WRITE(12,*) 'GridLines(1);'
WRITE(12,*) 'PenGrid(1);'
WRITE(12,*) 'DoubLines(0);'
WRITE(12,*) 'SetOriginAbsolute(0,0);'
!WRITE(12,*) 'Snap(TRUE,TRUE,FALSE);'
WRITE(12,*) 'OpenPoly;'
WRITE(12,*) 'SetDashStyle(TRUE,2,0.041656,0.041672);'
WRITE(12,*) 'SetDashStyle(TRUE,2,0.097214,0.041672);'
WRITE(12,*) 'SetDashStyle(TRUE,2,0.208328,0.041672);'
WRITE(12,*) 'SetDashStyle(TRUE,2,0.263885,0.041656);'
WRITE(12,*) 'SetDashStyle(TRUE,2,0.013885,0.027771);'
WRITE(12,*) 'SetDashStyle(TRUE,4,0.125,0.041656,0.013885,0.027786);'
WRITE(12,*) 'SetDashStyle(TRUE,6,0.125,0.041656,0.125,0.041672,0.013885,0.027786);'
WRITE(12,*) 'SetDashStyle(TRUE,6,0.125,0.041656,0.013885,0.027786,0.013885,0.027786);'
WRITE(12,*) 'SetDashStyle(TRUE,4,0.75,0.055557,0.138885,0.055557);'
WRITE(12,*) 'SetDashStyle(TRUE,6,0.75,0.055557,0.138885,0.055557,0.125,0.0625);'
WRITE(12,*) ''
WRITE(12,*) 'Layer(''要素分割図'');'
WRITE(12,*) 'SetScale(200);'
WRITE(12,*) ''
WRITE(12,*) '{End of Global Characteristics}'
WRITE(12,*) ''
WRITE(12,*) ''
WRITE(12,*) '{Record Format Entries}'
WRITE(12,*) ''
WRITE(12,*) '{End of Record Format Entries}'
WRITE(12,*) ''
WRITE(12,*) '{Worksheet Entries}'
WRITE(12,*) ''
WRITE(12,*) '{End of Worksheet Entries}'
WRITE(12,*) ''
WRITE(12,*) '{Hatch Definition Entries}'
WRITE(12,*) 'hatchName:= ''デフォルトハッチング'';'
!WRITE(12,*) 'BeginVectorFill(hatchName,TRUE,FALSE,0);'
WRITE(12,*) 'AddVectorFillLayer(0,0,1,1,0.1767767,-0.1767767,1,1,255);'
WRITE(12,*) 'EndVectorFill;'
WRITE(12,*) ''
WRITE(12,*) '{End of Hatch Definition Entries}'
WRITE(12,*) ''
WRITE(12,*) '{Symbol Library Entries}'
WRITE(12,*) ''
WRITE(12,*) '{End of Symbol Library Entries}'
WRITE(12,*) ''
RETURN
END
! ======================================================================
  SUBROUTINE vscript_element_division(INDEX,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
! ======================================================================
! INDEX=1の場合：膜要素の出力
! INDEX=2の場合：主軸方向ijを表す線分の出力
! INDEX=3の場合：線材要素の出力
! INDEX=4の場合：節点番号の出力
! INDEX=5の場合：膜要素番号の出力
! INDEX=6の場合：線材要素番号の出力
! PARAMETER(M_ARROWDIM=)        ! 要素分割図での主軸方向ijを２次元図形で表す（=2）、３次元（=3）
! PARAMETER(CMGNFY=)            ! CMGNFY：膜要素のIJ方向を表す矢印のIJ長さに対する倍率
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR0(IPOIN,IDIM)：初期形状時の節点座標
  DIMENSION NODEM(NELEM,4)        ! 膜要素のデータ（1～3：節点番号、4：部材番号）
  DIMENSION NODEC(NELEC,3)         !ケーブル要素の節点番号(1～2)、部材番号(3)
  DIMENSION NODEF(NELEF,3)        ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
  DIMENSION ARROW(3,3)            ! 膜要素のIJ方向を表す矢印の始点(1, )、終点(2, )、膜要素の重心点(3, )
  DIMENSION VEC1(2)
! *** 膜要素の出力
IF(INDEX.EQ.1)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''膜_0 要素'');'
    WRITE(12,*) 'FPatByClass;'
    WRITE(12,*) 'FillColorByClass;'
    WRITE(12,*) 'LSByClass;'
    WRITE(12,*) 'PenColorByClass;'
    WRITE(12,*) 'LWByClass;'
    WRITE(12,*) 'MarkerByClass;'
    WRITE(12,*) 'ClosePoly;'
  DO IELEM=1,NELEM
    WRITE(12,*) 'Poly3D('
    WRITE(12,1001)((COOR1(NODEM(IELEM,INODE),IDIM),IDIM=1,3),INODE=1,NNODE)
    WRITE(12,1002)(COOR1(NODEM(IELEM,3),IDIM),IDIM=1,3)
    1001 FORMAT(3(F10.4,','))
    1002 FORMAT(2(F10.4,','),F10.4,/,');')
  ENDDO
!==========================================================
! *** 主軸方向ijを表す線分の出力
ELSEIF(INDEX.EQ.2)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''主軸方向ij'');'
    WRITE(12,*) 'FPatByClass;'
    WRITE(12,*) 'FillColorByClass;'
    WRITE(12,*) 'LSByClass;'
    WRITE(12,*) 'PenColorByClass;'
    WRITE(12,*) 'LWByClass;'
    WRITE(12,*) 'MarkerByClass;'
  DO IELEM=1,NELEM
    DO IDIM=1,3
      ARROW(3,IDIM)=(COOR1(NODEM(IELEM,1),IDIM)+COOR1(NODEM(IELEM,2),IDIM)+COOR1(NODEM(IELEM,3),IDIM))/3    ! 重心点の座標Ｖ3
      ARROW(1,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*COOR1(NODEM(IELEM,1),IDIM)
      ARROW(2,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*COOR1(NODEM(IELEM,2),IDIM)
      ARROW(3,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*ARROW(2,IDIM)     ! 矢じりの長さをＶ32 の(1-m)倍にする
    ENDDO
    IF(M_ARROWDIM.EQ.2)THEN ! ２次元図形で表示する場合 
      WRITE(12,*) 'MoveTo(',ARROW(1,1),',',ARROW(1,2),');'
      WRITE(12,*) 'LineTo(',ARROW(2,1),',',ARROW(2,2),');'
    ELSE    ! ３次元図形で表示する場合（クラス設定が点線なので、Vector Scriptを取込んだ後に属性を実線に変更すれば表示される）
      WRITE(12,*) 'OpenPoly;'
      WRITE(12,*) 'Poly3D('
      WRITE(12,1001)((ARROW(INODE,IDIM),IDIM=1,3),INODE=1,3)
      WRITE(12,1002)(ARROW(3,IDIM),IDIM=1,3)
    ENDIF
  ENDDO
!==========================================================
! *** 線材要素の出力
ELSEIF(INDEX.EQ.3)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''線材-0 要素'');'
    WRITE(12,*) 'FPatByClass;'
    WRITE(12,*) 'FillColorByClass;'
    WRITE(12,*) 'LSByClass;'
    WRITE(12,*) 'PenColorByClass;'
    WRITE(12,*) 'LWByClass;'
    WRITE(12,*) 'MarkerByClass;'
  DO IELEC=1,NELEC
    WRITE(12,*) 'Poly3D('
    WRITE(12,1001)((COOR1(NODEC(IELEC,INODE),IDIM),IDIM=1,3),INODE=1,2)
    WRITE(12,1001)(COOR1(NODEC(IELEC,2),IDIM),IDIM=1,3)
    WRITE(12,1002)(COOR1(NODEC(IELEC,2),IDIM),IDIM=1,3)
  ENDDO
  DO IELEF=1,NELEF
    WRITE(12,*) 'Poly3D('
    WRITE(12,1001)((COOR1(NODEF(IELEF,INODE),IDIM),IDIM=1,3),INODE=1,2)
    WRITE(12,1001)(COOR1(NODEF(IELEF,2),IDIM),IDIM=1,3)
    WRITE(12,1002)(COOR1(NODEF(IELEF,2),IDIM),IDIM=1,3)
  ENDDO
!==========================================================
! *** 節点番号の出力
ELSEIF(INDEX.EQ.4)THEN
  WRITE(12,*) 'NameClass(''番号-0 節点'');'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillFore(0,0,0);'
  WRITE(12,*) 'FillBack(65535,65535,13107);'
  WRITE(12,*) 'PenFore(0,0,0);'
  WRITE(12,*) 'PenBack(65535,65535,65535);'
  WRITE(12,*) 'TextFont(GetFontID(''ＭＳ ゴシック''));'
  WRITE(12,*) 'TextSize(',INT(TEXTSIZE/1.5),');'
  WRITE(12,*) 'TextFace([]);'
  WRITE(12,*) 'TextFlip(0);'
  WRITE(12,*) 'TextRotate(0);'
  WRITE(12,*) 'TextSpace(2);'
  WRITE(12,*) 'TextJust(1);'
  WRITE(12,*) 'TextVerticalAlign(1);'
  DO IPOIN=1,NPOIN
    WRITE(12,*) 'TextOrigin(',COOR1(IPOIN,1),',',COOR1(IPOIN,2),');'
    WRITE(12,*) 'BeginText;'
    IF(IPOIN.LT.100)THEN
      WRITE(12,1012) IPOIN
    ELSEIF(IPOIN.LT.1000)THEN
      WRITE(12,1013) IPOIN
    ELSE
      WRITE(12,1014) IPOIN
    ENDIF
    1012 FORMAT(1H',I2,1H')
    1013 FORMAT(1H',I3,1H')
    1014 FORMAT(1H',I4,1H')
    WRITE(12,*) 'EndText;'
  ENDDO
!==========================================================
! *** 膜要素番号の出力
ELSEIF(INDEX.EQ.5)THEN
  WRITE(12,*) 'NameClass(''番号-1 膜要素'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'PenFore(0,0,54272);'
  DO IELEM=1,NELEM
    VEC1(1)=(COOR1(NODEM(IELEM,1),1)+COOR1(NODEM(IELEM,2),1)+COOR1(NODEM(IELEM,3),1))/3.0   !VEC1()： 三角形の重心の座標
    VEC1(2)=(COOR1(NODEM(IELEM,1),2)+COOR1(NODEM(IELEM,2),2)+COOR1(NODEM(IELEM,3),2))/3.0   !VEC1()： 三角形の重心の座標
    WRITE(12,*) 'TextOrigin(',VEC1(1),',',VEC1(2),');'
    WRITE(12,*) 'BeginText;'
    IF(IELEM.LT.100)THEN
      WRITE(12,1012) IELEM
    ELSEIF(IELEM.LT.1000)THEN
      WRITE(12,1013) IELEM
    ELSE
      WRITE(12,1014) IELEM
    ENDIF
    WRITE(12,*) 'EndText;'
  ENDDO
!========================================================== 
! *** 線材要素番号の出力
ELSEIF(INDEX.EQ.6)THEN
  WRITE(12,*) 'NameClass(''番号-2 線材'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'PenFore(56797,0,0);'
! ケーブル要素（軸力要素）の場合
  DO IELEC=1,NELEC
    VEC1(1)=(COOR1(NODEC(IELEC,1),1)+COOR1(NODEC(IELEC,2),1))/2.0   !VEC1()： ケーブル要素の中点の座標
    VEC1(2)=(COOR1(NODEC(IELEC,1),2)+COOR1(NODEC(IELEC,2),2))/2.0   !VEC1()： ケーブル要素の中点の座標
    WRITE(12,*) 'TextOrigin(',VEC1(1),',',VEC1(2),');'
    WRITE(12,*) 'BeginText;'
    IF(IELEC.LT.100)THEN
      WRITE(12,1012) IELEC
    ELSEIF(IELEC.LT.1000)THEN
      WRITE(12,1013) IELEC
    ELSE
      WRITE(12,1014) IELEC
    ENDIF
    WRITE(12,*) 'EndText;'
  ENDDO
! 曲げ要素の場合（カギ括弧つきの番号で表示 <1> ）
  DO IELEF=1,NELEF
    VEC1(1)=(COOR1(NODEF(IELEF,1),1)+COOR1(NODEF(IELEF,2),1))/2.0   !VEC1()： 曲げ要素の中点の座標
    VEC1(2)=(COOR1(NODEF(IELEF,1),2)+COOR1(NODEF(IELEF,2),2))/2.0   !VEC1()： 曲げ要素の中点の座標
    WRITE(12,*) 'TextOrigin(',VEC1(1),',',VEC1(2),');'
    WRITE(12,*) 'BeginText;'
    IF(IELEF.LT.100)THEN
      WRITE(12,1022) IELEF
    ELSEIF(IELEF.LT.1000)THEN
      WRITE(12,1023) IELEF
    ELSE
      WRITE(12,1024) IELEF
    ENDIF
    WRITE(12,*) 'EndText;'
    1022 FORMAT(2H'<,I2,2H>')
    1023 FORMAT(2H'<,I3,2H>')
    1024 FORMAT(2H'<,I4,2H>')
  ENDDO
ENDIF
RETURN
END
! ======================================================================
    SUBROUTINE vscript_palette(M_PALETTE,MSTRE,NFILL,FILCOLOR,NELEM,STDIV,STRSM)
! ======================================================================
! PARAMETER(M_PALETTE=0)          ! 膜応力の表示色パレットの設定（=0：均等に区分、=1：初期張力状態を想定して区分、=2：安全率4.0を想定して区分）
IMPLICIT DOUBLEPRECISION (A-H,O-Z)
DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)：膜要素の張力（非抗圧縮性評価）
DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)：応力の区分、STDIV(0,ISTRE)：最小値、STDIV(9,ISTRE)：最大値、ISTRE=4：主応力の場合
INTEGER NFILL(NELEM,4)          ! NFILL(IELEM,ISTRE)=IDIV：パレットの番号０～９（０はリンクリングを表す）、ISTRE=4：主応力の場合
CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)：カラーパレット
!
! カラーパレットの設定
  IF(M_PALETTE.NE.1)THEN    ! M_PALETTE=0：均等に区分 OR INDEX=2：安全率4.0を想定して区分
    FILCOLOR(9) = '(56797,00000,00000)'     ! STDIV(8, )＜    ≦STDIV(9, )
    FILCOLOR(8) = '(65535,26214,13107)'     ! STDIV(7, )＜    ≦STDIV(8, )
    FILCOLOR(7) = '(65535,65535,13107)'     ! STDIV(6, )＜    ≦STDIV(7, )
    FILCOLOR(6) = '(00000,34952,00000)'     ! STDIV(5, )＜    ≦STDIV(6, )
    FILCOLOR(5) = '(00000,65535,00000)'     ! STDIV(4, )＜    ≦STDIV(5, )
    FILCOLOR(4) = '(00000,00000,65535)'     ! STDIV(3, )＜    ≦STDIV(4, )
    FILCOLOR(3) = '(42662,51914,61680)'     ! STDIV(2, )＜    ≦STDIV(3, )
    FILCOLOR(2) = '(34952,34952,34952)'     ! STDIV(1, )＜    ≦STDIV(2, )
    FILCOLOR(1) = '(48059,48059,48059)'     ! STDIV(0, )≦    ≦STDIV(1, )
    FILCOLOR(0) = '(65535,65535,65535)'     ! リンクリング
  ELSE  ! M_PALETTE=1：初期張力状態を想定して区分
    FILCOLOR(9) = '(56797,00000,00000)'
    FILCOLOR(8) = '(65535,26214,13107)'
    FILCOLOR(7) = '(65535,65535,13107)'
    FILCOLOR(6) = '(00000,56797,00000)'
    FILCOLOR(5) = '(00000,43690,00000)'
    FILCOLOR(4) = '(00000,21845,00000)'
    FILCOLOR(3) = '(13107,21845,65535)'
    FILCOLOR(2) = '(00000,00000,65535)'
    FILCOLOR(1) = '(00000,00000,32896)'
    FILCOLOR(0) = '(65535,65535,65535)'
  ENDIF
!
! 膜応力の最大値、最小値を求める
  DO IELEM=1,NELEM
    DO ISTRE=1,2
      IF(IELEM.EQ.1)THEN
        STDIV(9,ISTRE)=STRSM(1,ISTRE)    ! STDIV(9,ISTRE)：膜応力の最大値
        STDIV(0,ISTRE)=STRSM(1,ISTRE)    ! STDIV(0,ISTRE)：膜応力の最小値
      ELSE
        IF(STDIV(9,ISTRE).LT.STRSM(IELEM,ISTRE))STDIV(9,ISTRE)=STRSM(IELEM,ISTRE)
        IF(STDIV(0,ISTRE).GT.STRSM(IELEM,ISTRE))STDIV(0,ISTRE)=STRSM(IELEM,ISTRE)
      ENDIF
    ENDDO
    IF(IELEM.EQ.1)THEN
      STDIV(9,3)=ABS(STRSM(1,3))    ! STDIV(9,3)：せん断力の最大値
      STDIV(0,3)=ABS(STRSM(1,3))    ! STDIV(0,3)：せん断力の最小値
    ELSE
      IF(STDIV(9,3).LT.ABS(STRSM(IELEM,3)))STDIV(9,3)=ABS(STRSM(IELEM,3))
      IF(STDIV(0,3).GT.ABS(STRSM(IELEM,3)))STDIV(0,3)=ABS(STRSM(IELEM,3))
    ENDIF
  ENDDO
!
! 膜応力の区分の基準張力を決める(9分割）
  IF(M_PALETTE.EQ.0)THEN
!   均等に区分
    DO IDIV=1,8
      DO ISTRE=1,3
        STDIV(IDIV,ISTRE)=STDIV(0,ISTRE)+(STDIV(9,ISTRE)-STDIV(0,ISTRE))/16.0*(2*IDIV-1)
      ENDDO
    ENDDO
  ELSEIF(M_PALETTE.EQ.1)THEN
!   初期張力状態をｘ，ｙ方向とも同じ区分で指定して表示したい場合の処理
    IF(STDIV(9,1).GT.STDIV(9,2))THEN
      STDIV(9,2)=STDIV(9,1)
    ELSE
      STDIV(9,1)=STDIV(9,2)
    ENDIF
    STDIV(1,1)=200.0*0.9     ! Tx(1) = 200# * 0.9
    STDIV(2,1)=200.0*1.1     ! Tx(2) = 200# * 1.1
    STDIV(1,2)=STDIV(1,1)
    STDIV(2,2)=STDIV(2,2)
    DO IDIV=3,8
      STDIV(IDIV,1)=STDIV(IDIV-1,1)+(STDIV(9,1)-200.0*1.1)/7.0
      STDIV(IDIV,2)=STDIV(IDIV,1)
    ENDDO
  ELSEIF(M_PALETTE.EQ.2)THEN
!   安全率4.0を区切りとして表示したい場合の処理
    STDIV(6,1)=3750.0    ! Tx(6) = 3750    '膜材料Ａ種：引張強度450kgf/3cm相当
    STDIV(6,2)=3000.0    ! Ty(6) = 3000    '膜材料Ａ種：引張強度360kgf/3cm相当
    IF(STDIV(0,1).LT.0.0.OR.STDIV(0,2).LT.0.0)THEN
      STDIV(1,1)=0.0
      STDIV(1,2)=0.0
    ELSE
      STDIV(1,1)=STDIV(0,1)+(STDIV(6,1)-STDIV(0,1))/6.0
      STDIV(1,2)=STDIV(0,2)+(STDIV(6,2)-STDIV(0,2))/6.0
    ENDIF
    DO IDIV=2,8
      IF(IDIV.LT.6)THEN
        STDIV(IDIV,1)=STDIV(1,1)+(STDIV(6,1)-STDIV(1,1))/5.0*(IDIV-1)
        STDIV(IDIV,2)=STDIV(1,2)+(STDIV(6,2)-STDIV(1,2))/5.0*(IDIV-1)
      ELSEIF(IDIV.GT.6)THEN
        STDIV(IDIV,1)=STDIV(6,1)+(STDIV(9,1)-STDIV(6,1))/3.0*(IDIV-6)
        STDIV(IDIV,2)=STDIV(6,2)+(STDIV(9,2)-STDIV(6,2))/3.0*(IDIV-6)
      ENDIF
    ENDDO
  ENDIF
!
! 段階ごとに膜要素の色を決める(10段階）
  DO IELEM=1,NELEM
    DO ISTRE=1,3
      STRESS=STRSM(IELEM,ISTRE)
      IF(ISTRE.EQ.3)STRESS=ABS(STRESS)
      DO IDIV=0,8
        IF(STRESS.EQ.0.0)THEN
          NFILL(IELEM,ISTRE)=0
          EXIT
        ELSEIF(STRESS.EQ.STDIV(0,ISTRE))THEN
          NFILL(IELEM,ISTRE)=1
          EXIT
        ELSEIF(STRESS.GT.STDIV(IDIV,ISTRE).AND.STRESS.LE.STDIV(IDIV+1,ISTRE))THEN
          NFILL(IELEM,ISTRE)=IDIV+1
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!
  RETURN
  END
! ======================================================================
    SUBROUTINE vscript_stress(ISTRE,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,WIDTH_Y)
! ======================================================================
! 膜応力（x方向、y方向張力、せん断力）の出力
! PARAMETER(VALMARGIN=)            ! VALMARGIN：図形相互の間隔
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)：節点座標
  DIMENSION NODEM(NELEM,4)        ! 膜要素のデータ（1～3：節点番号、4：部材番号）
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)：カラーパレット
  INTEGER NFILL(NELEM,4)        ! NFILL(IELEM,ISTRE)=IDIV：パレットの番号０～９（０はリンクリングを表す）、ISTRE=4：主応力の場合
  DIMENSION VEC(3,3)              ! VEC(INUM,IDIM)：表示用の節点座標
! ISTRE：=1 張力（x方向）、=2 張力（y方向）、=3 せん断力xy
! WIDTH_Y：図形のＹ方向の幅
!
  WRITE(12,*) ''
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''膜_1 張力（x方向）'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''膜_2 張力（y方向）'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''膜_3 せん断力xy'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'ClosePoly;'
  DO IELEM=1,NELEM
    WRITE(12,*) 'FillBack',FILCOLOR(NFILL(IELEM,ISTRE)),';'
    WRITE(12,*) 'Poly3D('
    DO INODE=1,NNODE
      IPOIN=NODEM(IELEM,INODE)
      VEC(INODE,1)=COOR1(IPOIN,1)
      VEC(INODE,2)=COOR1(IPOIN,2)-(WIDTH_Y+VALMARGIN)*(ISTRE-1)
      VEC(INODE,3)=COOR1(IPOIN,3)
    ENDDO
    WRITE(12,1001)((VEC(INODE,IDIM),IDIM=1,3),INODE=1,NNODE)
    WRITE(12,1002)(VEC(3,IDIM),IDIM=1,3)
    1001 FORMAT(3(F10.4,','))
    1002 FORMAT(2(F10.4,','),F10.4,/,');')
  ENDDO
  RETURN
  END
! ======================================================================
  SUBROUTINE vscript_stress_principal(M_ANGLE,M_PRNSTRESS,VAL_THICK,COD_MM,COOR1,FILCOLOR,MSTRE,NELEM,NFILL,NNODE,NODEM,NPOIN,STDIV,STRSM,STRSM1,TMGNFY,VALMARGIN,VALST_MAX_FIX)
! ======================================================================
! *** 膜の主応力および張力矢印の出力
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)：節点座標
  DIMENSION NODEM(NELEM,4)        ! 膜要素のデータ（1～3：節点番号、4：部材番号）
  DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)：膜要素の張力（非抗圧縮性評価）
  DIMENSION STRSM1(NELEM,MSTRE)   ! STRSM1(IELEM,ISTRE)：膜要素の主応力（非抗圧縮性評価）
  DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)：応力の区分、STDIV(0,ISTRE)：最小値、STDIV(9,ISTRE)：最大値、ISTRE=4：主応力の場合
  DIMENSION VAL(3),VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC(3,3)
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )：座標の最小値、COD_MM(2, )：最大値、、COD_MM(3, )：幅、IDIM=1,2（Ｘ、Ｙ座標）
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)：カラーパレット
  INTEGER NFILL(NELEM,4)          ! NFILL(IELEM,ISTRE)=IDIV：パレットの番号０～９（０はリンクリングを表す）、ISTRE=4：主応力の場合
! PARAMETER(M_PRNSTRESS=)         ! 主応力の表示形式（=0：膜材料Ａ種想定許容応力度区分、=1：ETFEフィルム想定の降伏応力による区分）
! PARAMETER(M_ANGLE=)             ! 応力の矢印表示方向（=0：IJ方向に平行、=1：主応力方向）
! PARAMETER(VAL_THICK=)           ! ETFEフィルム想定（M_PRNSTRESS=1）の場合の膜厚[μm]
! PARAMETER(VALMARGIN=)           ! VALMARGIN：図形相互の間隔
!
  OFFSET_Y=-(COD_MM(3,2)+VALMARGIN)*3     ! OFFSET_Y：図形表示位置のＹ方向移動量
!================================================
! 膜要素の主応力STRSM1(IELEM,ISTRE)を求める
  DO IELEM=1,NELEM
    PX=STRSM(IELEM,1)
    PY=STRSM(IELEM,2)
    PXY=STRSM(IELEM,3)
    PXPY=(PX+PY)*0.5
    SQ=SQRT(((PX-PY)*(PX-PY))*0.25+PXY*PXY)
    STRSM1(IELEM,1)=PXPY+SQ
    STRSM1(IELEM,2)=PXPY-SQ
    IF(PX.LT.0.0.OR.PY.LT.0.0.OR.PX.EQ.PY)THEN
      STRSM1(IELEM,3)=0.0
    ELSE
      STRSM1(IELEM,3)=28.6*ATAN((2.0*PXY)/(PX-PY))
    ENDIF
    IF(PX.LE.PY)THEN
      IF(PXY.GE.0.0)STRSM1(IELEM,3)=STRSM1(IELEM,3)+90.0
      IF(PXY.LT.0.0)STRSM1(IELEM,3)=STRSM1(IELEM,3)-90.0
    ENDIF
  ENDDO
!================================================
! VALR_MAX：内接円の半径の最大値の計算
  VALR_MAX=0.0
  DO IELEM=1,NELEM
    DO ISIDE=1,3
      N1=NODEM(IELEM,ISIDE)
      IF(ISIDE.LT.3)N2=NODEM(IELEM,ISIDE+1)
      IF(ISIDE.EQ.3)N2=NODEM(IELEM,1)
      VAL(ISIDE)=SQRT((COOR1(N2,1)-COOR1(N1,1))**2+(COOR1(N2,2)-COOR1(N1,2))**2+(COOR1(N2,3)-COOR1(N1,3))**2)
    ENDDO
    VALS=0.5*(VAL(1)+VAL(2)+VAL(3))
    VALSEC=SQRT(VALS*(VALS-VAL(1))*(VALS-VAL(2))*(VALS-VAL(3)))     ! VALSEC：三角形の面積
    VALR=VALSEC/VALS   ! VALR：内接円の半径
    IF(VALR.GT.VALR_MAX)VALR_MAX=VALR
  ENDDO
! VALST_MAX：膜張力の最大値の計算
  IF(VALST_MAX_FIX.NE.0.0)THEN
    VALST_MAX=VALST_MAX_FIX     ! 複数のデータについて矢印の長さを統一したい場合にはVALST_MAXを指定値とする。
  ELSE
    VALST_MAX=0.0
    DO IELEM=1,NELEM
      IF(M_ANGLE.EQ.0)THEN      ! IJ方向に平行に矢印を出力する場合
        DO ISTRE=1,2
          IF(VALST_MAX.LT.STRSM(IELEM,ISTRE))VALST_MAX=STRSM(IELEM,ISTRE)
        ENDDO
      ELSEIF(M_ANGLE.EQ.1)THEN  ! 主応力方向の矢印を出力する場合
        DO ISTRE=1,2
          IF(VALST_MAX.LT.STRSM1(IELEM,ISTRE))VALST_MAX=STRSM1(IELEM,ISTRE)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
! VALMFY：張力矢印の倍率の計算
  IF(VALST_MAX.NE.0.0)VALMFY=VALR_MAX/VALST_MAX*TMGNFY
! 張力矢印の出力
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''膜_4 張力矢印'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  DO IELEM=1,NELEM
    WRITE(12,*) 'MarkerByClass;'
!   WRITE(12,*) 'Marker(11,0.047241,8);'    ! 両矢印←→、長さ1.2mm/25.4=0.047241、角度8°
    N1=NODEM(IELEM,1)
    N2=NODEM(IELEM,2)
    N3=NODEM(IELEM,3)
    VALLNG=SQRT((COOR1(N2,1)-COOR1(N1,1))**2+(COOR1(N2,2)-COOR1(N1,2))**2+(COOR1(N2,3)-COOR1(N1,3))**2) ! VALLNG：辺１２の長さ
    DO IDIM=1,3
      VEC1(IDIM)=(COOR1(N1,IDIM)+COOR1(N2,IDIM)+COOR1(N3,IDIM))/3.0   !VEC1()： 三角形の重心の座標
      VEC2(IDIM)=(COOR1(N2,IDIM)-COOR1(N1,IDIM))/VALLNG     ! VEC2()：辺１２の単位ベクトル
    ENDDO
    IF(M_ANGLE.EQ.0)THEN      ! IJ方向に平行に矢印を出力する場合
      DO IDIM=1,3
        VEC3(IDIM)=VEC2(IDIM)*STRSM(IELEM,1)*VALMFY   ! VEC3()：辺１２方向の張力矢印ベクトル
      ENDDO
     VEC4(1)=-VEC2(2)*STRSM(IELEM,2)*VALMFY          ! VEC4()：辺１２に直交方向の張力矢印ベクトル
     VEC4(2)= VEC2(1)*STRSM(IELEM,2)*VALMFY
    ELSEIF(M_ANGLE.EQ.1)THEN  ! 主応力方向の矢印を出力する場合
      ANGLE1=STRSM1(IELEM,3)*3.14/180.0
      ANGLE2=(STRSM1(IELEM,3)+90.0)*3.14/180.0
      VEC3(1)=(VEC2(1)*COS(ANGLE1)-VEC2(2)*SIN(ANGLE1))*STRSM1(IELEM,1)*VALMFY   ! VEC3()：辺１２方向の張力矢印ベクトル
      VEC3(2)=(VEC2(1)*SIN(ANGLE1)+VEC2(2)*COS(ANGLE1))*STRSM1(IELEM,1)*VALMFY   ! VEC3()：辺１２方向の張力矢印ベクトル
      VEC4(1)=(VEC2(1)*COS(ANGLE2)-VEC2(2)*SIN(ANGLE2))*STRSM1(IELEM,2)*VALMFY   ! VEC4()：辺１２に直交方向の張力矢印ベクトル
      VEC4(2)=(VEC2(1)*SIN(ANGLE2)+VEC2(2)*COS(ANGLE2))*STRSM1(IELEM,2)*VALMFY   ! VEC4()：辺１２に直交方向の張力矢印ベクトル
    ENDIF
    WRITE(12,*) 'MoveTo(',VEC1(1)-VEC3(1)*0.5,',',VEC1(2)-VEC3(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'LineTo(',VEC1(1)+VEC3(1)*0.5,',',VEC1(2)+VEC3(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'MoveTo(',VEC1(1)-VEC4(1)*0.5,',',VEC1(2)-VEC4(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'LineTo(',VEC1(1)+VEC4(1)*0.5,',',VEC1(2)+VEC4(2)*0.5+OFFSET_Y,');'
  ENDDO
!================================================
! ETFEフィルム想定の場合>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  IF(M_PRNSTRESS.EQ.1)THEN
    DO IELEM=1,NELEM
!     相当応力を求める
      PX=STRSM(IELEM,1)
      PY=STRSM(IELEM,2)
      PXY=STRSM(IELEM,3)
      STRSM1(IELEM,1)=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)
!     張力[kg/m]を応力[N/mm2]に変換する
      STRSM1(IELEM,1)=STRSM1(IELEM,1)*9.81/1000.0/(VAL_THICK*1.0E-3)
      STRSM1(IELEM,2)=STRSM1(IELEM,2)*9.81/1000.0/(VAL_THICK*1.0E-3)
    ENDDO
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! 膜の主応力STRSM1(IELEM,1)の最大値、最小値を求める
  DO IELEM=1,NELEM
    IF(IELEM.EQ.1)THEN
      STDIV(9,4)=STRSM1(1,1)        ! STDIV(9,4)：主応力の最大値
      STDIV(0,4)=STRSM1(1,1)        ! STDIV(0,4)：主応力の最小値
    ELSE
      IF(STDIV(9,4).LT.STRSM1(IELEM,1))STDIV(9,4)=STRSM1(IELEM,1)
      IF(STDIV(0,4).GT.STRSM1(IELEM,1))STDIV(0,4)=STRSM1(IELEM,1)
    ENDIF
  ENDDO
!
! 主応力の区分を求める
    IF(M_PRNSTRESS.EQ.0)THEN
!       安全率4.0を区切りとして表示したい場合の処理
        STDIV(6,4)=3000.0    ! 膜材料Ａ種：引張強度360kgf/3cm相当
      IF(STDIV(0,4).LT.0.0)THEN
        STDIV(1,4)=0.0
      ELSE
        STDIV(1,4)=STDIV(0,4)+(STDIV(6,4)-STDIV(0,4))/6.0
      ENDIF
      DO IDIV=2,8
        IF(IDIV.LT.6)THEN
          STDIV(IDIV,4)=STDIV(1,4)+(STDIV(6,4)-STDIV(1,4))/5.0*(IDIV-1)
        ELSEIF(IDIV.GT.6)THEN
          STDIV(IDIV,4)=STDIV(6,4)+(STDIV(9,4)-STDIV(6,4))/3.0*(IDIV-6)
        ENDIF
      ENDDO
    ELSEIF(M_PRNSTRESS.EQ.1)THEN
!     ETFEフィルム想定の場合
      STDIV(4,4)=13.0    ! 13 [N/mm2]
      STDIV(8,4)=20.0    ! 13 [N/mm2]
      IF(STDIV(0,4).LT.0.0)THEN
        STDIV(1,4)=0.0
      ELSE
        STDIV(1,4)=STDIV(0,4)+(STDIV(4,4)-STDIV(0,4))/4.0
      ENDIF
      DO IDIV=2,8
        IF(IDIV.LT.4)THEN
          STDIV(IDIV,4)=STDIV(1,4)+(STDIV(4,4)-STDIV(1,4))/3.0*(IDIV-1)
        ELSEIF(IDIV.LT.8)THEN
          STDIV(IDIV,4)=STDIV(4,4)+(STDIV(8,4)-STDIV(4,4))/4.0*(IDIV-4)
        ENDIF
      ENDDO
    ELSE
!     均等に区分する
      DO IDIV=1,8
        STDIV(IDIV,4)=STDIV(0,4)+(STDIV(9,4)-STDIV(0,4))/16.0*(2*IDIV-1)
      ENDDO
    ENDIF
!
! 段階ごとに膜要素の色を決める(10段階）
  DO IELEM=1,NELEM
    STRESS=STRSM1(IELEM,1)
    DO IDIV=0,8
      IF(STRESS.EQ.0.0)THEN
        NFILL(IELEM,4)=0
        EXIT
      ELSEIF(STRESS.EQ.STDIV(0,4))THEN
        NFILL(IELEM,4)=1
        EXIT
      ELSEIF(STRESS.GT.STDIV(IDIV,4).AND.STRESS.LE.STDIV(IDIV+1,4))THEN
        NFILL(IELEM,4)=IDIV+1
        EXIT
      ENDIF
    ENDDO
  ENDDO
!
! 膜要素を出力する
  DO IELEM=1,NELEM
    WRITE(12,*) 'FillBack',FILCOLOR(NFILL(IELEM,4)),';'
    WRITE(12,*) 'Poly3D('
    DO INODE=1,NNODE
      IPOIN=NODEM(IELEM,INODE)
      VEC(INODE,1)=COOR1(IPOIN,1)
      VEC(INODE,2)=COOR1(IPOIN,2)+OFFSET_Y
      VEC(INODE,3)=COOR1(IPOIN,3)
    ENDDO
    WRITE(12,1001)((VEC(INODE,IDIM),IDIM=1,3),INODE=1,NNODE)
    WRITE(12,1002)(VEC(3,IDIM),IDIM=1,3)
    1001 FORMAT(3(F10.4,','))
    1002 FORMAT(2(F10.4,','),F10.4,/,');')
  ENDDO
  RETURN
  END
! ======================================================================
  SUBROUTINE vscript_stress_div(ISTRE,FILCOLOR,COD_MM,STDIV,VALMARGIN)
! ======================================================================
! 応力区分の出力
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)：応力の区分、STDIV(0,ISTRE)：最小値、STDIV(9,ISTRE)：最大値、ISTRE=4：主応力の場合
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )：座標の最小値、COD_MM(2, )：最大値、、COD_MM(3, )：幅、IDIM=1,2（Ｘ、Ｙ座標）
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)：カラーパレット
  DIMENSION VEC(2,2)
! ISTRE：=1 張力（x方向）、=2 張力（y方向）、=3 せん断力xy、=4 主応力
! VALMARGIN：図形相互の間隔
!
! ＢＯＸサイズの設定
  TEXTSIZE=6.0
  IF(VALMARGIN*10.0.LE.COD_MM(3,2))THEN ! ＢＯＸの高さが図形の縦方向長さよりも小さい場合
    BOXSIZE=VALMARGIN*0.5
  ELSE
    BOXSIZE=COD_MM(3,2)/10.0
    TEXTSIZE=TEXTSIZE*BOXSIZE/VALMARGIN
  ENDIF
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''区分-膜_1 張力x'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''区分-膜_2 張力y'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''区分-膜_3 せん断力xy'');'
  IF(ISTRE.EQ.4)WRITE(12,*) 'NameClass(''膜_4 張力矢印'');'
  WRITE(12,*) 'BeginGroup;'
! 応力区分のＢＯＸ
  WRITE(12,*) 'FillPat(1);'
  DO IDIV=0,9   
    VEC(1,1)=COD_MM(2,1)+VALMARGIN
    VEC(2,1)=COD_MM(2,1)+VALMARGIN+BOXSIZE
    VEC(1,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*IDIV
    VEC(2,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*(IDIV+1)
    WRITE(12,*) 'FillBack',FILCOLOR(IDIV),';'
    WRITE(12,*) 'Rect(',VEC(1,1),',',VEC(1,2),',',VEC(2,1),',',VEC(2,2),');'
  ENDDO
! 応力区分の数字 
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillBack(0,0,0);'
  WRITE(12,*) 'TextSize(',INT(TEXTSIZE/1.5),');'
  WRITE(12,*) 'TextRotate(0);'
  DO IDIV=0,9
    VEC(1,1)=COD_MM(2,1)+VALMARGIN+BOXSIZE*1.2
    VEC(1,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*(IDIV+1)
    WRITE(12,*) 'TextOrigin(',VEC(1,1),',',VEC(1,2),');'
    WRITE(12,*) 'BeginText;'
    IF(IDIV.EQ.0)THEN
      WRITE(12,*) '''  WRINKLING'''
    ELSE
      WRITE(12,1030) INT(STDIV(IDIV-1,ISTRE)),INT(STDIV(IDIV,ISTRE))
      1030 FORMAT(1H',I4,' ～',I4,1H')
    ENDIF
    WRITE(12,*) 'EndText;'
  ENDDO
! タイトルの出力
  VEC(1,1)=COD_MM(2,1)+VALMARGIN+BOXSIZE*0.8
  VEC(1,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*11
  WRITE(12,*) 'TextOrigin(',VEC(1,1),',',VEC(1,2),');'
  WRITE(12,*) 'BeginText;'
  IF(ISTRE.EQ.1)WRITE(12,*) ''' TENSILE STRESS -x- '''
  IF(ISTRE.EQ.2)WRITE(12,*) ''' TENSILE STRESS -y- '''
  IF(ISTRE.EQ.3)WRITE(12,*) ''' SHEAR STRESS -xy- '''
  IF(ISTRE.EQ.4)WRITE(12,*) ''' PRINCIPAL STRESS '''
  WRITE(12,*) 'EndText;'
  WRITE(12,*) 'EndGroup;'
  RETURN
  END
! ======================================================================
  SUBROUTINE vscript_stress_axialforce(COD_MM,COOR1,NODEC,NDOFN,NELEF,NELEC,NODEF,NPOIN,STRSF,STRSC,TEXTSIZE,VALMARGIN)
! ======================================================================
! *** 線材応力（軸力)の出力
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)：図形を表示するときに使用する座標値
  DIMENSION NODEC(NELEC,3)        !ケーブル要素の節点番号(1～2)、部材番号(3)
  DIMENSION STRSC(NELEC)          ! STRSC(IELEC)：ケーブルの張力
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )：座標の最小値、COD_MM(2, )：最大値、、COD_MM(3, )：幅、IDIM=1,2（Ｘ、Ｙ座標）
  DIMENSION VEC1(2)
! PARAMETER(VALMARGIN=)        ! VALMARGIN：図形相互の間隔
!
  OFFSET_Y=-(COD_MM(3,2)+VALMARGIN)*4     ! OFFSET_Y：図形表示位置のＹ方向移動量
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''線材-1 x軸力'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
! 部材の出力
  DO INUM=1,NELEC
    WRITE(12,*) 'Poly3D('
    DO INODE=1,2
      IPOIN=NODEC(INUM,INODE)
      WRITE(12,1001)COOR1(IPOIN,1),COOR1(IPOIN,2)+OFFSET_Y,COOR1(IPOIN,3)
    ENDDO
    WRITE(12,1001)COOR1(IPOIN,1),COOR1(IPOIN,2)+OFFSET_Y,COOR1(IPOIN,3)
    WRITE(12,1002)COOR1(IPOIN,1),COOR1(IPOIN,2)+OFFSET_Y,COOR1(IPOIN,3)
    1001 FORMAT(3(F10.4,','))
    1002 FORMAT(2(F10.4,','),F10.4,/,');')
  ENDDO
! 軸力の数値出力
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillFore(0,0,0);'
  WRITE(12,*) 'FillBack(65535,26214,13107);'
  WRITE(12,*) 'PenFore(0,0,0);'
  WRITE(12,*) 'PenBack(65535,65535,65535);'
  WRITE(12,*) 'TextSize(',INT(TEXTSIZE),');'
  DO INUM=1,NELEC
    IPOIN=NODEC(INUM,1)
    JPOIN=NODEC(INUM,2)
    FORCE=STRSC(INUM)
    POSX=COOR1(JPOIN,1)-COOR1(IPOIN,1)
    POSY=COOR1(JPOIN,2)-COOR1(IPOIN,2)
    !軸力の数字を書く方向ROTANGLEの計算
    IF(POSX.EQ.0.0)THEN
      IF(POSY.EQ.0.0)ROTANGLE=0.0
      IF(POSY.GT.0.0)ROTANGLE=90.0
      IF(POSY.LT.0.0)ROTANGLE=-90.0
    ELSE
      IF(POSY.EQ.0.0)THEN
        IF(POSX.GT.0.0)ROTANGLE=0.0
        IF(POSX.LT.0.0)ROTANGLE=180.0
      ELSE
        IF(POSX.GT.0.0)ROTANGLE=ATAN(POSY/(+POSX))*180.0/3.1415
        IF(POSX.LT.0.0)ROTANGLE=180.0+ATAN(POSY/(+POSX))*180.0/3.1415
      ENDIF
    ENDIF
    WRITE(12,*) 'TextRotate(',ROTANGLE,');'
    VEC1(1)=(COOR1(IPOIN,1)+COOR1(JPOIN,1))/2.0
    VEC1(2)=(COOR1(IPOIN,2)+COOR1(JPOIN,2))/2.0+OFFSET_Y
    WRITE(12,*) 'TextOrigin(',VEC1(1),',',VEC1(2),');'
    WRITE(12,*) 'BeginText;'
    !軸力の数値の表示桁数の設定
        IF(FORCE.EQ.0.0)THEN
            WRITE(12,1004) 0
        ELSEIF(FORCE.LT.1000.0.AND.FORCE.GT.-99.9999)THEN
            WRITE(12,1005)INT(FORCE)
        ELSE
            WRITE(12,1006)INT(FORCE)
        ENDIF
        1004 FORMAT(1H',I1,1H')
        1005 FORMAT(1H',I3,1H')
        1006 FORMAT(1H',I6,1H')
    WRITE(12,*) 'EndText;'
  END DO
! タイトルの出力
  WRITE(12,*) 'TextRotate(',0,');'
  WRITE(12,*) 'TextOrigin(',COD_MM(1,1)+COD_MM(3,1)/2,',',COD_MM(1,2)+OFFSET_Y-VALMARGIN/2,');'
  WRITE(12,*) 'BeginText;'
  WRITE(12,*) ''' NODAL FORCES -x- '''
  WRITE(12,*) 'EndText;'
  RETURN
  END
! ======================================================================
  SUBROUTINE vscript_stress_frame(ISTRE,COD_MM,COOR0,COOR1,NVARF,NDOFN,NELEF,NODEF,NPOIN,PROPF,STRSF,STRSF_MAX,VALLNG_MAX,TEXTSIZE,TTDIS,VALMARGIN)
! ======================================================================
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  PARAMETER(SMGNFY=0.10)          ! SMGNFY：部材力の長さを調節する係数、0.25の場合にはもっとも長い部材の1/4の長さになる。
  DIMENSION COOR0(NPOIN,3)        ! COOR0(IPOIN,IDIM)：初期形状時の節点座標
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)：節点座標
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )：座標の最小値、COD_MM(2, )：最大値、、COD_MM(3, )：幅、IDIM=1,2（Ｘ、Ｙ座標）
  DIMENSION NODEF(NELEF,3)        ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
  DIMENSION PROPF(NVARF,10)       ! 曲げ要素の材料特性（種類番号、特性項目(1～10)）
  DIMENSION STRSF(NELEF,NDOFN*2)  ! STRSF(IELEF,IEVAB)：曲げ要素IELEFの部材力（部材座標系における部材端荷重)
  DIMENSION TTDIS(NPOIN*NDOFN)    ! TTDIS(ITOTV)：初期形状に対する変位
  DIMENSION VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(3),VEC6(3),RTANG(3)
  DIMENSION OFFSET(3)             ! OFFSET()：図形表示位置のＸ、Ｙ、Ｚ方向移動量
! ISTRE：ISTRE=１：x軸力、２：yせん断力、３：zせん断力、４：xxねじり、５：yy曲げ、６：zz曲げ
! STRSF_MAX：応力の最大値
! VALLNG_MAX：線材の長さの最大値
!
  WRITE(12,*) ''
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''線材-1 x軸力'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''線材-2 yせん断力'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''線材-3 zせん断力'');'
  IF(ISTRE.EQ.4)WRITE(12,*) 'NameClass(''線材-4 xxねじり'');'
  IF(ISTRE.EQ.5)WRITE(12,*) 'NameClass(''線材-5 yy曲げ'');'
  IF(ISTRE.EQ.6)WRITE(12,*) 'NameClass(''線材-6 zz曲げ'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
! 文字の書式設定
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillFore(0,0,0);'
  WRITE(12,*) 'FillBack(65535,26214,13107);'
  WRITE(12,*) 'PenFore(0,0,0);'
  WRITE(12,*) 'PenBack(65535,65535,65535);'
  WRITE(12,*) 'TextSize(',INT(TEXTSIZE),');'
!
  OFFSET(1)=COD_MM(3,1)+VALMARGIN+5              ! OFFSET(1)：図形表示位置のＸ方向移動量
  OFFSET(2)=-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)   ! OFFSET(2)：図形表示位置のＹ方向移動量
  OFFSET(3)=0.0                                  ! OFFSET(3)：図形表示位置のＺ方向移動量
  IEVAB1=ISTRE      ! IEVAB1：部材のｉ端側の応力番号
  IEVAB2=ISTRE+6    ! IEVAB2：部材のｊ端側の応力番号
! タイトルの出力
  WRITE(12,*) 'TextRotate(',0,');'
  WRITE(12,*) 'TextOrigin(',COD_MM(2,1)+OFFSET(1)/2+VALMARGIN,',',COD_MM(1,2)+OFFSET(2)-VALMARGIN/2,');'
  WRITE(12,*) 'BeginText;'
  IF(ISTRE.EQ.1)WRITE(12,*) ''' NODAL FORCES [N] -x- '''
  IF(ISTRE.EQ.2)WRITE(12,*) ''' SHEAR STRESS [N] -y- '''
  IF(ISTRE.EQ.3)WRITE(12,*) ''' SHEAR STRESS [N] -z- '''
  IF(ISTRE.EQ.4)WRITE(12,*) ''' TORSION MOMENT [Nm] -xx- '''
  IF(ISTRE.EQ.5)WRITE(12,*) ''' BENDING MOMENT [Nm] -yy- '''
  IF(ISTRE.EQ.6)WRITE(12,*) ''' BENDING MOMENT [Nm] -zz- '''
  WRITE(12,*) 'EndText;'
!
  DO IELEF=1,NELEF
    INFRM=NODEF(IELEF,3)  ! INFRM：曲げ要素IELEFの部材番号
!   部材の出力
    WRITE(12,*) 'Poly3D('
    IPOIN=NODEF(IELEF,1)
    JPOIN=NODEF(IELEF,2)
    VALLNG=SQRT((COOR1(JPOIN,1)-COOR1(IPOIN,1))**2+(COOR1(JPOIN,2)-COOR1(IPOIN,2))**2+(COOR1(JPOIN,3)-COOR1(IPOIN,3))**2)
    DO IDIM=1,3
      VEC1(IDIM)=(COOR0(JPOIN,IDIM)-COOR0(IPOIN,IDIM))/VALLNG       !VEC1( )：要素ｉｊ方向の単位ベクトル
      VEC2(IDIM)=PROPF(INFRM,IDIM+7)  !VEC2( )：弱軸方向の単位ベクトル
    ENDDO
    CALL VECTPRD(VEC1,VEC2,VEC3)    !VECT3( )：強軸方向の単位ベクトル＝VECT1×VECT2
    DO INODE=1,2
      DO IDIM=1,3
        KPOIN=NODEF(IELEF,INODE)
        ITOTV=(KPOIN-1)*NDOFN+IDIM+3
        RTANG(IDIM)=TTDIS(ITOTV)    ! RTANG(IDIM)：Ｘ、Ｙ、Ｚ軸（IDIM=1～3）周りの変形角TTDIS(ITOTV)
      ENDDO
      IF(ISTRE.EQ.2.OR.ISTRE.EQ.6)THEN
        CALL ROTXYZ(RTANG,VEC2,VEC4)   ! VEC4( )：VEC3( )をＸ、Ｙ、Ｚ軸周りに回転したベクトル
      ELSE
        CALL ROTXYZ(RTANG,VEC3,VEC4)   ! VEC4( )：VEC2( )をＸ、Ｙ、Ｚ軸周りに回転したベクトル
      ENDIF
CONTINUE
      DO IDIM=1,3
        IF(INODE.EQ.1)VEC5(IDIM)=COOR1(IPOIN,IDIM)-STRSF(IELEF,IEVAB1)*VEC4(IDIM)*(VALLNG_MAX/STRSF_MAX)*SMGNFY+OFFSET(IDIM)
        IF(INODE.EQ.2)VEC6(IDIM)=COOR1(JPOIN,IDIM)+STRSF(IELEF,IEVAB2)*VEC4(IDIM)*(VALLNG_MAX/STRSF_MAX)*SMGNFY+OFFSET(IDIM)
        CONTINUE
      ENDDO
    ENDDO
    WRITE(12,1001)COOR1(IPOIN,1)+OFFSET(1),COOR1(IPOIN,2)+OFFSET(2),COOR1(IPOIN,3)+OFFSET(3)
    WRITE(12,1001)(VEC5(IDIM),IDIM=1,3)
    WRITE(12,1001)(VEC6(IDIM),IDIM=1,3)
    WRITE(12,1002)COOR1(JPOIN,1)+OFFSET(1),COOR1(JPOIN,2)+OFFSET(2),COOR1(JPOIN,3)+OFFSET(3)
    1001 FORMAT(3(F10.4,','))
    1002 FORMAT(2(F10.4,','),F10.4,/,');')
!   応力の数値の出力
    POSX=VEC1(1)
    POSY=VEC1(2)
    !応力の数字を書く方向ROTANGLEの計算
    IF(POSX.EQ.0.0)THEN
      IF(POSY.EQ.0.0)ROTANGLE=0.0
      IF(POSY.GT.0.0)ROTANGLE=90.0
      IF(POSY.LT.0.0)ROTANGLE=-90.0
    ELSE
      IF(POSY.EQ.0.0)THEN
        IF(POSX.GT.0.0)ROTANGLE=0.0
        IF(POSX.LT.0.0)ROTANGLE=180.0
      ELSE
        IF(POSX.GT.0.0)ROTANGLE=ATAN(POSY/(+POSX))*180.0/3.1415
        IF(POSX.LT.0.0)ROTANGLE=180.0+ATAN(POSY/(+POSX))*180.0/3.1415
      ENDIF
    ENDIF
ROTANGLE=0.0    !数字の方向を揃える
    WRITE(12,*) 'TextRotate(',ROTANGLE,');'
    DO INODE=1,1
!   DO INODE=1,2
      IPOIN=NODEF(IELEF,INODE)
      IEVAB=ISTRE+(INODE-1)*6
      WRITE(12,*) 'TextOrigin(',COOR1(IPOIN,1)+OFFSET(1),',',COOR1(IPOIN,2)+OFFSET(2),');'
      WRITE(12,*) 'BeginText;'
      WRITE(12,1006)INT(STRSF(IELEF,IEVAB))
      1006 FORMAT(1H',I5,1H')
      WRITE(12,*) 'EndText;'
    ENDDO
  ENDDO
  RETURN
  END
! ======================================================================
  SUBROUTINE ROTXYZ(RTANG,VEC1,VEC4)
! ======================================================================
! VEC1をＸ、Ｙ、Ｚ軸周りにRTANG（1～3)回転して、VEC4を作成するサブルーチン
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION VEC1(3),VEC4(3),RTANG(3)
  DIMENSION RX(3,3),RY(3,3),RZ(3,3)
  RX=0.0
  RY=0.0
  RZ=0.0
  RX(1,1)= 1.0
  RX(2,2)= COS(RTANG(1))
  RX(2,3)=-SIN(RTANG(1))
  RX(3,2)=-RX(2,3)
  RX(3,3)= RX(2,2)
  RY(1,1)= COS(RTANG(2))
  RY(1,3)=-SIN(RTANG(2))
  RY(2,2)= 1.0
  RY(3,1)=-RY(1,3)
  RY(3,3)= RY(1,1)
  RZ(1,1)= COS(RTANG(3))
  RZ(1,2)=-SIN(RTANG(3))
  RZ(2,1)=-RZ(1,2)
  RZ(2,2)= RZ(1,1)
  RZ(3,3)= 1.0
!  CALL MLTPLY_ABC(3,3,3,1,RY,RZ,VEC1,VEC4)     ! VEC4( )←RY(RZ(VEC1))
  CALL MLTPLY(3,3,1,RX,VEC1,VEC4)   ! VEC4( )←RX(VEC4)
  CALL MLTPLY(3,3,1,RY,VEC4,VEC4)   ! VEC4( )←RX(VEC4)
  CALL MLTPLY(3,3,1,RZ,VEC4,VEC4)   ! VEC4( )←RX(VEC4)
  RETURN
  END

! ======================================================================
    SUBROUTINE vscript_class()    ! クラスの出力
! ======================================================================
WRITE(12,*) '{End of Creation Code}'
WRITE(12,*) ''
WRITE(12,*) '{Classes}'
WRITE(12,*) ''
WRITE(12,*) 'NameClass(''一般'');'
WRITE(12,*) 'SetClFillFore(''一般'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''一般'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''一般'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''一般'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''一般'',1);'
WRITE(12,*) 'SetClLS(''一般'',2);'
WRITE(12,*) 'SetClLW(''一般'',1);'
WRITE(12,*) 'SetClUseGraphic(''一般'',FALSE);'
WRITE(12,*) 'NameClass(''寸法'');'
WRITE(12,*) 'SetClFillFore(''寸法'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''寸法'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''寸法'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''寸法'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''寸法'',1);'
WRITE(12,*) 'SetClLS(''寸法'',2);'
WRITE(12,*) 'SetClLW(''寸法'',1);'
WRITE(12,*) 'SetClUseGraphic(''寸法'',FALSE);'
WRITE(12,*) 'NameClass(''線材-0 要素'');'
WRITE(12,*) 'SetClFillFore(''線材-0 要素'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-0 要素'',0,0,0);'
WRITE(12,*) 'SetClPenFore(''線材-0 要素'',56797,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-0 要素'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-0 要素'',0);'
WRITE(12,*) 'SetClLS(''線材-0 要素'',2);'
WRITE(12,*) 'SetClLW(''線材-0 要素'',6);'
WRITE(12,*) 'SetClUseGraphic(''線材-0 要素'',TRUE);'
WRITE(12,*) 'NameClass(''膜_0 要素'');'
WRITE(12,*) 'SetClFillFore(''膜_0 要素'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''膜_0 要素'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''膜_0 要素'',53520,41058,29490);'
WRITE(12,*) 'SetClPenBack(''膜_0 要素'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''膜_0 要素'',1);'
WRITE(12,*) 'SetClLS(''膜_0 要素'',2);'
WRITE(12,*) 'SetClLW(''膜_0 要素'',3);'
WRITE(12,*) 'SetClUseGraphic(''膜_0 要素'',TRUE);'
WRITE(12,*) 'NameClass(''膜_1 張力（x方向）'');'
WRITE(12,*) 'SetClFillFore(''膜_1 張力（x方向）'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''膜_1 張力（x方向）'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''膜_1 張力（x方向）'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''膜_1 張力（x方向）'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''膜_1 張力（x方向）'',1);'
WRITE(12,*) 'SetClLS(''膜_1 張力（x方向）'',2);'
WRITE(12,*) 'SetClLW(''膜_1 張力（x方向）'',3);'
WRITE(12,*) 'SetClUseGraphic(''膜_1 張力（x方向）'',FALSE);'
WRITE(12,*) 'NameClass(''線材-2 yせん断力'');'
WRITE(12,*) 'SetClFillFore(''線材-2 yせん断力'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-2 yせん断力'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-2 yせん断力'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-2 yせん断力'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-2 yせん断力'',0);'
WRITE(12,*) 'SetClLS(''線材-2 yせん断力'',2);'
WRITE(12,*) 'SetClLW(''線材-2 yせん断力'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-2 yせん断力'',TRUE);'
WRITE(12,*) 'NameClass(''区分-膜_1 張力x'');'
WRITE(12,*) 'SetClFillFore(''区分-膜_1 張力x'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''区分-膜_1 張力x'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''区分-膜_1 張力x'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''区分-膜_1 張力x'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''区分-膜_1 張力x'',1);'
WRITE(12,*) 'SetClLS(''区分-膜_1 張力x'',2);'
WRITE(12,*) 'SetClLW(''区分-膜_1 張力x'',1);'
WRITE(12,*) 'SetClUseGraphic(''区分-膜_1 張力x'',FALSE);'
WRITE(12,*) 'NameClass(''区分-膜_2 張力y'');'
WRITE(12,*) 'SetClFillFore(''区分-膜_2 張力y'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''区分-膜_2 張力y'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''区分-膜_2 張力y'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''区分-膜_2 張力y'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''区分-膜_2 張力y'',1);'
WRITE(12,*) 'SetClLS(''区分-膜_2 張力y'',2);'
WRITE(12,*) 'SetClLW(''区分-膜_2 張力y'',1);'
WRITE(12,*) 'SetClUseGraphic(''区分-膜_2 張力y'',FALSE);'
WRITE(12,*) 'NameClass(''区分-膜_3 せん断力xy'');'
WRITE(12,*) 'SetClFillFore(''区分-膜_3 せん断力xy'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''区分-膜_3 せん断力xy'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''区分-膜_3 せん断力xy'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''区分-膜_3 せん断力xy'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''区分-膜_3 せん断力xy'',1);'
WRITE(12,*) 'SetClLS(''区分-膜_3 せん断力xy'',2);'
WRITE(12,*) 'SetClLW(''区分-膜_3 せん断力xy'',1);'
WRITE(12,*) 'SetClUseGraphic(''区分-膜_3 せん断力xy'',FALSE);'
WRITE(12,*) 'NameClass(''主軸方向ij'');'
WRITE(12,*) 'SetClFillFore(''主軸方向ij'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''主軸方向ij'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''主軸方向ij'',22102,11308,1285);'
WRITE(12,*) 'SetClPenBack(''主軸方向ij'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''主軸方向ij'',0);'
WRITE(12,*) 'SetClLS(''主軸方向ij'',-2);'
WRITE(12,*) 'SetClLW(''主軸方向ij'',3);'
WRITE(12,*) 'SetClUseGraphic(''主軸方向ij'',TRUE);'
WRITE(12,*) 'NameClass(''膜_3 せん断力xy'');'
WRITE(12,*) 'SetClFillFore(''膜_3 せん断力xy'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''膜_3 せん断力xy'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''膜_3 せん断力xy'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''膜_3 せん断力xy'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''膜_3 せん断力xy'',1);'
WRITE(12,*) 'SetClLS(''膜_3 せん断力xy'',2);'
WRITE(12,*) 'SetClLW(''膜_3 せん断力xy'',3);'
WRITE(12,*) 'SetClUseGraphic(''膜_3 せん断力xy'',FALSE);'
WRITE(12,*) 'NameClass(''線材-3 zせん断力'');'
WRITE(12,*) 'SetClFillFore(''線材-3 zせん断力'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-3 zせん断力'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-3 zせん断力'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-3 zせん断力'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-3 zせん断力'',0);'
WRITE(12,*) 'SetClLS(''線材-3 zせん断力'',2);'
WRITE(12,*) 'SetClLW(''線材-3 zせん断力'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-3 zせん断力'',TRUE);'
WRITE(12,*) 'NameClass(''膜_2 張力（y方向）'');'
WRITE(12,*) 'SetClFillFore(''膜_2 張力（y方向）'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''膜_2 張力（y方向）'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''膜_2 張力（y方向）'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''膜_2 張力（y方向）'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''膜_2 張力（y方向）'',1);'
WRITE(12,*) 'SetClLS(''膜_2 張力（y方向）'',2);'
WRITE(12,*) 'SetClLW(''膜_2 張力（y方向）'',3);'
WRITE(12,*) 'SetClUseGraphic(''膜_2 張力（y方向）'',FALSE);'
WRITE(12,*) 'NameClass(''線材-1 x軸力'');'
WRITE(12,*) 'SetClFillFore(''線材-1 x軸力'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-1 x軸力'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-1 x軸力'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-1 x軸力'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-1 x軸力'',0);'
WRITE(12,*) 'SetClLS(''線材-1 x軸力'',2);'
WRITE(12,*) 'SetClLW(''線材-1 x軸力'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-1 x軸力'',TRUE);'
WRITE(12,*) 'NameClass(''膜_4 張力矢印'');'
WRITE(12,*) 'SetClFillFore(''膜_4 張力矢印'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''膜_4 張力矢印'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''膜_4 張力矢印'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''膜_4 張力矢印'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''膜_4 張力矢印'',1);'
WRITE(12,*) 'SetClLS(''膜_4 張力矢印'',2);'
WRITE(12,*) 'SetClLW(''膜_4 張力矢印'',3);'
WRITE(12,*) 'SetClUseGraphic(''膜_4 張力矢印'',TRUE);'
WRITE(12,*) 'NameClass(''番号-0 節点'');'
WRITE(12,*) 'SetClFillFore(''番号-0 節点'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''番号-0 節点'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''番号-0 節点'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''番号-0 節点'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''番号-0 節点'',0);'
WRITE(12,*) 'SetClLS(''番号-0 節点'',2);'
WRITE(12,*) 'SetClLW(''番号-0 節点'',1);'
WRITE(12,*) 'SetClUseGraphic(''番号-0 節点'',TRUE);'
WRITE(12,*) 'NameClass(''番号-1 膜要素'');'
WRITE(12,*) 'SetClFillFore(''番号-1 膜要素'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''番号-1 膜要素'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''番号-1 膜要素'',0,0,54272);'
WRITE(12,*) 'SetClPenBack(''番号-1 膜要素'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''番号-1 膜要素'',0);'
WRITE(12,*) 'SetClLS(''番号-1 膜要素'',2);'
WRITE(12,*) 'SetClLW(''番号-1 膜要素'',1);'
WRITE(12,*) 'SetClUseGraphic(''番号-1 膜要素'',TRUE);'
WRITE(12,*) 'NameClass(''番号-2 線材'');'
WRITE(12,*) 'SetClFillFore(''番号-2 線材'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''番号-2 線材'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''番号-2 線材'',56797,0,0);'
WRITE(12,*) 'SetClPenBack(''番号-2 線材'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''番号-2 線材'',0);'
WRITE(12,*) 'SetClLS(''番号-2 線材'',2);'
WRITE(12,*) 'SetClLW(''番号-2 線材'',1);'
WRITE(12,*) 'SetClUseGraphic(''番号-2 線材'',TRUE);'
WRITE(12,*) 'NameClass(''線材-4 xxねじり'');'
WRITE(12,*) 'SetClFillFore(''線材-4 xxねじり'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-4 xxねじり'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-4 xxねじり'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-4 xxねじり'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-4 xxねじり'',0);'
WRITE(12,*) 'SetClLS(''線材-4 xxねじり'',2);'
WRITE(12,*) 'SetClLW(''線材-4 xxねじり'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-4 xxねじり'',TRUE);'
WRITE(12,*) 'NameClass(''線材-5 yy曲げ'');'
WRITE(12,*) 'SetClFillFore(''線材-5 yy曲げ'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-5 yy曲げ'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-5 yy曲げ'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-5 yy曲げ'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-5 yy曲げ'',0);'
WRITE(12,*) 'SetClLS(''線材-5 yy曲げ'',2);'
WRITE(12,*) 'SetClLW(''線材-5 yy曲げ'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-5 yy曲げ'',TRUE);'
WRITE(12,*) 'NameClass(''線材-6 zz曲げ'');'
WRITE(12,*) 'SetClFillFore(''線材-6 zz曲げ'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''線材-6 zz曲げ'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''線材-6 zz曲げ'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''線材-6 zz曲げ'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''線材-6 zz曲げ'',0);'
WRITE(12,*) 'SetClLS(''線材-6 zz曲げ'',2);'
WRITE(12,*) 'SetClLW(''線材-6 zz曲げ'',3);'
WRITE(12,*) 'SetClUseGraphic(''線材-6 zz曲げ'',TRUE);'
WRITE(12,*) ''
WRITE(12,*) '{End of Class Entries}'
WRITE(12,*) ''
WRITE(12,*) 'END;'
WRITE(12,*) ''
WRITE(12,*) 'Run(LoadFile);'
RETURN
END
! ======================================================================
!    SUBROUTINE vscript_
! ======================================================================
!  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
!  RETURN
!  END
END SUBROUTINE WRITE_VSCRIPT
! ... (他のサブルーチンや関数があれば、それらのEND SUBROUTINE/FUNCTIONの後に)
END MODULE WRITE_VSCRIPT_MOD