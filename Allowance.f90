! =====================================================================
    SUBROUTINE ALLOWANCE(COORD,NPOIN,NVARF,NDOFN,NELEF,NODEF,PROPF,STRSF,VALK)
! =====================================================================
! 鉄骨部材（鋼管）の許容応力度検定
    IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    DIMENSION COORD(NPOIN,3)    ! COORD(IPOIN,IDIR) 節点座標
    DIMENSION NODEF(NELEF,3)    ! 曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
    DIMENSION STRSF(NELEF,NDOFN*2)                                    ! STRSF(IELEF,IEVAB)：曲げ要素の応力（部材座標系における部材端荷重)
    DIMENSION VALK(NELEF,2)     ! VALK(IELEF, ) 材端1,2の接合条件（＝1.0E+20：剛接、＝0.0：ピン接、その他：半剛接）
    DIMENSION PROPF(NVARF,10)   ! 曲げ要素の材料特性（種類番号、特性項目(1～10)）
!
    DIMENSION PRFRM(NVARF,8)    ! PRFRM(IVARF,INUM) 曲げ部材の材料諸元
    DOUBLEPRECISION J,JT,LBY,LBZ,LKY,LKZ,LLMD,LMDY,LMDZ,LT,MJ,MY,MZ,QY,QZ,N

    OPEN(14,FILE='allow.txt')   ! 検定結果の出力先ファイル

    PI=4.0*ATAN(1.0)    ! PI π

DO IE=1,NELEF
    N1=NODEF(IE,1)  ! N1 部材１端の節点番号
    N2=NODEF(IE,2)  ! N2 部材２端の節点番号
    IV=NODEF(IE,3)  ! IV 部材リスト番号
    EA=PROPF(IV,1)  ! EA 軸剛性[N]
    EIy=PROPF(IV,3) ! EIy ｙ軸周り曲げ剛性[Nm2]
    PRFRM(IV,1)=1.0     ! PRFRM(IV,1) 部材タイプ、=1円形鋼管、=2角形鋼管、=3Ｈ形鋼、=4平鋼、=5丸鋼
    PRFRM(IV,2)=139.8   ! PRFRM(IV,2) WD1 外法1 mm
    PRFRM(IV,3)=0.0     ! PRFRM(IV,3) WD2 外法2 mm
    PRFRM(IV,4)=4.5     ! PRFRM(IV,4) TH1 厚さ1 mm
    PRFRM(IV,5)=0.0     ! PRFRM(IV,5) TH2 厚さ2 mm
    PRFRM(IV,6)=7.8     ! PRFRM(IV,6) DNS 密度  g/cm3, t/m3
    PRFRM(IV,7)=235.4   ! PRFRM(IV,7) SST 基準強度 MPa,N/mm2
    PRFRM(IV,8)=2.06E+5 ! PRFRM(IV,4) YNG ヤング率 MPa,N/mm2

    ITFRM=INT(PRFRM(IV,1))   ! ITFRM 部材タイプ、=1円形鋼管、=2角形鋼管、=3Ｈ形鋼、=4平鋼、=5丸鋼
    SST=PRFRM(IV,7)     ! SST 基準強度
    YNG=PRFRM(IV,8)     ! YNG ヤング率
    POI=0.3D0           ! POI ポアソン比
    LLMD=SQRT(PI**2*YNG/0.6D0/SST)  ! 限界細長比Λ
    IF(ITFRM.EQ.1)THEN  ! 円形鋼管
      RAD=SQRT(2.0*EIy/EA)*1.0D+3   ! RAD 半径 mm
      THK=EA*9.81/YNG/2.0/PI/RAD    ! THK 肉厚 mm
      DMT=(RAD+THK/2.0)*2.0         ! DMT 外径 mm
      DMT1=(RAD-THK/2.0)*2.0        ! DMT1 内径 mm
    ENDIF

    IEBND=3 ! IEBND 材端の拘束条件 =1：ピン＋ピン、=2：ピン＋剛、=3：剛＋剛
    IF(VALK(IE,1).LE.1.0)IEBND=IEBND-1
    IF(VALK(IE,2).LE.1.0)IEBND=IEBND-1

!デバッグ用のサンプルデータ
!IEBND=3
!LT=5000.0D0
!STRSF(1,1)= -4703.669725
!STRSF(1,2)=  48.92966361
!STRSF(1,3)=  360.8562691
!STRSF(1,4)=  60.95820591
!STRSF(1,5)=  706.9317023
!STRSF(1,6)=  -88.37920489
!STRSF(1,7)=  4703.669725
!STRSF(1,8)=  48.92966361
!STRSF(1,9)=  360.8562691
!STRSF(1,10)=  60.95820591
!STRSF(1,11)=  706.9317023
!STRSF(1,12)=  -88.37920489

    IF(IEBND.EQ.1)THEN   ! FLK 座屈長さ（圧縮）の係数
      FLK=1.0D0
    ELSEIF(IEBND.EQ.2)THEN
      FLK=0.65D0
    ELSEIF(IEBND.EQ.3)THEN
      FLK=0.8D0
    ENDIF

    LOAD_TYPE=1 ! LOAD_TYPE 荷重の種類、=0 長期荷重、=1 短期荷重
    IF(LOAD_TYPE.EQ.0)THEN  ! LOAD_TYPE 荷重の種類、=0 長期荷重、=1 短期荷重
      SAFE_FACT=1.5     ! SAFE_FACT 長期安全率1.5
    ELSE
      SAFE_FACT=1.0     ! SAFE_FACT 短期安全率1.0
    ENDIF

    LT=0.0  ! LT 変形後部材長
    DO IDIM=1,3
      LT=LT+(COORD(N2,IDIM)-COORD(N1,IDIM))**2
    ENDDO
    LT=SQRT(LT)*1.0D+3

    AN=PI*(DMT**2-DMT1**2)/4.0          ! AN 圧縮用断面積 mm2
    AT=AN                               ! AT 引張用断面積 mm2
    ASZ=AN/1910D0*1013D0    ! 算定式不明につき暫定 ! ASZ せん断断面積Ｚ mm2
    ASY=ASZ                 ! 算定式不明につき暫定 ! ASY せん断断面積Ｙ mm2
    J=PI*RAD**3*THK/(1.0D0+POI)         ! J ねじり剛性係数
    JT=J/RAD                            ! JT ねじり抵抗係数
    ZZT=PI*(DMT**4-DMT1**4)/32.0/DMT    ! ZZT 断面係数Ｚ上部 mm3
    ZZB=ZZT                             ! ZZB 断面係数Ｚ下部 mm3
    ZYL=ZZT                             ! ZYL 断面係数Ｙ左部 mm3
    ZYR=ZYL                             ! ZYR 断面係数Ｙ右部 mm3
    SIZ=SQRT(DMT**2+DMT1**2)/4.0        ! SIZ 断面2次半径Ｚ mm
    SIY=SIZ                             ! SIY 断面2次半径Ｙ mm
    LKZ=LT*FLK                          ! LKZ 座屈長さＺ mm
    LKY=LKZ                             ! LKY 座屈長さＹ mm
    !LBZ=?mm                            ! LBZ 曲げ座屈計算用長さＺ mm
    !LBY=?mm                            ! LBY 曲げ座屈計算用長さＹ mm

IF(IV.EQ.1)THEN ! 円弧アーチの部材長さ
   LKZ=2500.0D0
   LKY=2500.0D0
ENDIF

    LMDZ=LKZ/SIZ                        ! LMDZ 細長比Ｚ λz
    LMDY=LKY/SIZ                        ! LMDY 細長比Ｙ λy
    RATIO=MIN(LMDZ,LMDY)/LLMD           ! λ/Λ
    IF(RATIO.LE.1.0)THEN                ! λ≦Λ
      VALNY=3.0/2.0+2.0/3.0*(RATIO)**2  ! VALNY 安全率 ν
      FC=(1.0-0.4*RATIO**2)*SST/VALNY*(1.5/SAFE_FACT)   ! FC 許容圧縮応力度 MPa
    ELSE                                ! λ＞Λ
      VALNY=2.17
      FC=0.277*SST/RATIO**2
    ENDIF
    FT=SST/SAFE_FACT                    ! FT 許容引張応力度 MPa
    FBZ=SST                             ! FBZ 許容曲げ応力度Ｚ MPa
    FBY=SST                             ! FBY 許容曲げ応力度Ｙ MPa
    FSZ=SST/SQRT(3.0)/SAFE_FACT         ! FSZ 許容せん断応力度Ｚ MPa
    FSY=SST/SQRT(3.0)/SAFE_FACT         ! FSY=許容せん断応力度Ｙ MPa

    N=-STRSF(IE,1)*9.81*1.0D-3          ! N 軸力 kN
    SND=N*1.0D+3/AN                     ! SND 発生軸応力度 MPa
    IF(SND.LE.0.0)THEN
      RND=-SND/FC                       ! RND 軸力判定値（圧縮)
    ELSE
      RND=SND/FT                        ! RND 軸力判定値（引張)
    ENDIF
    MZ=MAX(ABS(STRSF(IE,5)),ABS(STRSF(IE,11)))*9.81*1.0D-3      ! MZ 曲げモーメントＺ kNm
    SBZ=MZ*1.0D+6/MIN(ZZT,ZZB)          ! SBZ 曲げ応力度Ｚ MPa
    RBZ=SBZ/FBZ                         ! RBZ 曲げ応力度Ｚ判定値
    MY=MAX(ABS(STRSF(IE,6)),ABS(STRSF(IE,12)))*9.81*1.0D-3      ! MY 曲げモーメントＹ kNm
    SBY=MY*1.0D+6/MIN(ZYL,ZYR)          ! SBY 曲げ応力度Ｙ MPa
    RBY=SBY/FBY                         ! RBY 曲げ応力度Ｙ判定値
    QZ=MAX(ABS(STRSF(IE,3)),ABS(STRSF(IE,9)))*9.81*1.0D-3       ! QZ せん断力Ｚ kN
    TAZ=QZ*1.0D+3/ASZ                   ! TAZ せん断応力度Ｚ MPa
    RTZ=TAZ/FSZ                         ! RTZ せん断応力度Ｚ判定値
    QY=MAX(ABS(STRSF(IE,2)),ABS(STRSF(IE,8)))*9.81*1.0D-3       ! QY せん断力Ｙ kN
    TAY=QY*1.0D+3/ASY                   ! TAY せん断応力度Ｙ MPa
    RTY=TAY/FSY                         ! RTY せん断応力度Ｙ判定値
    MJ=MAX(ABS(STRSF(IE,4)),ABS(STRSF(IE,10)))*9.81*1.0D-3      ! MJ ねじりモーメント kNm
    TAJ=MJ*1.0D+6/JT                    ! TAJ 捩りせん断応力度 MPa
    RTJ=TAJ/JT                          ! RTJ ねじりせん断応力度判定値    
    RCB=SQRT((ABS(SND)+SQRT(SBZ**2+SBY**2))**2+3.0*(SQRT(TAZ**2+TAY**2)+TAJ)**2)/FT
!   RCB=SQRT((SND+SQRT(SBZ**2+SBY**2))**2)/FT   ! RCB 組合せ応力度

! 検定結果のファイル出力

    IF(IE.EQ.1)CALL TEXT    ! 各変数の説明をファイルに出力

! 部材名称=CONCATENATE("P-",E55,"φ×",E56) ! 部材名称

    WRITE(14,401) IE
    WRITE(14,*)'鉄骨部材（鋼管）の検定'
    WRITE(14,*)''
    WRITE(14,402) DMT,THK
    WRITE(14,403) DMT
    WRITE(14,404) THK
    WRITE(14,*)''
    WRITE(14,*)'材質      STK400'
    WRITE(14,405) SST
    WRITE(14,406) YNG
    WRITE(14,407) LLMD
401 FORMAT(///,'=== IELEF = ',I4,/)
402 FORMAT('部材名称',5X,'P-',F6.1,'φ×',F4.1)
403 FORMAT('            外径        ',F10.1,' mm')
404 FORMAT('            肉厚        ',F10.1,' mm')
405 FORMAT('            基準強度    ',F10.1,' MPa')
406 FORMAT('            ヤング率    ',E10.2,' MPa')
407 FORMAT('            限界細長比Λ',F10.1,'    ')
    WRITE(14,*)''
    WRITE(14,*)'材端の拘束条件'
!   WRITE(14,'('            材端拘束条件Z =',I2,10X,'材端拘束条件Y =',I2)') 
    WRITE(14,411) IEBND
411 FORMAT('            材端拘束条件  =',I2)
    WRITE(14,*)''
    WRITE(14,*)'荷重の種類'
    IF(LOAD_TYPE.EQ.0) WRITE(14,412) SAFE_FACT
    IF(LOAD_TYPE.EQ.1) WRITE(14,413) SAFE_FACT
412 FORMAT('          長期荷重',/,'          安全率',F15.1)
413 FORMAT('          短期荷重',/,'          安全率',F15.1)
    WRITE(14,*)''
    WRITE(14,421) LT
    WRITE(14,422) AN,AT
    WRITE(14,423) ASZ,ASY
    WRITE(14,424) J,JT
    WRITE(14,425) ZZT,ZZB
    WRITE(14,426) ZYL,ZYR
    WRITE(14,427) SIZ,SIY
    WRITE(14,428) LKZ,LKY
!   WRITE(14,429) LBZ,LNY
421 FORMAT('LT      = ',F15.1,' [mm ]')
422 FORMAT('AN      = ',F15.1,' [mm2]',5X,'AT      =',F15.1,' [mm2]')
423 FORMAT('ASZ     = ',F15.1,' [mm2]',5X,'ASY     =',F15.1,' [mm2]')
424 FORMAT('J       = ',F15.1,' [mm4]',5X,'JT      =',F15.1,' [mm3]')
425 FORMAT('ZZT     = ',F15.1,' [mm3]',5X,'ZZB     =',F15.1,' [mm3]')
426 FORMAT('ZYL     = ',F15.1,' [mm3]',5X,'ZYR     =',F15.1,' [mm3]')
427 FORMAT('SIZ     = ',F15.1,' [mm ]',5X,'SIY     =',F15.1,' [mm ]')
428 FORMAT('LKZ     = ',F15.1,' [mm ]',5X,'LKY     =',F15.1,' [mm ]')
429 FORMAT('LBZ     = ',F15.1,' [mm ]',5X,'LBY     =',F15.1,' [mm ]')
    WRITE(14,*)''
    WRITE(14,431) LMDZ,LMDY
    WRITE(14,432) FC,FT
    WRITE(14,433) FBZ,FBY
    WRITE(14,434) FSZ,FSY
431 FORMAT('LMDZ    = LKZ/SIZ = ',F5.1,11X,'LMDY    = LKY/SIY = ',F5.1)
432 FORMAT('FC      = ',F15.1,' [MPa]',5X,'FT      = ',F15.1,' [MPa]')
433 FORMAT('FBZ     = ',F15.1,' [MPa]',5X,'FBY     = ',F15.1,' [MPa]')
434 FORMAT('FSZ     = ',F15.1,' [MPa]',5X,'FSY     = ',F15.1,' [MPa]')
    WRITE(14,*)''
    IF(N.LT.0.0)WRITE(14,440) N,SND
    IF(N.GE.0.0)WRITE(14,441) N,SND
    WRITE(14,442) MZ,SBZ
    WRITE(14,443) MY,SBY
    WRITE(14,444) QZ,TAZ 
    WRITE(14,445) QY,TAY
    WRITE(14,446) MJ,TAJ
440 FORMAT('N       = ',F15.2,' [kN ]',5X,'ST, SC  = N/AN     = ',F15.1,' [MPa]')
441 FORMAT('N       = ',F15.2,' [kN ]',5X,'ST, SC  = N/AT     = ',F15.1,' [MPa]')
442 FORMAT('MZ      = ',F15.2,' [kNm]',5X,'SBZ     = MZ/ZZ    = ',F15.1,' [MPa]')
443 FORMAT('MY      = ',F15.2,' [kNm]',5X,'SBY     = MY/ZY    = ',F15.1,' [MPa]')
444 FORMAT('QZ      = ',F15.2,' [kN ]',5X,'TAZ     = QZ/ASZ   = ',F15.1,' [MPa]')
445 FORMAT('QY      = ',F15.2,' [kN ]',5X,'TAY     = QY/ASY   = ',F15.1,' [MPa]')
446 FORMAT('MJ      = ',F15.2,' [kNm]',5X,'TAJ     = MJ/JT    = ',F15.1,' [MPa]')
    WRITE(14,*)''
    IF(RND.LT.0.0)WRITE(14,450)  RND
    IF(RND.GE.0.0)WRITE(14,451)  RND
    WRITE(14,452)  RBZ
    WRITE(14,453)  RBY
    WRITE(14,454)  RTZ
    WRITE(14,455)  RTY
    WRITE(14,456)  RCB
450 FORMAT('ST/FC   = ',F15.4)
451 FORMAT('ST/FT   = ',F15.4)
452 FORMAT('SBZ/FBZ = ',F15.4)
453 FORMAT('SBY/FBY = ',F15.4)
454 FORMAT('TAZ/FSZ = ',F15.4)
455 FORMAT('TAY/FSY = ',F15.4)
456 FORMAT('SQRT((SND+SQRT(SBZ**2+SBY**2))**2+3.0*(SQRT(TAZ**2+TAY**2)+TAJ)**2)/FT=',F15.3)
    WRITE(14,*)''
    IF(RND.LE.1.0.AND.RBZ.LE.1.0.AND.RBY.LE.1.0.AND.RTZ.LE.1.0.AND.RTY.LE.1.0.AND.RCB.LE.1.0)THEN
      WRITE(14,457) DMT,THK
    ELSE
      WRITE(14,458) DMT,THK
    ENDIF
457 FORMAT('USE :  ','P-',F6.1,'φ×',F4.1,' -----> OK')
458 FORMAT('USE :  ','P-',F6.1,'φ×',F4.1,' -----> NG')

ENDDO ! LOOP OVER NELEF

RETURN
END

    SUBROUTINE TEXT     ! 変数の説明を出力する
    WRITE(14,*)'鉄骨部材（鋼管）の検定'
    WRITE(14,*)''
    WRITE(14,*)'部材名称'
    WRITE(14,*)'            外径'
    WRITE(14,*)'            肉厚'
    WRITE(14,*)''
    WRITE(14,*)'材質      STK400'
    WRITE(14,*)'            基準強度'
    WRITE(14,*)'            ヤング率'
    WRITE(14,*)'            ポアソン比'
    WRITE(14,*)'          限界細長比Λ'
    WRITE(14,*)''
    WRITE(14,*)'材端の拘束条件'
    WRITE(14,*)'          材端拘束条件Z                  材端拘束条件Z'
    WRITE(14,*)'             =1：ピン＋ピン、=2：剛＋剛、=3：ピン＋剛、=4：剛＋剛ローラー、=5：ピン＋剛ローラー'
    WRITE(14,*)'          座屈長さ（圧縮）の係数Z        座屈長さ（圧縮）の係数Y'
    WRITE(14,*)''
    WRITE(14,*)'荷重の種類'
    WRITE(14,*)'          短期or長期          1       =0 長期荷重、 =1 短期荷重'
    WRITE(14,*)'          安全率           1.0'
    WRITE(14,*)''
    WRITE(14,*)'L0      = 変形前部材長         ＬＴ    = 変形後部材長'
    WRITE(14,*)'AN      = 圧縮用断面積         AT      = 引張用断面積'
    WRITE(14,*)'ASZ     = せん断断面積Z        ASY     = せん断断面積Ｙ'
    WRITE(14,*)'J       = ねじり剛性係数       JT      = ねじり抵抗係数'
    WRITE(14,*)'ZZT     = 断面係数Z上部        ZXB     = 断面係数Ｚ下部'
    WRITE(14,*)'ZYL     = 断面係数Ｙ左部       ZYR     = 断面係数Ｙ右部'
    WRITE(14,*)'SIZ     = 断面2次半径Ｚ        SIY     = 断面2次半径Ｙ'
    WRITE(14,*)'LKZ     = 座屈長さＺ           LKY     = 座屈長さＹ'
    WRITE(14,*)'LBZ     = 曲げ座屈計算用長さＺ LBY     = 曲げ座屈計算用長さＹ'
    WRITE(14,*)''
    WRITE(14,*)'LMDZ    = LKZ/SIZ =細長比Ｚ    LMDY    = LKY/SIY =細長比Ｙ'
    WRITE(14,*)'FC      = 許容圧縮応力度       FT      = 許容引張応力度'
    WRITE(14,*)'FBZ     = 許容曲げ応力度Ｚ     FBY     = 許容曲げ応力度Ｙ'
    WRITE(14,*)'FSZ     = 許容せん断応力度Ｚ   FSY     = 許容せん断応力度Ｙ'
    WRITE(14,*)''
    WRITE(14,*)'N       = 軸力                 ST, SC  = N/AN     = 発生軸応力度'
    WRITE(14,*)'MZ      = 曲げモーメントＺ     SBZ     = MZ/ZZ    = 発生曲げ応力度Ｚ'
    WRITE(14,*)'MY      = 曲げモーメントＹ     SBY     = MY/ZY    = 発生曲げ応力度Ｙ'
    WRITE(14,*)'QXY     = せん断力Ｚ           TAZ     = QZ/ASZ  = 発生せん断応力度Ｚ'
    WRITE(14,*)'QY      = せん断力Ｙ           TAY     = QY/ASY  = 発生せん断応力度Ｙ'
    WRITE(14,*)'MJ      = ねじりモーメント     TAJ     = MJ/JT    = 発生ねじりせん断応力度'
    WRITE(14,*)''
    WRITE(14,*)'算定式'
    WRITE(14,*)'算定式＝判定値＜１．０'
    WRITE(14,*)''
    WRITE(14,*)'USE：使用部材　---->　判定'
    RETURN
    END