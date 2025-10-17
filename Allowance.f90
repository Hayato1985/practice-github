! =====================================================================
    SUBROUTINE ALLOWANCE(COORD,NPOIN,NVARF,NDOFN,NELEF,NODEF,PROPF,STRSF,VALK)
! =====================================================================
! �S�����ށi�|�ǁj�̋��e���͓x����
    IMPLICIT DOUBLEPRECISION(A-H,O-Z)
    DIMENSION COORD(NPOIN,3)    ! COORD(IPOIN,IDIR) �ߓ_���W
    DIMENSION NODEF(NELEF,3)    ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
    DIMENSION STRSF(NELEF,NDOFN*2)                                    ! STRSF(IELEF,IEVAB)�F�Ȃ��v�f�̉��́i���ލ��W�n�ɂ����镔�ޒ[�׏d)
    DIMENSION VALK(NELEF,2)     ! VALK(IELEF, ) �ޒ[1,2�̐ڍ������i��1.0E+20�F���ځA��0.0�F�s���ځA���̑��F�����ځj
    DIMENSION PROPF(NVARF,10)   ! �Ȃ��v�f�̍ޗ������i��ޔԍ��A��������(1�`10)�j
!
    DIMENSION PRFRM(NVARF,8)    ! PRFRM(IVARF,INUM) �Ȃ����ނ̍ޗ�����
    DOUBLEPRECISION J,JT,LBY,LBZ,LKY,LKZ,LLMD,LMDY,LMDZ,LT,MJ,MY,MZ,QY,QZ,N

    OPEN(14,FILE='allow.txt')   ! ���茋�ʂ̏o�͐�t�@�C��

    PI=4.0*ATAN(1.0)    ! PI ��

DO IE=1,NELEF
    N1=NODEF(IE,1)  ! N1 ���ނP�[�̐ߓ_�ԍ�
    N2=NODEF(IE,2)  ! N2 ���ނQ�[�̐ߓ_�ԍ�
    IV=NODEF(IE,3)  ! IV ���ރ��X�g�ԍ�
    EA=PROPF(IV,1)  ! EA ������[N]
    EIy=PROPF(IV,3) ! EIy ��������Ȃ�����[Nm2]
    PRFRM(IV,1)=1.0     ! PRFRM(IV,1) ���ރ^�C�v�A=1�~�`�|�ǁA=2�p�`�|�ǁA=3�g�`�|�A=4���|�A=5�ۍ|
    PRFRM(IV,2)=139.8   ! PRFRM(IV,2) WD1 �O�@1 mm
    PRFRM(IV,3)=0.0     ! PRFRM(IV,3) WD2 �O�@2 mm
    PRFRM(IV,4)=4.5     ! PRFRM(IV,4) TH1 ����1 mm
    PRFRM(IV,5)=0.0     ! PRFRM(IV,5) TH2 ����2 mm
    PRFRM(IV,6)=7.8     ! PRFRM(IV,6) DNS ���x  g/cm3, t/m3
    PRFRM(IV,7)=235.4   ! PRFRM(IV,7) SST ����x MPa,N/mm2
    PRFRM(IV,8)=2.06E+5 ! PRFRM(IV,4) YNG �����O�� MPa,N/mm2

    ITFRM=INT(PRFRM(IV,1))   ! ITFRM ���ރ^�C�v�A=1�~�`�|�ǁA=2�p�`�|�ǁA=3�g�`�|�A=4���|�A=5�ۍ|
    SST=PRFRM(IV,7)     ! SST ����x
    YNG=PRFRM(IV,8)     ! YNG �����O��
    POI=0.3D0           ! POI �|�A�\����
    LLMD=SQRT(PI**2*YNG/0.6D0/SST)  ! ���E�ג��䃩
    IF(ITFRM.EQ.1)THEN  ! �~�`�|��
      RAD=SQRT(2.0*EIy/EA)*1.0D+3   ! RAD ���a mm
      THK=EA*9.81/YNG/2.0/PI/RAD    ! THK ���� mm
      DMT=(RAD+THK/2.0)*2.0         ! DMT �O�a mm
      DMT1=(RAD-THK/2.0)*2.0        ! DMT1 ���a mm
    ENDIF

    IEBND=3 ! IEBND �ޒ[�̍S������ =1�F�s���{�s���A=2�F�s���{���A=3�F���{��
    IF(VALK(IE,1).LE.1.0)IEBND=IEBND-1
    IF(VALK(IE,2).LE.1.0)IEBND=IEBND-1

!�f�o�b�O�p�̃T���v���f�[�^
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

    IF(IEBND.EQ.1)THEN   ! FLK ���������i���k�j�̌W��
      FLK=1.0D0
    ELSEIF(IEBND.EQ.2)THEN
      FLK=0.65D0
    ELSEIF(IEBND.EQ.3)THEN
      FLK=0.8D0
    ENDIF

    LOAD_TYPE=1 ! LOAD_TYPE �׏d�̎�ށA=0 �����׏d�A=1 �Z���׏d
    IF(LOAD_TYPE.EQ.0)THEN  ! LOAD_TYPE �׏d�̎�ށA=0 �����׏d�A=1 �Z���׏d
      SAFE_FACT=1.5     ! SAFE_FACT �������S��1.5
    ELSE
      SAFE_FACT=1.0     ! SAFE_FACT �Z�����S��1.0
    ENDIF

    LT=0.0  ! LT �ό`�㕔�ޒ�
    DO IDIM=1,3
      LT=LT+(COORD(N2,IDIM)-COORD(N1,IDIM))**2
    ENDDO
    LT=SQRT(LT)*1.0D+3

    AN=PI*(DMT**2-DMT1**2)/4.0          ! AN ���k�p�f�ʐ� mm2
    AT=AN                               ! AT �����p�f�ʐ� mm2
    ASZ=AN/1910D0*1013D0    ! �Z�莮�s���ɂ��b�� ! ASZ ����f�f�ʐςy mm2
    ASY=ASZ                 ! �Z�莮�s���ɂ��b�� ! ASY ����f�f�ʐςx mm2
    J=PI*RAD**3*THK/(1.0D0+POI)         ! J �˂��荄���W��
    JT=J/RAD                            ! JT �˂����R�W��
    ZZT=PI*(DMT**4-DMT1**4)/32.0/DMT    ! ZZT �f�ʌW���y�㕔 mm3
    ZZB=ZZT                             ! ZZB �f�ʌW���y���� mm3
    ZYL=ZZT                             ! ZYL �f�ʌW���x���� mm3
    ZYR=ZYL                             ! ZYR �f�ʌW���x�E�� mm3
    SIZ=SQRT(DMT**2+DMT1**2)/4.0        ! SIZ �f��2�����a�y mm
    SIY=SIZ                             ! SIY �f��2�����a�x mm
    LKZ=LT*FLK                          ! LKZ ���������y mm
    LKY=LKZ                             ! LKY ���������x mm
    !LBZ=?mm                            ! LBZ �Ȃ������v�Z�p�����y mm
    !LBY=?mm                            ! LBY �Ȃ������v�Z�p�����x mm

IF(IV.EQ.1)THEN ! �~�ʃA�[�`�̕��ޒ���
   LKZ=2500.0D0
   LKY=2500.0D0
ENDIF

    LMDZ=LKZ/SIZ                        ! LMDZ �ג���y ��z
    LMDY=LKY/SIZ                        ! LMDY �ג���x ��y
    RATIO=MIN(LMDZ,LMDY)/LLMD           ! ��/��
    IF(RATIO.LE.1.0)THEN                ! �Ɂ���
      VALNY=3.0/2.0+2.0/3.0*(RATIO)**2  ! VALNY ���S�� ��
      FC=(1.0-0.4*RATIO**2)*SST/VALNY*(1.5/SAFE_FACT)   ! FC ���e���k���͓x MPa
    ELSE                                ! �Ɂ���
      VALNY=2.17
      FC=0.277*SST/RATIO**2
    ENDIF
    FT=SST/SAFE_FACT                    ! FT ���e�������͓x MPa
    FBZ=SST                             ! FBZ ���e�Ȃ����͓x�y MPa
    FBY=SST                             ! FBY ���e�Ȃ����͓x�x MPa
    FSZ=SST/SQRT(3.0)/SAFE_FACT         ! FSZ ���e����f���͓x�y MPa
    FSY=SST/SQRT(3.0)/SAFE_FACT         ! FSY=���e����f���͓x�x MPa

    N=-STRSF(IE,1)*9.81*1.0D-3          ! N ���� kN
    SND=N*1.0D+3/AN                     ! SND ���������͓x MPa
    IF(SND.LE.0.0)THEN
      RND=-SND/FC                       ! RND ���͔���l�i���k)
    ELSE
      RND=SND/FT                        ! RND ���͔���l�i����)
    ENDIF
    MZ=MAX(ABS(STRSF(IE,5)),ABS(STRSF(IE,11)))*9.81*1.0D-3      ! MZ �Ȃ����[�����g�y kNm
    SBZ=MZ*1.0D+6/MIN(ZZT,ZZB)          ! SBZ �Ȃ����͓x�y MPa
    RBZ=SBZ/FBZ                         ! RBZ �Ȃ����͓x�y����l
    MY=MAX(ABS(STRSF(IE,6)),ABS(STRSF(IE,12)))*9.81*1.0D-3      ! MY �Ȃ����[�����g�x kNm
    SBY=MY*1.0D+6/MIN(ZYL,ZYR)          ! SBY �Ȃ����͓x�x MPa
    RBY=SBY/FBY                         ! RBY �Ȃ����͓x�x����l
    QZ=MAX(ABS(STRSF(IE,3)),ABS(STRSF(IE,9)))*9.81*1.0D-3       ! QZ ����f�͂y kN
    TAZ=QZ*1.0D+3/ASZ                   ! TAZ ����f���͓x�y MPa
    RTZ=TAZ/FSZ                         ! RTZ ����f���͓x�y����l
    QY=MAX(ABS(STRSF(IE,2)),ABS(STRSF(IE,8)))*9.81*1.0D-3       ! QY ����f�͂x kN
    TAY=QY*1.0D+3/ASY                   ! TAY ����f���͓x�x MPa
    RTY=TAY/FSY                         ! RTY ����f���͓x�x����l
    MJ=MAX(ABS(STRSF(IE,4)),ABS(STRSF(IE,10)))*9.81*1.0D-3      ! MJ �˂��胂�[�����g kNm
    TAJ=MJ*1.0D+6/JT                    ! TAJ ���肹��f���͓x MPa
    RTJ=TAJ/JT                          ! RTJ �˂��肹��f���͓x����l    
    RCB=SQRT((ABS(SND)+SQRT(SBZ**2+SBY**2))**2+3.0*(SQRT(TAZ**2+TAY**2)+TAJ)**2)/FT
!   RCB=SQRT((SND+SQRT(SBZ**2+SBY**2))**2)/FT   ! RCB �g�������͓x

! ���茋�ʂ̃t�@�C���o��

    IF(IE.EQ.1)CALL TEXT    ! �e�ϐ��̐������t�@�C���ɏo��

! ���ޖ���=CONCATENATE("P-",E55,"�Ӂ~",E56) ! ���ޖ���

    WRITE(14,401) IE
    WRITE(14,*)'�S�����ށi�|�ǁj�̌���'
    WRITE(14,*)''
    WRITE(14,402) DMT,THK
    WRITE(14,403) DMT
    WRITE(14,404) THK
    WRITE(14,*)''
    WRITE(14,*)'�ގ�      STK400'
    WRITE(14,405) SST
    WRITE(14,406) YNG
    WRITE(14,407) LLMD
401 FORMAT(///,'=== IELEF = ',I4,/)
402 FORMAT('���ޖ���',5X,'P-',F6.1,'�Ӂ~',F4.1)
403 FORMAT('            �O�a        ',F10.1,' mm')
404 FORMAT('            ����        ',F10.1,' mm')
405 FORMAT('            ����x    ',F10.1,' MPa')
406 FORMAT('            �����O��    ',E10.2,' MPa')
407 FORMAT('            ���E�ג��䃩',F10.1,'    ')
    WRITE(14,*)''
    WRITE(14,*)'�ޒ[�̍S������'
!   WRITE(14,'('            �ޒ[�S������Z =',I2,10X,'�ޒ[�S������Y =',I2)') 
    WRITE(14,411) IEBND
411 FORMAT('            �ޒ[�S������  =',I2)
    WRITE(14,*)''
    WRITE(14,*)'�׏d�̎��'
    IF(LOAD_TYPE.EQ.0) WRITE(14,412) SAFE_FACT
    IF(LOAD_TYPE.EQ.1) WRITE(14,413) SAFE_FACT
412 FORMAT('          �����׏d',/,'          ���S��',F15.1)
413 FORMAT('          �Z���׏d',/,'          ���S��',F15.1)
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
457 FORMAT('USE :  ','P-',F6.1,'�Ӂ~',F4.1,' -----> OK')
458 FORMAT('USE :  ','P-',F6.1,'�Ӂ~',F4.1,' -----> NG')

ENDDO ! LOOP OVER NELEF

RETURN
END

    SUBROUTINE TEXT     ! �ϐ��̐������o�͂���
    WRITE(14,*)'�S�����ށi�|�ǁj�̌���'
    WRITE(14,*)''
    WRITE(14,*)'���ޖ���'
    WRITE(14,*)'            �O�a'
    WRITE(14,*)'            ����'
    WRITE(14,*)''
    WRITE(14,*)'�ގ�      STK400'
    WRITE(14,*)'            ����x'
    WRITE(14,*)'            �����O��'
    WRITE(14,*)'            �|�A�\����'
    WRITE(14,*)'          ���E�ג��䃩'
    WRITE(14,*)''
    WRITE(14,*)'�ޒ[�̍S������'
    WRITE(14,*)'          �ޒ[�S������Z                  �ޒ[�S������Z'
    WRITE(14,*)'             =1�F�s���{�s���A=2�F���{���A=3�F�s���{���A=4�F���{�����[���[�A=5�F�s���{�����[���['
    WRITE(14,*)'          ���������i���k�j�̌W��Z        ���������i���k�j�̌W��Y'
    WRITE(14,*)''
    WRITE(14,*)'�׏d�̎��'
    WRITE(14,*)'          �Z��or����          1       =0 �����׏d�A =1 �Z���׏d'
    WRITE(14,*)'          ���S��           1.0'
    WRITE(14,*)''
    WRITE(14,*)'L0      = �ό`�O���ޒ�         �k�s    = �ό`�㕔�ޒ�'
    WRITE(14,*)'AN      = ���k�p�f�ʐ�         AT      = �����p�f�ʐ�'
    WRITE(14,*)'ASZ     = ����f�f�ʐ�Z        ASY     = ����f�f�ʐςx'
    WRITE(14,*)'J       = �˂��荄���W��       JT      = �˂����R�W��'
    WRITE(14,*)'ZZT     = �f�ʌW��Z�㕔        ZXB     = �f�ʌW���y����'
    WRITE(14,*)'ZYL     = �f�ʌW���x����       ZYR     = �f�ʌW���x�E��'
    WRITE(14,*)'SIZ     = �f��2�����a�y        SIY     = �f��2�����a�x'
    WRITE(14,*)'LKZ     = ���������y           LKY     = ���������x'
    WRITE(14,*)'LBZ     = �Ȃ������v�Z�p�����y LBY     = �Ȃ������v�Z�p�����x'
    WRITE(14,*)''
    WRITE(14,*)'LMDZ    = LKZ/SIZ =�ג���y    LMDY    = LKY/SIY =�ג���x'
    WRITE(14,*)'FC      = ���e���k���͓x       FT      = ���e�������͓x'
    WRITE(14,*)'FBZ     = ���e�Ȃ����͓x�y     FBY     = ���e�Ȃ����͓x�x'
    WRITE(14,*)'FSZ     = ���e����f���͓x�y   FSY     = ���e����f���͓x�x'
    WRITE(14,*)''
    WRITE(14,*)'N       = ����                 ST, SC  = N/AN     = ���������͓x'
    WRITE(14,*)'MZ      = �Ȃ����[�����g�y     SBZ     = MZ/ZZ    = �����Ȃ����͓x�y'
    WRITE(14,*)'MY      = �Ȃ����[�����g�x     SBY     = MY/ZY    = �����Ȃ����͓x�x'
    WRITE(14,*)'QXY     = ����f�͂y           TAZ     = QZ/ASZ  = ��������f���͓x�y'
    WRITE(14,*)'QY      = ����f�͂x           TAY     = QY/ASY  = ��������f���͓x�x'
    WRITE(14,*)'MJ      = �˂��胂�[�����g     TAJ     = MJ/JT    = �����˂��肹��f���͓x'
    WRITE(14,*)''
    WRITE(14,*)'�Z�莮'
    WRITE(14,*)'�Z�莮������l���P�D�O'
    WRITE(14,*)''
    WRITE(14,*)'USE�F�g�p���ށ@---->�@����'
    RETURN
    END