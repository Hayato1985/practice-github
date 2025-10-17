MODULE WRITE_VSCRIPT_MOD
    ! ======================================================================
    SUBROUTINE WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                             NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,NPROB,NTIMES,STRSF,STRSM,STRSC,TTDIS)
! ======================================================================
! �v�f�����}�Ɖ��͕��z�}�� VECTOR SCRIPT �`���̃t�@�C���ŏo�͂��邽�߂̃T�u���[�`��
!
! �T�u���[�`���̏����𐧌䂷�邽�߂̃p�����[�^
    IMPLICIT DOUBLEPRECISION (A-H,O-Z)
    PARAMETER(M_DIMSWITCH=0)        ! =1�̂Ƃ��͂x���W�Ƃy���W��������ĕ\������B
!    PARAMETER(M_PRNSTRESS=1)        ! �剞�͂̕\���`���i=0�F���ޗ��`��z�苖�e���͓x�敪�A=1�FETFE�t�B�����z��̍~�����͂ɂ��敪�j
    PARAMETER(M_ANGLE=1)            ! ���̖͂��\�������i=0�FIJ�����ɕ��s�A=1�F�剞�͕����j
    PARAMETER(M_PALETTE=0)          ! �����͂̕\���F�p���b�g�̐ݒ�i=0�F�ϓ��ɋ敪�A=1�F�������͏�Ԃ�z�肵�ċ敪�A=2�F���S��4.0��z�肵�ċ敪�j
    PARAMETER(M_ARROWDIM=3)         ! �v�f�����}�ł̎厲����ij���Q�����}�`�ŕ\���i=2�j�A�R�����i=3�j
!    PARAMETER(VAL_THICK=100.0)      ! ETFE�t�B�����z��iM_PRNSTRESS=1�j�̏ꍇ�̖���[��m]
    PARAMETER(VALST_MAX_FIX=0.0)    ! �����̃f�[�^�ɂ��Ē��͖��̒����𓝈ꂵ�����ꍇ�ɂ�VALST_MAX_FIX��0.0�ȊO�̐�����ݒ肷��
    PARAMETER(TMGNFY=2.0)           ! TMGNFY�F���͖��̒����𒲐߂���W���A2.0�̏ꍇ�ɂ͍ŏ����ډ~�̒��a�Ɉ�v����B
    PARAMETER(FMGNFY=1.0000)        ! FMGNFY�F�}�`�̊g��{���A�Ⴆ�΂��P�ʂ������P�ʂɕϊ�����ꍇ�ɂ�1000.0�Ƃ���΂悢�B
    PARAMETER(CMGNFY=0.8E+00)       ! CMGNFY�F���v�f��IJ������\������IJ�����ɑ΂���{��
!   PARAMETER(VALMARGIN=2.0)        ! VALMARGIN�F�}�`���݂̊Ԋu
    PARAMETER( NNODE = 3      )     ! NNODE �F���v�f���\������ߓ_�̐�
    DIMENSION COORD(NPOIN,3)        ! COORD(IPOIN,IDIM)�F�ߓ_���W
    DIMENSION COOR0(NPOIN,3)        ! COOR0(IPOIN,IDIM)�F�����`�󎞂̐ߓ_���W
    DIMENSION NODEM(NELEM,4)        ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
    DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)�F���v�f�̉��́i��R���k���]���j
    DIMENSION NODEC(NELEC,3)        !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
    DIMENSION STRSC(NELEC)          ! STRSC(IELEC)�F�P�[�u���̒���
    DIMENSION NODEF(NELEF,3)        ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
    DIMENSION STRSF(NELEF,NDOFN*2)  ! STRSF(IELEF,IEVAB)�F�Ȃ��v�fIELEF�̕��ޗ́i���ލ��W�n�ɂ����镔�ޒ[�׏d)
    DIMENSION TTDIS(NPOIN*NDOFN)    ! TTDIS(ITOTV)�F�����`��ɑ΂���ψ�
    DIMENSION VAL(3),VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(2,3)
    DIMENSION COD_MM(3,2)           ! COD_MM(1, )�F���W�̍ŏ��l�ACOD_MM(2, )�F�ő�l�A�ACOD_MM(3, )�F���AIDIM=1,2�i�w�A�x���W�j
    DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)�F���͂̋敪�ASTDIV(0,ISTRE)�F�ŏ��l�ASTDIV(9,ISTRE)�F�ő�l�AISTRE=4�F�剞�͂̏ꍇ
    CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)�F�J���[�p���b�g
    INTEGER,ALLOCATABLE :: NFILL(:,:)   ! NFILL(IELEM,ISTRE)=IDIV�F�p���b�g�̔ԍ��O�`�X�i�O�̓����N�����O��\���j
    DOUBLE PRECISION,ALLOCATABLE :: COOR1(:,:)  ! COOR1(IPOIN,IDIM)�F�}�`��\������Ƃ��Ɏg�p������W�l
    DOUBLE PRECISION,ALLOCATABLE :: STRSM1(:,:)    ! STRSM1(IELEM,ISTRE)�F���v�f�̎剞�́i��R���k���]���j
!
IF(NPROB.GE.50000.AND.NPROB.LE.50999)THEN   ! ETFE�t�B�����̍ޗ�����`���ɑΉ�����ꍇ�iNPROB=50000�`50999�Ԃ̏ꍇ�j
  M_PRNSTRESS=1        ! �剞�͂̕\���`���i=1�FETFE�t�B�����z��̍~�����͂ɂ��敪�j
  VAL_THICK=FLOAT(NPROB-50000)
ELSE
  M_PRNSTRESS=0        ! �剞�͂̕\���`���i=0�F���ޗ��`��z�苖�e���͓x�敪�j
ENDIF
!    OPEN(7,FILE='sample.txt')
!    READ(7,*)NPOIN,NELEM,NELEC,NELEF
!    LOADCREM=1  ! LOADCREM
!    NTIMES=1    ! NTIMES
!==========================================================
ALLOCATE(NFILL(NELEM,4),STAT=IERROR)        ! NFILL(IELEM,ISTRE)=IDIV�F�p���b�g�̔ԍ��O�`�X�i�O�̓����N�����O��\���j�AISTRE=4�F�剞�͂̏ꍇ
ALLOCATE(COOR1(NPOIN,3),STAT=IERROR)        ! COOR1(IPOIN,IDIM)�F�}�`��\������Ƃ��Ɏg�p������W�l
ALLOCATE(STRSM1(NELEM,MSTRE),STAT=IERROR)   ! STRSM1(IELEM,ISTRE)�F���v�f�̎剞�́i��R���k���]���j
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
! *** �}�`��\������Ƃ��Ɏg�p������W�lCOOR1��ݒ肷��
  DO IPOIN=1,NPOIN
    DO IDIM=1,3
!     COOR1(IPOIN,IDIM)=COOR0(IPOIN,IDIM)*FMGNFY
      COOR1(IPOIN,IDIM)=COORD(IPOIN,IDIM)*FMGNFY
    ENDDO
    ! M_DIMSWITCH=1�̂Ƃ��͂x���W�Ƃy���W���������B
    IF(M_DIMSWITCH.EQ.1)THEN
      COOR_TEMP=COOR1(IPOIN,2)
      COOR1(IPOIN,2)=COOR1(IPOIN,3)     ! �y���x
      COOR1(IPOIN,3)=-COOR_TEMP         ! �x���|�y�i�E��n��ۂ��߂Ƀ}�C�i�X�Ƃ���j
    ENDIF
  ENDDO
!==========================================================
! *** �w�A�x���W�̍ŏ��l�A�ő�l�A���̌v�Z
  DO IDIM=1,2
    COD_MM(1,IDIM)=COOR1(1,IDIM)    ! COD_MM(1, )�F���W�̍ŏ��l
    COD_MM(2,IDIM)=COOR1(1,IDIM)    ! COD_MM(2, )�F���W�̍ő�l
    DO IPOIN=2,NPOIN
      IF(COOR1(IPOIN,IDIM).LT.COD_MM(1,IDIM))COD_MM(1,IDIM)=COOR1(IPOIN,IDIM)
      IF(COOR1(IPOIN,IDIM).GT.COD_MM(2,IDIM))COD_MM(2,IDIM)=COOR1(IPOIN,IDIM)
    ENDDO
    COD_MM(3,IDIM)=COD_MM(2,IDIM)-COD_MM(1,IDIM)    ! COD_MM(3, )�F��
  ENDDO
  DO IDIM=1,3
    DO JDIM=1,2
      COD_MM(IDIM,JDIM)=COD_MM(IDIM,JDIM)*FMGNFY
    ENDDO
  ENDDO
!==========================================================
! *** �}�`���݂̊Ԋu�A�����T�C�Y�̐ݒ�
  VALMARGIN=COD_MM(3,1)/6.0
  TEXTSIZE=6.0
!==========================================================
! *** LOADCREM=0 or 1 �̏ꍇ�́A�w�b�_�[�Ɨv�f�����}���o�͂���
IF(LOADCREM.EQ.0.OR.LOADCREM.EQ.1)THEN
!==========================================================
! *** �f�[�^���o�͂���t�@�C���̐ݒ�
! OPEN(12,FILE='vector_script.txt')
!==========================================================
! *** �w�b�_�[�̏o��
  CALL vscript_header()
!==========================================================
! *** LAYER�i�v�f�����}�j�̏o��
  WRITE(12,*) '{Layer Characteristics}'
  WRITE(12,*) ''
  WRITE(12,*) 'Layer(''�v�f�����}'');'
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
! *** ���v�f�̏o��
  CALL vscript_element_division(1,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** �厲����ij��\�������̏o��
  CALL vscript_element_division(2,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** ���ޗv�f�̏o��
  CALL vscript_element_division(3,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** �ߓ_�ԍ��̏o��
  CALL vscript_element_division(4,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
! *** ���v�f�ԍ��̏o��
  CALL vscript_element_division(5,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!========================================================== 
! *** ���ޗv�f�ԍ��̏o��
  CALL vscript_element_division(6,CMGNFY,COOR1,NODEC,M_ARROWDIM,NODEM,NODEF,NELEF,NELEM,NELEC,NNODE,NPOIN,TEXTSIZE)
!==========================================================
ENDIF
!
!
!==========================================================
! *** LOADCREM��1 �̏ꍇ�́A���͕��z�}���o�͂���
IF(LOADCREM.GE.1) THEN
  DO IPOIN=1,NPOIN
    DO IDIM=1,3
      COOR1(IPOIN,IDIM)=COORD(IPOIN,IDIM)*FMGNFY
    ENDDO
    ! M_DIMSWITCH=1�̂Ƃ��͂x���W�Ƃy���W���������B
    IF(M_DIMSWITCH.EQ.1)THEN
      COOR_TEMP=COOR1(IPOIN,2)
      COOR1(IPOIN,2)=COOR1(IPOIN,3)     ! �y���x
      COOR1(IPOIN,3)=-COOR_TEMP         ! �x���|�y�i�E��n��ۂ��߂Ƀ}�C�i�X�Ƃ���j
    ENDIF
  ENDDO
!==========================================================
! *** LAYER�icase-xx�j�̏o��
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
  WRITE(12,*) 'SetView(#0� 0'' 0" ,#0� 0'' 0" ,#0� 0'' 0" ,0,0,0);'
  WRITE(12,*) ''
  WRITE(12,*) '{End of Layer Characteristics}'
  WRITE(12,*) ''
  WRITE(12,*) '{Object Creation Code}'
  WRITE(12,*) ''
!==========================================================
IF(NELEM.GT.0)THEN
! *** �����͂̕\���F�p���b�g�̐ݒ�
  CALL vscript_palette(M_PALETTE,MSTRE,NFILL,FILCOLOR,NELEM,STDIV,STRSM)
!==========================================================
! *** �����́i�������j�̏o��
  CALL vscript_stress(1,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** �����́i�������j�̏o��
  CALL vscript_stress(2,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** ���̂���f�́i�����j�̏o��
  CALL vscript_stress(3,COOR1,FILCOLOR,NPOIN,NELEM,NFILL,NNODE,NODEM,VALMARGIN,COD_MM(3,2))
!==========================================================
! *** �敪�Q�����́i�������j�̏o��
  CALL vscript_stress_div(1,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** �敪�Q�����́i�������j�̏o��
  CALL vscript_stress_div(2,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** �敪�Q���̂���f�́i�����j�̏o��
  CALL vscript_stress_div(3,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
! *** ���̎剞�͂���ђ��͖��̏o��
  CALL vscript_stress_principal(M_ANGLE,M_PRNSTRESS,VAL_THICK,COD_MM,COOR1,FILCOLOR,MSTRE,NELEM,NFILL,NNODE,NODEM,NPOIN,STDIV,STRSM,STRSM1,TMGNFY,VALMARGIN,VALST_MAX_FIX)
!==========================================================
! *** �敪�Q���̎剞�͂̏o��
  CALL vscript_stress_div(4,FILCOLOR,COD_MM,STDIV,VALMARGIN)
!==========================================================
ENDIF
! *** ���މ��́i����)�̏o��
IF(NELEC.GT.0)&
  CALL vscript_stress_axialforce(COD_MM,COOR1,NODEC,NDOFN,NELEF,NELEC,NODEF,NPOIN,STRSF,STRSC,TEXTSIZE,VALMARGIN)
!==========================================================
IF(NELEF.GT.0)THEN
! *** ���މ��́i�P�Fx���́A�Q�Fy����f�́A�R�Fz����f�́A�S�Fxx�˂���A�T�Fyy�Ȃ��A�U�Fzz�Ȃ�)�̏o��
  !���ނ̒����̍ő�lVALLNG_MAX�̌v�Z
  VALLNG_MAX=0.0
  DO IELEF=1,NELEF
    IPOIN=NODEF(IELEF,1)
    JPOIN=NODEF(IELEF,2)
    VALLNG=SQRT((COOR1(JPOIN,1)-COOR1(IPOIN,1))**2+(COOR1(JPOIN,2)-COOR1(IPOIN,2))**2+(COOR1(JPOIN,3)-COOR1(IPOIN,3))**2)
    IF(VALLNG_MAX.LT.VALLNG)VALLNG_MAX=VALLNG
  ENDDO
!���ޗ͂̍ő�lFORCE_MEX,MOMENT_MAX�̌v�Z
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
    IF(STRSF_MAX.EQ.0.0)STRSF_MAX=1.0   ! STRSF_MAX=0�̏ꍇ�̃G���[�iZero Divide�j�΍�
IF(ISTRE.EQ.1.OR.ISTRE.EQ.5.OR.ISTRE.EQ.6)THEN  !���́A�����Ȃ��A�㎲�Ȃ��̂ݏo�͂���
    CALL vscript_stress_frame(ISTRE,COD_MM,COOR0,COOR1,NVARF,NDOFN,NELEF,NODEF,NPOIN,PROPF,STRSF,STRSF_MAX,VALLNG_MAX,TEXTSIZE,TTDIS,VALMARGIN)
ENDIF
  ENDDO
ENDIF
!==========================================================
ENDIF
!==========================================================
! *** LOADCASE=0 or NTIMES(�Ō�)�̏ꍇ�́A�N���X���o�͂��ăt�@�C�������
IF(LOADCREM.EQ.0.OR.LOADCREM.EQ.NTIMES)THEN
  CALL vscript_class()    ! �N���X�̏o��
  CLOSE(12)
ENDIF
!
RETURN
END
!
! ======================================================================
    SUBROUTINE vscript_header()    ! �w�b�_�[�̏o��
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
WRITE(12,*) 'Layer(''�v�f�����}'');'
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
WRITE(12,*) 'hatchName:= ''�f�t�H���g�n�b�`���O'';'
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
! INDEX=1�̏ꍇ�F���v�f�̏o��
! INDEX=2�̏ꍇ�F�厲����ij��\�������̏o��
! INDEX=3�̏ꍇ�F���ޗv�f�̏o��
! INDEX=4�̏ꍇ�F�ߓ_�ԍ��̏o��
! INDEX=5�̏ꍇ�F���v�f�ԍ��̏o��
! INDEX=6�̏ꍇ�F���ޗv�f�ԍ��̏o��
! PARAMETER(M_ARROWDIM=)        ! �v�f�����}�ł̎厲����ij���Q�����}�`�ŕ\���i=2�j�A�R�����i=3�j
! PARAMETER(CMGNFY=)            ! CMGNFY�F���v�f��IJ������\������IJ�����ɑ΂���{��
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR0(IPOIN,IDIM)�F�����`�󎞂̐ߓ_���W
  DIMENSION NODEM(NELEM,4)        ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
  DIMENSION NODEC(NELEC,3)         !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
  DIMENSION NODEF(NELEF,3)        ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
  DIMENSION ARROW(3,3)            ! ���v�f��IJ������\�����̎n�_(1, )�A�I�_(2, )�A���v�f�̏d�S�_(3, )
  DIMENSION VEC1(2)
! *** ���v�f�̏o��
IF(INDEX.EQ.1)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''��_0 �v�f'');'
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
! *** �厲����ij��\�������̏o��
ELSEIF(INDEX.EQ.2)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''�厲����ij'');'
    WRITE(12,*) 'FPatByClass;'
    WRITE(12,*) 'FillColorByClass;'
    WRITE(12,*) 'LSByClass;'
    WRITE(12,*) 'PenColorByClass;'
    WRITE(12,*) 'LWByClass;'
    WRITE(12,*) 'MarkerByClass;'
  DO IELEM=1,NELEM
    DO IDIM=1,3
      ARROW(3,IDIM)=(COOR1(NODEM(IELEM,1),IDIM)+COOR1(NODEM(IELEM,2),IDIM)+COOR1(NODEM(IELEM,3),IDIM))/3    ! �d�S�_�̍��W�u3
      ARROW(1,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*COOR1(NODEM(IELEM,1),IDIM)
      ARROW(2,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*COOR1(NODEM(IELEM,2),IDIM)
      ARROW(3,IDIM)=(1.0-CMGNFY)*ARROW(3,IDIM)+CMGNFY*ARROW(2,IDIM)     ! ���̒������u32 ��(1-m)�{�ɂ���
    ENDDO
    IF(M_ARROWDIM.EQ.2)THEN ! �Q�����}�`�ŕ\������ꍇ 
      WRITE(12,*) 'MoveTo(',ARROW(1,1),',',ARROW(1,2),');'
      WRITE(12,*) 'LineTo(',ARROW(2,1),',',ARROW(2,2),');'
    ELSE    ! �R�����}�`�ŕ\������ꍇ�i�N���X�ݒ肪�_���Ȃ̂ŁAVector Script���捞�񂾌�ɑ����������ɕύX����Ε\�������j
      WRITE(12,*) 'OpenPoly;'
      WRITE(12,*) 'Poly3D('
      WRITE(12,1001)((ARROW(INODE,IDIM),IDIM=1,3),INODE=1,3)
      WRITE(12,1002)(ARROW(3,IDIM),IDIM=1,3)
    ENDIF
  ENDDO
!==========================================================
! *** ���ޗv�f�̏o��
ELSEIF(INDEX.EQ.3)THEN
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''����-0 �v�f'');'
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
! *** �ߓ_�ԍ��̏o��
ELSEIF(INDEX.EQ.4)THEN
  WRITE(12,*) 'NameClass(''�ԍ�-0 �ߓ_'');'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillFore(0,0,0);'
  WRITE(12,*) 'FillBack(65535,65535,13107);'
  WRITE(12,*) 'PenFore(0,0,0);'
  WRITE(12,*) 'PenBack(65535,65535,65535);'
  WRITE(12,*) 'TextFont(GetFontID(''�l�r �S�V�b�N''));'
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
! *** ���v�f�ԍ��̏o��
ELSEIF(INDEX.EQ.5)THEN
  WRITE(12,*) 'NameClass(''�ԍ�-1 ���v�f'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'PenFore(0,0,54272);'
  DO IELEM=1,NELEM
    VEC1(1)=(COOR1(NODEM(IELEM,1),1)+COOR1(NODEM(IELEM,2),1)+COOR1(NODEM(IELEM,3),1))/3.0   !VEC1()�F �O�p�`�̏d�S�̍��W
    VEC1(2)=(COOR1(NODEM(IELEM,1),2)+COOR1(NODEM(IELEM,2),2)+COOR1(NODEM(IELEM,3),2))/3.0   !VEC1()�F �O�p�`�̏d�S�̍��W
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
! *** ���ޗv�f�ԍ��̏o��
ELSEIF(INDEX.EQ.6)THEN
  WRITE(12,*) 'NameClass(''�ԍ�-2 ����'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
  WRITE(12,*) 'PenFore(56797,0,0);'
! �P�[�u���v�f�i���͗v�f�j�̏ꍇ
  DO IELEC=1,NELEC
    VEC1(1)=(COOR1(NODEC(IELEC,1),1)+COOR1(NODEC(IELEC,2),1))/2.0   !VEC1()�F �P�[�u���v�f�̒��_�̍��W
    VEC1(2)=(COOR1(NODEC(IELEC,1),2)+COOR1(NODEC(IELEC,2),2))/2.0   !VEC1()�F �P�[�u���v�f�̒��_�̍��W
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
! �Ȃ��v�f�̏ꍇ�i�J�M���ʂ��̔ԍ��ŕ\�� <1> �j
  DO IELEF=1,NELEF
    VEC1(1)=(COOR1(NODEF(IELEF,1),1)+COOR1(NODEF(IELEF,2),1))/2.0   !VEC1()�F �Ȃ��v�f�̒��_�̍��W
    VEC1(2)=(COOR1(NODEF(IELEF,1),2)+COOR1(NODEF(IELEF,2),2))/2.0   !VEC1()�F �Ȃ��v�f�̒��_�̍��W
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
! PARAMETER(M_PALETTE=0)          ! �����͂̕\���F�p���b�g�̐ݒ�i=0�F�ϓ��ɋ敪�A=1�F�������͏�Ԃ�z�肵�ċ敪�A=2�F���S��4.0��z�肵�ċ敪�j
IMPLICIT DOUBLEPRECISION (A-H,O-Z)
DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)�F���v�f�̒��́i��R���k���]���j
DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)�F���͂̋敪�ASTDIV(0,ISTRE)�F�ŏ��l�ASTDIV(9,ISTRE)�F�ő�l�AISTRE=4�F�剞�͂̏ꍇ
INTEGER NFILL(NELEM,4)          ! NFILL(IELEM,ISTRE)=IDIV�F�p���b�g�̔ԍ��O�`�X�i�O�̓����N�����O��\���j�AISTRE=4�F�剞�͂̏ꍇ
CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)�F�J���[�p���b�g
!
! �J���[�p���b�g�̐ݒ�
  IF(M_PALETTE.NE.1)THEN    ! M_PALETTE=0�F�ϓ��ɋ敪 OR INDEX=2�F���S��4.0��z�肵�ċ敪
    FILCOLOR(9) = '(56797,00000,00000)'     ! STDIV(8, )��    ��STDIV(9, )
    FILCOLOR(8) = '(65535,26214,13107)'     ! STDIV(7, )��    ��STDIV(8, )
    FILCOLOR(7) = '(65535,65535,13107)'     ! STDIV(6, )��    ��STDIV(7, )
    FILCOLOR(6) = '(00000,34952,00000)'     ! STDIV(5, )��    ��STDIV(6, )
    FILCOLOR(5) = '(00000,65535,00000)'     ! STDIV(4, )��    ��STDIV(5, )
    FILCOLOR(4) = '(00000,00000,65535)'     ! STDIV(3, )��    ��STDIV(4, )
    FILCOLOR(3) = '(42662,51914,61680)'     ! STDIV(2, )��    ��STDIV(3, )
    FILCOLOR(2) = '(34952,34952,34952)'     ! STDIV(1, )��    ��STDIV(2, )
    FILCOLOR(1) = '(48059,48059,48059)'     ! STDIV(0, )��    ��STDIV(1, )
    FILCOLOR(0) = '(65535,65535,65535)'     ! �����N�����O
  ELSE  ! M_PALETTE=1�F�������͏�Ԃ�z�肵�ċ敪
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
! �����͂̍ő�l�A�ŏ��l�����߂�
  DO IELEM=1,NELEM
    DO ISTRE=1,2
      IF(IELEM.EQ.1)THEN
        STDIV(9,ISTRE)=STRSM(1,ISTRE)    ! STDIV(9,ISTRE)�F�����͂̍ő�l
        STDIV(0,ISTRE)=STRSM(1,ISTRE)    ! STDIV(0,ISTRE)�F�����͂̍ŏ��l
      ELSE
        IF(STDIV(9,ISTRE).LT.STRSM(IELEM,ISTRE))STDIV(9,ISTRE)=STRSM(IELEM,ISTRE)
        IF(STDIV(0,ISTRE).GT.STRSM(IELEM,ISTRE))STDIV(0,ISTRE)=STRSM(IELEM,ISTRE)
      ENDIF
    ENDDO
    IF(IELEM.EQ.1)THEN
      STDIV(9,3)=ABS(STRSM(1,3))    ! STDIV(9,3)�F����f�͂̍ő�l
      STDIV(0,3)=ABS(STRSM(1,3))    ! STDIV(0,3)�F����f�͂̍ŏ��l
    ELSE
      IF(STDIV(9,3).LT.ABS(STRSM(IELEM,3)))STDIV(9,3)=ABS(STRSM(IELEM,3))
      IF(STDIV(0,3).GT.ABS(STRSM(IELEM,3)))STDIV(0,3)=ABS(STRSM(IELEM,3))
    ENDIF
  ENDDO
!
! �����͂̋敪�̊���͂����߂�(9�����j
  IF(M_PALETTE.EQ.0)THEN
!   �ϓ��ɋ敪
    DO IDIV=1,8
      DO ISTRE=1,3
        STDIV(IDIV,ISTRE)=STDIV(0,ISTRE)+(STDIV(9,ISTRE)-STDIV(0,ISTRE))/16.0*(2*IDIV-1)
      ENDDO
    ENDDO
  ELSEIF(M_PALETTE.EQ.1)THEN
!   �������͏�Ԃ����C�������Ƃ������敪�Ŏw�肵�ĕ\���������ꍇ�̏���
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
!   ���S��4.0����؂�Ƃ��ĕ\���������ꍇ�̏���
    STDIV(6,1)=3750.0    ! Tx(6) = 3750    '���ޗ��`��F�������x450kgf/3cm����
    STDIV(6,2)=3000.0    ! Ty(6) = 3000    '���ޗ��`��F�������x360kgf/3cm����
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
! �i�K���Ƃɖ��v�f�̐F�����߂�(10�i�K�j
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
! �����́ix�����Ay�������́A����f�́j�̏o��
! PARAMETER(VALMARGIN=)            ! VALMARGIN�F�}�`���݂̊Ԋu
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)�F�ߓ_���W
  DIMENSION NODEM(NELEM,4)        ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)�F�J���[�p���b�g
  INTEGER NFILL(NELEM,4)        ! NFILL(IELEM,ISTRE)=IDIV�F�p���b�g�̔ԍ��O�`�X�i�O�̓����N�����O��\���j�AISTRE=4�F�剞�͂̏ꍇ
  DIMENSION VEC(3,3)              ! VEC(INUM,IDIM)�F�\���p�̐ߓ_���W
! ISTRE�F=1 ���́ix�����j�A=2 ���́iy�����j�A=3 ����f��xy
! WIDTH_Y�F�}�`�̂x�����̕�
!
  WRITE(12,*) ''
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''��_1 ���́ix�����j'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''��_2 ���́iy�����j'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''��_3 ����f��xy'');'
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
! *** ���̎剞�͂���ђ��͖��̏o��
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)�F�ߓ_���W
  DIMENSION NODEM(NELEM,4)        ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
  DIMENSION STRSM(NELEM,MSTRE)    ! STRSM(IELEM,ISTRE)�F���v�f�̒��́i��R���k���]���j
  DIMENSION STRSM1(NELEM,MSTRE)   ! STRSM1(IELEM,ISTRE)�F���v�f�̎剞�́i��R���k���]���j
  DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)�F���͂̋敪�ASTDIV(0,ISTRE)�F�ŏ��l�ASTDIV(9,ISTRE)�F�ő�l�AISTRE=4�F�剞�͂̏ꍇ
  DIMENSION VAL(3),VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC(3,3)
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )�F���W�̍ŏ��l�ACOD_MM(2, )�F�ő�l�A�ACOD_MM(3, )�F���AIDIM=1,2�i�w�A�x���W�j
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)�F�J���[�p���b�g
  INTEGER NFILL(NELEM,4)          ! NFILL(IELEM,ISTRE)=IDIV�F�p���b�g�̔ԍ��O�`�X�i�O�̓����N�����O��\���j�AISTRE=4�F�剞�͂̏ꍇ
! PARAMETER(M_PRNSTRESS=)         ! �剞�͂̕\���`���i=0�F���ޗ��`��z�苖�e���͓x�敪�A=1�FETFE�t�B�����z��̍~�����͂ɂ��敪�j
! PARAMETER(M_ANGLE=)             ! ���̖͂��\�������i=0�FIJ�����ɕ��s�A=1�F�剞�͕����j
! PARAMETER(VAL_THICK=)           ! ETFE�t�B�����z��iM_PRNSTRESS=1�j�̏ꍇ�̖���[��m]
! PARAMETER(VALMARGIN=)           ! VALMARGIN�F�}�`���݂̊Ԋu
!
  OFFSET_Y=-(COD_MM(3,2)+VALMARGIN)*3     ! OFFSET_Y�F�}�`�\���ʒu�̂x�����ړ���
!================================================
! ���v�f�̎剞��STRSM1(IELEM,ISTRE)�����߂�
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
! VALR_MAX�F���ډ~�̔��a�̍ő�l�̌v�Z
  VALR_MAX=0.0
  DO IELEM=1,NELEM
    DO ISIDE=1,3
      N1=NODEM(IELEM,ISIDE)
      IF(ISIDE.LT.3)N2=NODEM(IELEM,ISIDE+1)
      IF(ISIDE.EQ.3)N2=NODEM(IELEM,1)
      VAL(ISIDE)=SQRT((COOR1(N2,1)-COOR1(N1,1))**2+(COOR1(N2,2)-COOR1(N1,2))**2+(COOR1(N2,3)-COOR1(N1,3))**2)
    ENDDO
    VALS=0.5*(VAL(1)+VAL(2)+VAL(3))
    VALSEC=SQRT(VALS*(VALS-VAL(1))*(VALS-VAL(2))*(VALS-VAL(3)))     ! VALSEC�F�O�p�`�̖ʐ�
    VALR=VALSEC/VALS   ! VALR�F���ډ~�̔��a
    IF(VALR.GT.VALR_MAX)VALR_MAX=VALR
  ENDDO
! VALST_MAX�F�����͂̍ő�l�̌v�Z
  IF(VALST_MAX_FIX.NE.0.0)THEN
    VALST_MAX=VALST_MAX_FIX     ! �����̃f�[�^�ɂ��Ė��̒����𓝈ꂵ�����ꍇ�ɂ�VALST_MAX���w��l�Ƃ���B
  ELSE
    VALST_MAX=0.0
    DO IELEM=1,NELEM
      IF(M_ANGLE.EQ.0)THEN      ! IJ�����ɕ��s�ɖ����o�͂���ꍇ
        DO ISTRE=1,2
          IF(VALST_MAX.LT.STRSM(IELEM,ISTRE))VALST_MAX=STRSM(IELEM,ISTRE)
        ENDDO
      ELSEIF(M_ANGLE.EQ.1)THEN  ! �剞�͕����̖����o�͂���ꍇ
        DO ISTRE=1,2
          IF(VALST_MAX.LT.STRSM1(IELEM,ISTRE))VALST_MAX=STRSM1(IELEM,ISTRE)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
! VALMFY�F���͖��̔{���̌v�Z
  IF(VALST_MAX.NE.0.0)VALMFY=VALR_MAX/VALST_MAX*TMGNFY
! ���͖��̏o��
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''��_4 ���͖��'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  DO IELEM=1,NELEM
    WRITE(12,*) 'MarkerByClass;'
!   WRITE(12,*) 'Marker(11,0.047241,8);'    ! ����󁩁��A����1.2mm/25.4=0.047241�A�p�x8��
    N1=NODEM(IELEM,1)
    N2=NODEM(IELEM,2)
    N3=NODEM(IELEM,3)
    VALLNG=SQRT((COOR1(N2,1)-COOR1(N1,1))**2+(COOR1(N2,2)-COOR1(N1,2))**2+(COOR1(N2,3)-COOR1(N1,3))**2) ! VALLNG�F�ӂP�Q�̒���
    DO IDIM=1,3
      VEC1(IDIM)=(COOR1(N1,IDIM)+COOR1(N2,IDIM)+COOR1(N3,IDIM))/3.0   !VEC1()�F �O�p�`�̏d�S�̍��W
      VEC2(IDIM)=(COOR1(N2,IDIM)-COOR1(N1,IDIM))/VALLNG     ! VEC2()�F�ӂP�Q�̒P�ʃx�N�g��
    ENDDO
    IF(M_ANGLE.EQ.0)THEN      ! IJ�����ɕ��s�ɖ����o�͂���ꍇ
      DO IDIM=1,3
        VEC3(IDIM)=VEC2(IDIM)*STRSM(IELEM,1)*VALMFY   ! VEC3()�F�ӂP�Q�����̒��͖��x�N�g��
      ENDDO
     VEC4(1)=-VEC2(2)*STRSM(IELEM,2)*VALMFY          ! VEC4()�F�ӂP�Q�ɒ�������̒��͖��x�N�g��
     VEC4(2)= VEC2(1)*STRSM(IELEM,2)*VALMFY
    ELSEIF(M_ANGLE.EQ.1)THEN  ! �剞�͕����̖����o�͂���ꍇ
      ANGLE1=STRSM1(IELEM,3)*3.14/180.0
      ANGLE2=(STRSM1(IELEM,3)+90.0)*3.14/180.0
      VEC3(1)=(VEC2(1)*COS(ANGLE1)-VEC2(2)*SIN(ANGLE1))*STRSM1(IELEM,1)*VALMFY   ! VEC3()�F�ӂP�Q�����̒��͖��x�N�g��
      VEC3(2)=(VEC2(1)*SIN(ANGLE1)+VEC2(2)*COS(ANGLE1))*STRSM1(IELEM,1)*VALMFY   ! VEC3()�F�ӂP�Q�����̒��͖��x�N�g��
      VEC4(1)=(VEC2(1)*COS(ANGLE2)-VEC2(2)*SIN(ANGLE2))*STRSM1(IELEM,2)*VALMFY   ! VEC4()�F�ӂP�Q�ɒ�������̒��͖��x�N�g��
      VEC4(2)=(VEC2(1)*SIN(ANGLE2)+VEC2(2)*COS(ANGLE2))*STRSM1(IELEM,2)*VALMFY   ! VEC4()�F�ӂP�Q�ɒ�������̒��͖��x�N�g��
    ENDIF
    WRITE(12,*) 'MoveTo(',VEC1(1)-VEC3(1)*0.5,',',VEC1(2)-VEC3(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'LineTo(',VEC1(1)+VEC3(1)*0.5,',',VEC1(2)+VEC3(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'MoveTo(',VEC1(1)-VEC4(1)*0.5,',',VEC1(2)-VEC4(2)*0.5+OFFSET_Y,');'
    WRITE(12,*) 'LineTo(',VEC1(1)+VEC4(1)*0.5,',',VEC1(2)+VEC4(2)*0.5+OFFSET_Y,');'
  ENDDO
!================================================
! ETFE�t�B�����z��̏ꍇ>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  IF(M_PRNSTRESS.EQ.1)THEN
    DO IELEM=1,NELEM
!     �������͂����߂�
      PX=STRSM(IELEM,1)
      PY=STRSM(IELEM,2)
      PXY=STRSM(IELEM,3)
      STRSM1(IELEM,1)=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)
!     ����[kg/m]������[N/mm2]�ɕϊ�����
      STRSM1(IELEM,1)=STRSM1(IELEM,1)*9.81/1000.0/(VAL_THICK*1.0E-3)
      STRSM1(IELEM,2)=STRSM1(IELEM,2)*9.81/1000.0/(VAL_THICK*1.0E-3)
    ENDDO
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! ���̎剞��STRSM1(IELEM,1)�̍ő�l�A�ŏ��l�����߂�
  DO IELEM=1,NELEM
    IF(IELEM.EQ.1)THEN
      STDIV(9,4)=STRSM1(1,1)        ! STDIV(9,4)�F�剞�͂̍ő�l
      STDIV(0,4)=STRSM1(1,1)        ! STDIV(0,4)�F�剞�͂̍ŏ��l
    ELSE
      IF(STDIV(9,4).LT.STRSM1(IELEM,1))STDIV(9,4)=STRSM1(IELEM,1)
      IF(STDIV(0,4).GT.STRSM1(IELEM,1))STDIV(0,4)=STRSM1(IELEM,1)
    ENDIF
  ENDDO
!
! �剞�͂̋敪�����߂�
    IF(M_PRNSTRESS.EQ.0)THEN
!       ���S��4.0����؂�Ƃ��ĕ\���������ꍇ�̏���
        STDIV(6,4)=3000.0    ! ���ޗ��`��F�������x360kgf/3cm����
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
!     ETFE�t�B�����z��̏ꍇ
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
!     �ϓ��ɋ敪����
      DO IDIV=1,8
        STDIV(IDIV,4)=STDIV(0,4)+(STDIV(9,4)-STDIV(0,4))/16.0*(2*IDIV-1)
      ENDDO
    ENDIF
!
! �i�K���Ƃɖ��v�f�̐F�����߂�(10�i�K�j
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
! ���v�f���o�͂���
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
! ���͋敪�̏o��
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION STDIV(0:9,4)          ! STDIV(IDIV,ISTRE)�F���͂̋敪�ASTDIV(0,ISTRE)�F�ŏ��l�ASTDIV(9,ISTRE)�F�ő�l�AISTRE=4�F�剞�͂̏ꍇ
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )�F���W�̍ŏ��l�ACOD_MM(2, )�F�ő�l�A�ACOD_MM(3, )�F���AIDIM=1,2�i�w�A�x���W�j
  CHARACTER*19 FILCOLOR(0:9)      ! FILCOLOR(IDIV)�F�J���[�p���b�g
  DIMENSION VEC(2,2)
! ISTRE�F=1 ���́ix�����j�A=2 ���́iy�����j�A=3 ����f��xy�A=4 �剞��
! VALMARGIN�F�}�`���݂̊Ԋu
!
! �a�n�w�T�C�Y�̐ݒ�
  TEXTSIZE=6.0
  IF(VALMARGIN*10.0.LE.COD_MM(3,2))THEN ! �a�n�w�̍������}�`�̏c�������������������ꍇ
    BOXSIZE=VALMARGIN*0.5
  ELSE
    BOXSIZE=COD_MM(3,2)/10.0
    TEXTSIZE=TEXTSIZE*BOXSIZE/VALMARGIN
  ENDIF
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''�敪-��_1 ����x'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''�敪-��_2 ����y'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''�敪-��_3 ����f��xy'');'
  IF(ISTRE.EQ.4)WRITE(12,*) 'NameClass(''��_4 ���͖��'');'
  WRITE(12,*) 'BeginGroup;'
! ���͋敪�̂a�n�w
  WRITE(12,*) 'FillPat(1);'
  DO IDIV=0,9   
    VEC(1,1)=COD_MM(2,1)+VALMARGIN
    VEC(2,1)=COD_MM(2,1)+VALMARGIN+BOXSIZE
    VEC(1,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*IDIV
    VEC(2,2)=COD_MM(1,2)-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)+BOXSIZE*0.5*(IDIV+1)
    WRITE(12,*) 'FillBack',FILCOLOR(IDIV),';'
    WRITE(12,*) 'Rect(',VEC(1,1),',',VEC(1,2),',',VEC(2,1),',',VEC(2,2),');'
  ENDDO
! ���͋敪�̐��� 
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
      1030 FORMAT(1H',I4,' �`',I4,1H')
    ENDIF
    WRITE(12,*) 'EndText;'
  ENDDO
! �^�C�g���̏o��
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
! *** ���މ��́i����)�̏o��
  IMPLICIT DOUBLEPRECISION (A-H,O-Z)
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)�F�}�`��\������Ƃ��Ɏg�p������W�l
  DIMENSION NODEC(NELEC,3)        !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
  DIMENSION STRSC(NELEC)          ! STRSC(IELEC)�F�P�[�u���̒���
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )�F���W�̍ŏ��l�ACOD_MM(2, )�F�ő�l�A�ACOD_MM(3, )�F���AIDIM=1,2�i�w�A�x���W�j
  DIMENSION VEC1(2)
! PARAMETER(VALMARGIN=)        ! VALMARGIN�F�}�`���݂̊Ԋu
!
  OFFSET_Y=-(COD_MM(3,2)+VALMARGIN)*4     ! OFFSET_Y�F�}�`�\���ʒu�̂x�����ړ���
  WRITE(12,*) ''
  WRITE(12,*) 'NameClass(''����-1 x����'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
! ���ނ̏o��
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
! ���͂̐��l�o��
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
    !���͂̐�������������ROTANGLE�̌v�Z
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
    !���͂̐��l�̕\�������̐ݒ�
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
! �^�C�g���̏o��
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
  PARAMETER(SMGNFY=0.10)          ! SMGNFY�F���ޗ͂̒����𒲐߂���W���A0.25�̏ꍇ�ɂ͂����Ƃ��������ނ�1/4�̒����ɂȂ�B
  DIMENSION COOR0(NPOIN,3)        ! COOR0(IPOIN,IDIM)�F�����`�󎞂̐ߓ_���W
  DIMENSION COOR1(NPOIN,3)        ! COOR1(IPOIN,IDIM)�F�ߓ_���W
  DIMENSION COD_MM(3,2)           ! COD_MM(1, )�F���W�̍ŏ��l�ACOD_MM(2, )�F�ő�l�A�ACOD_MM(3, )�F���AIDIM=1,2�i�w�A�x���W�j
  DIMENSION NODEF(NELEF,3)        ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
  DIMENSION PROPF(NVARF,10)       ! �Ȃ��v�f�̍ޗ������i��ޔԍ��A��������(1�`10)�j
  DIMENSION STRSF(NELEF,NDOFN*2)  ! STRSF(IELEF,IEVAB)�F�Ȃ��v�fIELEF�̕��ޗ́i���ލ��W�n�ɂ����镔�ޒ[�׏d)
  DIMENSION TTDIS(NPOIN*NDOFN)    ! TTDIS(ITOTV)�F�����`��ɑ΂���ψ�
  DIMENSION VEC1(3),VEC2(3),VEC3(3),VEC4(3),VEC5(3),VEC6(3),RTANG(3)
  DIMENSION OFFSET(3)             ! OFFSET()�F�}�`�\���ʒu�̂w�A�x�A�y�����ړ���
! ISTRE�FISTRE=�P�Fx���́A�Q�Fy����f�́A�R�Fz����f�́A�S�Fxx�˂���A�T�Fyy�Ȃ��A�U�Fzz�Ȃ�
! STRSF_MAX�F���͂̍ő�l
! VALLNG_MAX�F���ނ̒����̍ő�l
!
  WRITE(12,*) ''
  IF(ISTRE.EQ.1)WRITE(12,*) 'NameClass(''����-1 x����'');'
  IF(ISTRE.EQ.2)WRITE(12,*) 'NameClass(''����-2 y����f��'');'
  IF(ISTRE.EQ.3)WRITE(12,*) 'NameClass(''����-3 z����f��'');'
  IF(ISTRE.EQ.4)WRITE(12,*) 'NameClass(''����-4 xx�˂���'');'
  IF(ISTRE.EQ.5)WRITE(12,*) 'NameClass(''����-5 yy�Ȃ�'');'
  IF(ISTRE.EQ.6)WRITE(12,*) 'NameClass(''����-6 zz�Ȃ�'');'
  WRITE(12,*) 'FPatByClass;'
  WRITE(12,*) 'FillColorByClass;'
  WRITE(12,*) 'LSByClass;'
  WRITE(12,*) 'PenColorByClass;'
  WRITE(12,*) 'LWByClass;'
  WRITE(12,*) 'MarkerByClass;'
! �����̏����ݒ�
  WRITE(12,*) 'FillPat(0);'
  WRITE(12,*) 'FillFore(0,0,0);'
  WRITE(12,*) 'FillBack(65535,26214,13107);'
  WRITE(12,*) 'PenFore(0,0,0);'
  WRITE(12,*) 'PenBack(65535,65535,65535);'
  WRITE(12,*) 'TextSize(',INT(TEXTSIZE),');'
!
  OFFSET(1)=COD_MM(3,1)+VALMARGIN+5              ! OFFSET(1)�F�}�`�\���ʒu�̂w�����ړ���
  OFFSET(2)=-(COD_MM(3,2)+VALMARGIN)*(ISTRE-1)   ! OFFSET(2)�F�}�`�\���ʒu�̂x�����ړ���
  OFFSET(3)=0.0                                  ! OFFSET(3)�F�}�`�\���ʒu�̂y�����ړ���
  IEVAB1=ISTRE      ! IEVAB1�F���ނ̂��[���̉��͔ԍ�
  IEVAB2=ISTRE+6    ! IEVAB2�F���ނ̂��[���̉��͔ԍ�
! �^�C�g���̏o��
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
    INFRM=NODEF(IELEF,3)  ! INFRM�F�Ȃ��v�fIELEF�̕��ޔԍ�
!   ���ނ̏o��
    WRITE(12,*) 'Poly3D('
    IPOIN=NODEF(IELEF,1)
    JPOIN=NODEF(IELEF,2)
    VALLNG=SQRT((COOR1(JPOIN,1)-COOR1(IPOIN,1))**2+(COOR1(JPOIN,2)-COOR1(IPOIN,2))**2+(COOR1(JPOIN,3)-COOR1(IPOIN,3))**2)
    DO IDIM=1,3
      VEC1(IDIM)=(COOR0(JPOIN,IDIM)-COOR0(IPOIN,IDIM))/VALLNG       !VEC1( )�F�v�f���������̒P�ʃx�N�g��
      VEC2(IDIM)=PROPF(INFRM,IDIM+7)  !VEC2( )�F�㎲�����̒P�ʃx�N�g��
    ENDDO
    CALL VECTPRD(VEC1,VEC2,VEC3)    !VECT3( )�F���������̒P�ʃx�N�g����VECT1�~VECT2
    DO INODE=1,2
      DO IDIM=1,3
        KPOIN=NODEF(IELEF,INODE)
        ITOTV=(KPOIN-1)*NDOFN+IDIM+3
        RTANG(IDIM)=TTDIS(ITOTV)    ! RTANG(IDIM)�F�w�A�x�A�y���iIDIM=1�`3�j����̕ό`�pTTDIS(ITOTV)
      ENDDO
      IF(ISTRE.EQ.2.OR.ISTRE.EQ.6)THEN
        CALL ROTXYZ(RTANG,VEC2,VEC4)   ! VEC4( )�FVEC3( )���w�A�x�A�y������ɉ�]�����x�N�g��
      ELSE
        CALL ROTXYZ(RTANG,VEC3,VEC4)   ! VEC4( )�FVEC2( )���w�A�x�A�y������ɉ�]�����x�N�g��
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
!   ���͂̐��l�̏o��
    POSX=VEC1(1)
    POSY=VEC1(2)
    !���͂̐�������������ROTANGLE�̌v�Z
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
ROTANGLE=0.0    !�����̕����𑵂���
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
! VEC1���w�A�x�A�y�������RTANG�i1�`3)��]���āAVEC4���쐬����T�u���[�`��
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
!  CALL MLTPLY_ABC(3,3,3,1,RY,RZ,VEC1,VEC4)     ! VEC4( )��RY(RZ(VEC1))
  CALL MLTPLY(3,3,1,RX,VEC1,VEC4)   ! VEC4( )��RX(VEC4)
  CALL MLTPLY(3,3,1,RY,VEC4,VEC4)   ! VEC4( )��RX(VEC4)
  CALL MLTPLY(3,3,1,RZ,VEC4,VEC4)   ! VEC4( )��RX(VEC4)
  RETURN
  END

! ======================================================================
    SUBROUTINE vscript_class()    ! �N���X�̏o��
! ======================================================================
WRITE(12,*) '{End of Creation Code}'
WRITE(12,*) ''
WRITE(12,*) '{Classes}'
WRITE(12,*) ''
WRITE(12,*) 'NameClass(''���'');'
WRITE(12,*) 'SetClFillFore(''���'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''���'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''���'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''���'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''���'',1);'
WRITE(12,*) 'SetClLS(''���'',2);'
WRITE(12,*) 'SetClLW(''���'',1);'
WRITE(12,*) 'SetClUseGraphic(''���'',FALSE);'
WRITE(12,*) 'NameClass(''���@'');'
WRITE(12,*) 'SetClFillFore(''���@'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''���@'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''���@'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''���@'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''���@'',1);'
WRITE(12,*) 'SetClLS(''���@'',2);'
WRITE(12,*) 'SetClLW(''���@'',1);'
WRITE(12,*) 'SetClUseGraphic(''���@'',FALSE);'
WRITE(12,*) 'NameClass(''����-0 �v�f'');'
WRITE(12,*) 'SetClFillFore(''����-0 �v�f'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-0 �v�f'',0,0,0);'
WRITE(12,*) 'SetClPenFore(''����-0 �v�f'',56797,0,0);'
WRITE(12,*) 'SetClPenBack(''����-0 �v�f'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-0 �v�f'',0);'
WRITE(12,*) 'SetClLS(''����-0 �v�f'',2);'
WRITE(12,*) 'SetClLW(''����-0 �v�f'',6);'
WRITE(12,*) 'SetClUseGraphic(''����-0 �v�f'',TRUE);'
WRITE(12,*) 'NameClass(''��_0 �v�f'');'
WRITE(12,*) 'SetClFillFore(''��_0 �v�f'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''��_0 �v�f'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''��_0 �v�f'',53520,41058,29490);'
WRITE(12,*) 'SetClPenBack(''��_0 �v�f'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''��_0 �v�f'',1);'
WRITE(12,*) 'SetClLS(''��_0 �v�f'',2);'
WRITE(12,*) 'SetClLW(''��_0 �v�f'',3);'
WRITE(12,*) 'SetClUseGraphic(''��_0 �v�f'',TRUE);'
WRITE(12,*) 'NameClass(''��_1 ���́ix�����j'');'
WRITE(12,*) 'SetClFillFore(''��_1 ���́ix�����j'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''��_1 ���́ix�����j'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''��_1 ���́ix�����j'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''��_1 ���́ix�����j'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''��_1 ���́ix�����j'',1);'
WRITE(12,*) 'SetClLS(''��_1 ���́ix�����j'',2);'
WRITE(12,*) 'SetClLW(''��_1 ���́ix�����j'',3);'
WRITE(12,*) 'SetClUseGraphic(''��_1 ���́ix�����j'',FALSE);'
WRITE(12,*) 'NameClass(''����-2 y����f��'');'
WRITE(12,*) 'SetClFillFore(''����-2 y����f��'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-2 y����f��'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-2 y����f��'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-2 y����f��'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-2 y����f��'',0);'
WRITE(12,*) 'SetClLS(''����-2 y����f��'',2);'
WRITE(12,*) 'SetClLW(''����-2 y����f��'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-2 y����f��'',TRUE);'
WRITE(12,*) 'NameClass(''�敪-��_1 ����x'');'
WRITE(12,*) 'SetClFillFore(''�敪-��_1 ����x'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�敪-��_1 ����x'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''�敪-��_1 ����x'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''�敪-��_1 ����x'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�敪-��_1 ����x'',1);'
WRITE(12,*) 'SetClLS(''�敪-��_1 ����x'',2);'
WRITE(12,*) 'SetClLW(''�敪-��_1 ����x'',1);'
WRITE(12,*) 'SetClUseGraphic(''�敪-��_1 ����x'',FALSE);'
WRITE(12,*) 'NameClass(''�敪-��_2 ����y'');'
WRITE(12,*) 'SetClFillFore(''�敪-��_2 ����y'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�敪-��_2 ����y'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''�敪-��_2 ����y'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''�敪-��_2 ����y'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�敪-��_2 ����y'',1);'
WRITE(12,*) 'SetClLS(''�敪-��_2 ����y'',2);'
WRITE(12,*) 'SetClLW(''�敪-��_2 ����y'',1);'
WRITE(12,*) 'SetClUseGraphic(''�敪-��_2 ����y'',FALSE);'
WRITE(12,*) 'NameClass(''�敪-��_3 ����f��xy'');'
WRITE(12,*) 'SetClFillFore(''�敪-��_3 ����f��xy'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�敪-��_3 ����f��xy'',56797,0,0);'
WRITE(12,*) 'SetClPenFore(''�敪-��_3 ����f��xy'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''�敪-��_3 ����f��xy'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�敪-��_3 ����f��xy'',1);'
WRITE(12,*) 'SetClLS(''�敪-��_3 ����f��xy'',2);'
WRITE(12,*) 'SetClLW(''�敪-��_3 ����f��xy'',1);'
WRITE(12,*) 'SetClUseGraphic(''�敪-��_3 ����f��xy'',FALSE);'
WRITE(12,*) 'NameClass(''�厲����ij'');'
WRITE(12,*) 'SetClFillFore(''�厲����ij'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�厲����ij'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''�厲����ij'',22102,11308,1285);'
WRITE(12,*) 'SetClPenBack(''�厲����ij'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�厲����ij'',0);'
WRITE(12,*) 'SetClLS(''�厲����ij'',-2);'
WRITE(12,*) 'SetClLW(''�厲����ij'',3);'
WRITE(12,*) 'SetClUseGraphic(''�厲����ij'',TRUE);'
WRITE(12,*) 'NameClass(''��_3 ����f��xy'');'
WRITE(12,*) 'SetClFillFore(''��_3 ����f��xy'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''��_3 ����f��xy'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''��_3 ����f��xy'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''��_3 ����f��xy'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''��_3 ����f��xy'',1);'
WRITE(12,*) 'SetClLS(''��_3 ����f��xy'',2);'
WRITE(12,*) 'SetClLW(''��_3 ����f��xy'',3);'
WRITE(12,*) 'SetClUseGraphic(''��_3 ����f��xy'',FALSE);'
WRITE(12,*) 'NameClass(''����-3 z����f��'');'
WRITE(12,*) 'SetClFillFore(''����-3 z����f��'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-3 z����f��'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-3 z����f��'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-3 z����f��'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-3 z����f��'',0);'
WRITE(12,*) 'SetClLS(''����-3 z����f��'',2);'
WRITE(12,*) 'SetClLW(''����-3 z����f��'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-3 z����f��'',TRUE);'
WRITE(12,*) 'NameClass(''��_2 ���́iy�����j'');'
WRITE(12,*) 'SetClFillFore(''��_2 ���́iy�����j'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''��_2 ���́iy�����j'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''��_2 ���́iy�����j'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''��_2 ���́iy�����j'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''��_2 ���́iy�����j'',1);'
WRITE(12,*) 'SetClLS(''��_2 ���́iy�����j'',2);'
WRITE(12,*) 'SetClLW(''��_2 ���́iy�����j'',3);'
WRITE(12,*) 'SetClUseGraphic(''��_2 ���́iy�����j'',FALSE);'
WRITE(12,*) 'NameClass(''����-1 x����'');'
WRITE(12,*) 'SetClFillFore(''����-1 x����'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-1 x����'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-1 x����'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-1 x����'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-1 x����'',0);'
WRITE(12,*) 'SetClLS(''����-1 x����'',2);'
WRITE(12,*) 'SetClLW(''����-1 x����'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-1 x����'',TRUE);'
WRITE(12,*) 'NameClass(''��_4 ���͖��'');'
WRITE(12,*) 'SetClFillFore(''��_4 ���͖��'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''��_4 ���͖��'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''��_4 ���͖��'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''��_4 ���͖��'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''��_4 ���͖��'',1);'
WRITE(12,*) 'SetClLS(''��_4 ���͖��'',2);'
WRITE(12,*) 'SetClLW(''��_4 ���͖��'',3);'
WRITE(12,*) 'SetClUseGraphic(''��_4 ���͖��'',TRUE);'
WRITE(12,*) 'NameClass(''�ԍ�-0 �ߓ_'');'
WRITE(12,*) 'SetClFillFore(''�ԍ�-0 �ߓ_'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�ԍ�-0 �ߓ_'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''�ԍ�-0 �ߓ_'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''�ԍ�-0 �ߓ_'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�ԍ�-0 �ߓ_'',0);'
WRITE(12,*) 'SetClLS(''�ԍ�-0 �ߓ_'',2);'
WRITE(12,*) 'SetClLW(''�ԍ�-0 �ߓ_'',1);'
WRITE(12,*) 'SetClUseGraphic(''�ԍ�-0 �ߓ_'',TRUE);'
WRITE(12,*) 'NameClass(''�ԍ�-1 ���v�f'');'
WRITE(12,*) 'SetClFillFore(''�ԍ�-1 ���v�f'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�ԍ�-1 ���v�f'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''�ԍ�-1 ���v�f'',0,0,54272);'
WRITE(12,*) 'SetClPenBack(''�ԍ�-1 ���v�f'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�ԍ�-1 ���v�f'',0);'
WRITE(12,*) 'SetClLS(''�ԍ�-1 ���v�f'',2);'
WRITE(12,*) 'SetClLW(''�ԍ�-1 ���v�f'',1);'
WRITE(12,*) 'SetClUseGraphic(''�ԍ�-1 ���v�f'',TRUE);'
WRITE(12,*) 'NameClass(''�ԍ�-2 ����'');'
WRITE(12,*) 'SetClFillFore(''�ԍ�-2 ����'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''�ԍ�-2 ����'',65535,65535,13107);'
WRITE(12,*) 'SetClPenFore(''�ԍ�-2 ����'',56797,0,0);'
WRITE(12,*) 'SetClPenBack(''�ԍ�-2 ����'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''�ԍ�-2 ����'',0);'
WRITE(12,*) 'SetClLS(''�ԍ�-2 ����'',2);'
WRITE(12,*) 'SetClLW(''�ԍ�-2 ����'',1);'
WRITE(12,*) 'SetClUseGraphic(''�ԍ�-2 ����'',TRUE);'
WRITE(12,*) 'NameClass(''����-4 xx�˂���'');'
WRITE(12,*) 'SetClFillFore(''����-4 xx�˂���'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-4 xx�˂���'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-4 xx�˂���'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-4 xx�˂���'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-4 xx�˂���'',0);'
WRITE(12,*) 'SetClLS(''����-4 xx�˂���'',2);'
WRITE(12,*) 'SetClLW(''����-4 xx�˂���'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-4 xx�˂���'',TRUE);'
WRITE(12,*) 'NameClass(''����-5 yy�Ȃ�'');'
WRITE(12,*) 'SetClFillFore(''����-5 yy�Ȃ�'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-5 yy�Ȃ�'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-5 yy�Ȃ�'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-5 yy�Ȃ�'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-5 yy�Ȃ�'',0);'
WRITE(12,*) 'SetClLS(''����-5 yy�Ȃ�'',2);'
WRITE(12,*) 'SetClLW(''����-5 yy�Ȃ�'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-5 yy�Ȃ�'',TRUE);'
WRITE(12,*) 'NameClass(''����-6 zz�Ȃ�'');'
WRITE(12,*) 'SetClFillFore(''����-6 zz�Ȃ�'',0,0,0);'
WRITE(12,*) 'SetClFillBack(''����-6 zz�Ȃ�'',65535,65535,65535);'
WRITE(12,*) 'SetClPenFore(''����-6 zz�Ȃ�'',0,0,0);'
WRITE(12,*) 'SetClPenBack(''����-6 zz�Ȃ�'',65535,65535,65535);'
WRITE(12,*) 'SetClFPat(''����-6 zz�Ȃ�'',0);'
WRITE(12,*) 'SetClLS(''����-6 zz�Ȃ�'',2);'
WRITE(12,*) 'SetClLW(''����-6 zz�Ȃ�'',3);'
WRITE(12,*) 'SetClUseGraphic(''����-6 zz�Ȃ�'',TRUE);'
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
! ... (���̃T�u���[�`����֐�������΁A������END SUBROUTINE/FUNCTION�̌��)
END MODULE WRITE_VSCRIPT_MOD