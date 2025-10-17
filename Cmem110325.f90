   
    
    ! ======================================================================
!     NONLINEAR ANALYSIS OF CABLE REINFORCED MEMBRANE STRUCTURE    
!                                      - EXACT SOL. BY F.E.M.
!                                      - THE NEWTON RAPHSON METHOD
!                                      - SKYLINE MATRIX
!                                      - STRESS TRANSFER METHOD
! ======================================================================
! �v���O�����̍\��
! Cmem7.f90	         ���C���v���O�����A������уP�[�u���̍����}�g���N�X�쐬
! Allowance.f90      �S�����ށi�|�ǁj�̋��e���͓x����
! Arclm.f90          �ʒ������@�ɂ���͂̂��߂̃��[�`��
! EIGRS.FOR          �ŗL�l��̓��[�`��
! Frame.f90          �Ȃ��v�f�̍����}�g���N�X�쐬�A�Ȃ��v�f�֘A�f�[�^�̃t�@�C�����o��
! IOcntl.f90         INPUT�t�@�C������̓��́AOUTPUT�t�@�C���ւ̏o��
! loads.f90          ���d�A�O�́i���׏d�A��׏d�A�ߓ_�W���׏d�j�A�����̓����ߓ_�͂̌v�Z
! skylne2M.for       �X�J�C���C���@�ɂ��A���ꎟ�������̉�@���[�`��
! TOOLS.F90          �⏕�I�ȃ��[�`��
! write_vsctipt.f90  �v�f�����}�Ɖ��͕��z�}�� VECTOR SCRIPT �`���̃t�@�C���ŏo�͂��邽�߂̃T�u���[�`��
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ 121214 GRANPA DOME  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ IOcntl ���͌W���̌v�Z @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ Cmem�P�[�u���̃����O�W��PROPC @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
USE WRITE_VSCRIPT_MOD
implicit doubleprecision(A-H,O-Z)
character*15 text
!
!*** �����I�z��̐錾��
!
REAL*8 ATMOS_P PARAMETER (ATMOS_P=101325.0D0)
REAL*8 V_INIT
 
doubleprecision, allocatable :: CELL1(:,:)    !Excel�̓��̓f�[�^
doubleprecision, allocatable :: PMGNF(:,:)    !�׏d�{��
doubleprecision, allocatable :: PPRES(:)      !����[N/m2]
doubleprecision, allocatable :: WWIND(:)      !���x��[N/m2]
doubleprecision, allocatable :: SSNOW(:)      !��׏d[N/m2]
doubleprecision, allocatable :: F(:)          !���׏d�A��׏d�x�N�g��
doubleprecision, allocatable :: FF(:)         !�����̓x�N�g��
doubleprecision, allocatable :: FFF(:)        !�����׏d�x�N�g��
doubleprecision, allocatable :: FINIT(:)      !�Œ�׏d�x�N�g��
doubleprecision, allocatable :: FCOE(:)       !���ڐߓ_�׏d�̉׏d�W��
doubleprecision, allocatable :: Q(:)          !�c���̓x�N�g��
doubleprecision, allocatable :: PSUM(:)       !���ڐߓ_�O�͂̊e�����̍��v

doubleprecision, allocatable :: COORD(:,:)    !�ߓ_���W
doubleprecision, allocatable :: COOR0(:,:)    !�����`�󎞂̐ߓ_���W
doubleprecision, allocatable :: COOR1(:,:)    !�j���[�g���@�ɂ��J�Ԃ��v�Z�O�� �ߓ_���W
doubleprecision, allocatable :: DIS(:)        !�ψ�
doubleprecision, allocatable :: TTDIS(:)      !�����`��ɑ΂���ψ�

integer, allocatable         :: IFFIX(:)      !�S���Ɋւ�����i0:���R�A1:�S��)
integer, allocatable         :: NOFIX(:)      !�S���ߓ_�̔ԍ�
doubleprecision, allocatable :: PRESC(:,:)    !IVFIX�Ԗڂ̍S���ߓ_��IDOFN�Ԗڂ̎��R�x�̍S����
doubleprecision, allocatable :: PRESC0(:,:)   !IVFIX�Ԗڂ̍S���ߓ_��IDOFN�Ԗڂ̎��R�x�̍S���ʁi�ݒ�l���i�[���邽�߂̔z��j
integer, allocatable         :: IFPRE(:)      !�S����Ԃ�\���ϐ��i��@101:x,z�����S���Ay�������R�j

integer, allocatable         :: NODEM(:,:)    !���v�f�̐ߓ_�ԍ�
doubleprecision, allocatable :: PROPM(:,:,:)  !���v�f�̍ޗ������i��ޔԍ��A�������ځA�i�K�j
doubleprecision, allocatable :: STRSM(:,:)    !���v�f�̒���
doubleprecision, allocatable :: STRM1(:,:)    !�����N�����O�������s���O�̖��v�f�̒���
integer, allocatable         :: NANG(:)       !���v�fIJK��IJ�����ƍޗ��̎厲�������X���Ă���v�f�̔ԍ�
doubleprecision, allocatable :: ANG(:)        !I���v�fIJK��IJ�����ƍޗ��̎厲�����̂Ȃ��p�x
integer, allocatable         :: MEME1(:)      !�剞�͕����w�̃����N�����O�����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
integer, allocatable         :: MEME2(:)      !�剞�͕����x�̃����N�����O�����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
doubleprecision, allocatable :: AREA(:)       !���v�f�̏����`�󎞂̖ʐ�
doubleprecision, allocatable :: WCF(:)        !�����W��
doubleprecision, allocatable :: SCF(:)        !��׏d�W��

integer, allocatable         :: NODEC(:,:)    !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
doubleprecision, allocatable :: BL(:)         !�P�[�u���̒���
doubleprecision, allocatable :: CAL(:)        !�P�[�u���̖��Ђ��ݒ���
doubleprecision, allocatable :: STRSC(:)      !�P�[�u���̒���
doubleprecision, allocatable :: STRC1(:)      !���݂̏������s���O�̃P�[�u���̒���
integer, allocatable         :: NEA(:)        !�P�[�u���̂��݂̔����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
doubleprecision, allocatable :: PROPC(:,:)    !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj

doubleprecision, allocatable :: BETA(:)       !�Ȃ����ނ̉�]�p
doubleprecision, allocatable :: CALFR(:)      !���ނ̒���
doubleprecision, allocatable :: VALK(:,:)     !�ޒ[1,2�̐ڍ������i��1.0E+20�F���ځA��0.0�F�s���ځA���̑��F�����ځj
integer, allocatable         :: NODEF(:,:)    !�Ȃ��v�f���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
integer, allocatable         :: NPIN(:)       !���ޒ[���̐ڍ������i=0�F�����A=1�F�s�����A=2�F���s���A=3�F�s���s���j
doubleprecision, allocatable :: PROPF(:,:)    !�Ȃ��v�f�̍ޗ������i��ޔԍ��A�������ځj
doubleprecision, allocatable :: STRSF(:,:)    !�Ȃ��v�fIELEF�̕��ޗ́i����W�n�ɂ����镔�ޒ[�׏d)
doubleprecision, allocatable :: STRF1(:,:)    !�Ȃ��v�fIELEF�̕��ޗ́i���ލ��W�n�ɂ����镔�ޒ[�׏d)

doubleprecision, allocatable :: EE(:,:)       !�e���萔�}�g���N�X�c
doubleprecision, allocatable :: TSTIF(:)      !�S�̍����}�g���N�X�i�X�J�C���C���}�g���N�X�j
integer, allocatable         :: NORDR(:)      !�S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ����ĕ������鎞�̐V�����S�̎��R�x�ԍ�
integer, allocatable         :: NWDTH(:)      !�S�̍����}�g���N�X�̑�ITOTV �s�̑�P��[���v�f����Ίp�v�f�܂ł̌�
integer, allocatable         :: NWSUM(:)      !�S�̍����}�g���N�X�̑�ITOTV �s�̑�P��[���v�f����Ίp�v�f�܂ł̌�

! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>
doubleprecision, allocatable :: ELOAD(:)      !�c����
doubleprecision, allocatable :: TLOAD(:)      !���׏d
doubleprecision, allocatable :: ADIS1(:)      !���t'[k,j1]�ATLOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
doubleprecision, allocatable :: ADIS2(:)      !���t'[k,j2]�AELOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
doubleprecision, allocatable :: ASDIS(:)      !�t[k,(j)]�A��IINCS�񑝕��ɂ��ψʗ�
doubleprecision, allocatable :: RPRVS(:)      !��(IITER-1)�񔽕���̂��x�N�g��
doubleprecision, allocatable :: D(:)          !SUB SKYLNE�̍�Ɨp�z��
doubleprecision, allocatable :: T(:)          !SUB SKYLNE�̍�Ɨp�z��
integer, allocatable         :: NWK(:)        !SUB SKYLNE�̍�Ɨp�z��
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
doubleprecision, allocatable :: CTEMP1(:)     !CTEMP1(LOADCREM)�F�P�[�u���ɗ^���鉷�x����[��]�AMCTEMP=1�i�P�[�u���̉��x���͉�͂��s���j�̏ꍇ
doubleprecision, allocatable :: CTEMP2(:,:)   !CTEMP2(1,IELEC)�F���x���͂�^����(=1.0)�A�^���Ȃ�(=0.0)�ACTEMP2(2,IELEC)�F���x�Ђ���
!
! *** �����o�̓t�@�C���̐ݒ聄 ***
!
!open(7,file='J:\Source\Excel Cmem\Cmem110325\data\Cmem110325.prn',status='unknown')
!open(8,file='J:\Source\Excel Cmem\Cmem110325\data\Cmem110325_out.txt',status='unknown')
open( 7,file='Cmem110325.prn',status='old')
open( 8,file='Cmem110325_out.txt',status='unknown')
open( 9,file='Cmem110325.cod',status='unknown')
open(10,file='Cmem110325.ten',status='unknown')
open(11,file='Cmem110325_pro.txt',status='unknown')
open(12,file='Cmem110325_vector_script.txt',status='unknown')
open(13,file='Cmem110325_arm.txt',status='unknown')
!
!*** �����I�z��̊��蓖�ĂɊւ���ϐ��̃f�[�^���t�@�C������ǂݍ��ށ�
!
!��{�I�ȕϐ��̐���
!NPOIN  �ߓ_��
!NFPOIN ���R�ߓ_���iNFPOIN+1�`NPOIN�̐ߓ_�͑S�S���j
!NVFIX  �S���ߓ_��
!NELEM  ���v�f��
!NVARM  ���v�f�̎�ސ�
!NPRET  ���v�f�̏������̓f�[�^�A=0�F�Ȃ��A=1�F����
!NELEC  �P�[�u���v�f��
!NVARC  �P�[�u���v�f�̎�ސ�
!NPRSTR �P�[�u���̏������̓f�[�^�A=0�Ȃ��A=1����
!NELEF  �Ȃ��v�f��
!NVARF  �Ȃ��v�f�̎�ސ�
!NPFRM  �Ȃ��v�f�̏������̓f�[�^�A=0�Ȃ��A=1����
!NTIMES �׏d�����@�ŉ�͂���׏d�P�[�X��
!NCONC  �׏d�l�𒼐ړ��͂���ߓ_�̐�

!LODFIX �׏d�x�N�g���̈����A=0 �ό`�Ǐ]�A=1�����l�ŌŒ�
!MAXCYL �j���[�g���@�̍ő�J��Ԃ���
!DISMAX �j���[�g���@�̎�������
!ERR    �ߓ_�s�ލ��͂̎�������
!NREPT  �s�ލ��͎����̂��߂̍ő�J��Ԃ��v�Z��

read(7,*)
read(7,*)text,NPOIN
read(7,*)text,NFPOIN
read(7,*)text,NVFIX
read(7,*)
read(7,*)text,NELEM
read(7,*)text,NVARM
read(7,*)text,NPRET
read(7,*)
read(7,*)text,NELEC
read(7,*)text,NVARC
read(7,*)text,NPRSTR
read(7,*)
read(7,*)text,NELEF
read(7,*)text,NVARF
read(7,*)text,NPFRM
read(7,*)
read(7,*)text,NTIMES
read(7,*)text,NCONC
read(7,*)
read(7,*)
read(7,*)text,LODFIX
read(7,*)text,MAXCYL
read(7,*)text,DISMAX
read(7,*)text,ERR
read(7,*)text,NREPT
read(7,*)text,NMONT

write(8,501)
write(8,502)NPOIN
write(8,503)NFPOIN
write(8,504)NVFIX
write(8,*)
write(8,505)NELEM
write(8,506)NVARM
write(8,507)NPRET
write(8,*)
write(8,508)NELEC
write(8,509)NVARC
write(8,510)NPRSTR
write(8,*)
write(8,511)NELEF
write(8,512)NVARF
write(8,513)NPFRM
write(8,*)
write(8,514)NTIMES
write(8,515)NCONC
write(8,516)
write(8,517)LODFIX
write(8,518)MAXCYL
write(8,519)DISMAX
write(8,520)ERR
write(8,521)NREPT
write(8,522)NMONT

501 format(' ****************** NONLINEAR ANALYSIS OF MEMBRANE STRUCTURE (SOL. BY F.E.M.) *****************,'&
           //,67X,'-- GEOMETRIC NONLINEARITY --',//,64X,'-- THE NEWTON-RAPHSON METHOD --',&
           //,47X,'-- WRINKLING ANALYSIS(STRESS TRANSFER METHOD) --',&
           ///,'���̓f�[�^�̐���ϐ�')
502 format(9X,'NPOIN ',I10,5X,'�ߓ_��')
503 format(9X,'NFPOIN',I10,5X,'���R�ߓ_���iNFPOIN+1�`NPOIN�̐ߓ_�͑S�S���j')
504 format(9X,'NVFIX ',I10,5X,'�S���ߓ_��')
505 format(9X,'NELEM ',I10,5X,'���v�f��')
506 format(9X,'NVARM ',I10,5X,'���v�f�̎�ސ�')
507 format(9X,'NPRET ',I10,5X,'���v�f�̏������̓f�[�^�A=0�F�Ȃ��A=1�F����')
508 format(9X,'NELEC ',I10,5X,'�P�[�u���v�f��')
509 format(9X,'NVARC ',I10,5X,'�P�[�u���v�f�̎�ސ�')
510 format(9X,'NPRSTR',I10,5X,'�P�[�u���̏������̓f�[�^�A=0�Ȃ��A=1����')
511 format(9X,'NELEF ',I10,5X,'�Ȃ��v�f��')
512 format(9X,'NVARF ',I10,5X,'�Ȃ��v�f�̎�ސ�')
513 format(9X,'NPFRM ',I10,5X,'�Ȃ��v�f�̏������̓f�[�^�A=0�Ȃ��A=1����')
514 format(9X,'NTIMES',I10,5X,'�׏d�����@�ŉ�͂���׏d�P�[�X��')
515 format(9X,'NCONC ',I10,5X,'���ڐߓ_�׏d�f�[�^�A=0�Ȃ��A=1����')
516 format(/,'�v���O���������̐���ϐ�')
517 format(9X,'LODFIX',I10,5X,'�׏d�x�N�g���̈����A=0 �ό`�Ǐ]�A=1�����l�ŌŒ�')
518 format(9X,'MAXCYL',I10,5X,'�j���[�g���@�̍ő�J��Ԃ���')
519 format(9X,'DISMAX',F10.4,5X,'�j���[�g���@�̎�������')
520 format(9X,'ERR   ',F10.4,5X,'�ߓ_�s�ލ��͂̎�������')
521 format(9X,'NREPT ',I10,5X,'�s�ލ��͎����̂��߂̍ő�J��Ԃ��v�Z��')
522 format(9X,'NMONT ',I10,5X,'���W�l�����j�^�����O����iWRITE(10,*)�ŏo�͂���j�ߓ_�̔ԍ�')

nflag=0
if(NTIMES.eq.0)write(8,*)'NTIMES��1�ȏ�ɂ��Ă�������'
if(NTIMES.eq.0)nflag=1
if(nflag.eq.1)call echo

!�z��T�C�Y�Ɋ֌W����p�����[�^
NDOFN = 6                          !�ߓ_���R�x��
MSTRE = 3                          !���v�f�̉��͐����̐�

!��̓I�v�V����
NLOAD  = 0
PRETEN = 0.0                       !NPRET=0�̂Ƃ��A�����}�g���N�X���ق�����邽�߂̔����ȏ�������
DDV    = 1.0                       !���v�f�����N�����O���̍����ጸ��
DV     = 1.0                       !�P�[�u�����ɂ񂾏ꍇ�ɍ�����(1/DV)�{����
NPROB  = NMONT                     !���̔ԍ��̐ߓ_���W�f�[�^���t�@�C���ɏo�͂���WRITE(11,*)

!�B���@�\�i���ʂȉ�͂��s���Ƃ��ɐݒ肷��p�����[�^�j
MINCR  = 1                         !�׏d�����@�i=1�j�A�ʒ������@�i=2�j
MCTEMP = 0                         !�P�[�u���̉��x���͉�͐�p�̐ݒ�
MPOND  = 0                         !MPOND�F�|���f�B���O��͂��s��Ȃ��i=0�j�A�s���i=1�j
MEIGN  = 0                         !�ŗL�l�̌v�Z�i=0�F���Ȃ��A=1�F����j
NONISO = 1                         !���ޗ��̓��� =0�������A=1�ٕ����A=2�ٕ������@�ە����Ɩ��v�fij�������X��
NANGLE = 0                         !NONISO=2�̏ꍇ�ɁA�Y�����閌�v�f�̐�
NPRINT = 2                         !�o�͂̐���i=0	�v�Z�ߒ����ׂāA=1	�v�Z�ߒ��̍ő�ψʁA�ő�s�ލ��́A��������A�׏d���o�́A=2	�v�Z�ߒ��̍ő�ψʁA�ő�s�ލ��͂��o�́A=3	�v�Z���ʂ̂݁j

!�ł���΍폜�������ϐ�
NSTEP  = 1

!
!*** �����I�z��̊��蓖�ā�

allocate (PMGNF(NDOFN,NTIMES))     !�׏d�{��
allocate (PPRES(NTIMES))           !����[N/m2]
allocate (WWIND(NTIMES))           !���x��[N/m2]
allocate (SSNOW(NTIMES))           !��׏d[N/m2]
allocate (F(NPOIN*NDOFN))          !���׏d�A��׏d�x�N�g��
allocate (FF(NPOIN*NDOFN))         !�����̓x�N�g��
allocate (FFF(NPOIN*NDOFN))        !�����׏d�x�N�g��
allocate (FINIT(NPOIN*NDOFN))      !�Œ�׏d�x�N�g��
allocate (FCOE(NPOIN*NDOFN))       !���ڐߓ_�׏d�̉׏d�W��
allocate (Q(NPOIN*NDOFN))          !�c���̓x�N�g��
allocate (PSUM(NDOFN))             !���ڐߓ_�O�͂̊e�����̍��v

allocate (COORD(NPOIN,3))          !�ߓ_���W
allocate (COOR0(NPOIN,3))          !�����`�󎞂̐ߓ_���W
allocate (COOR1(NPOIN,3))          !�j���[�g���@�ɂ��J�Ԃ��v�Z�O�� �ߓ_���W
allocate (DIS(NPOIN*NDOFN))        !�ψ�
allocate (TTDIS(NPOIN*NDOFN))      !�����`��ɑ΂���ψ�

allocate (IFFIX(NPOIN*NDOFN))      !�S���Ɋւ�����i0:���R�A1:�S��)
allocate (NOFIX(NVFIX))            !�S���ߓ_�̔ԍ�
allocate (PRESC(NVFIX,NDOFN))      !IVFIX�Ԗڂ̍S���ߓ_��IDOFN�Ԗڂ̎��R�x�̍S����
allocate (PRESC0(NVFIX,NDOFN))     !IVFIX�Ԗڂ̍S���ߓ_��IDOFN�Ԗڂ̎��R�x�̍S���ʁi�ݒ�l���i�[���邽�߂̔z��j
allocate (IFPRE(NVFIX))            !�S����Ԃ�\���ϐ��i��@101:x,z�����S���Ay�������R�j

allocate (NODEM(NELEM,4))          !���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
allocate (PROPM(NVARM,7,3))        !���v�f�̍ޗ������i��ޔԍ��A��������(1�`7)�A�i�KTri-Linear(1�`3)�j
allocate (STRSM(NELEM,MSTRE))      !���v�f�̒���
allocate (STRM1(NELEM,MSTRE))      !�����N�����O�������s���O�̖��v�f�̒���
allocate (NANG(NELEM))             !���v�fIJK��IJ�����ƍޗ��̎厲�������X���Ă���v�f�̔ԍ�
allocate (ANG(NELEM))              !I���v�fIJK��IJ�����ƍޗ��̎厲�����̂Ȃ��p�x
allocate (MEME1(NELEM))            !�剞�͕����w�̃����N�����O�����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
allocate (MEME2(NELEM))            !�剞�͕����x�̃����N�����O�����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
allocate (AREA(NELEM))             !���v�f�̏����`�󎞂̖ʐ�
allocate (WCF(NELEM))              !�����W��
allocate (SCF(NELEM))              !��׏d�W��

allocate (NODEC(NELEC,3))          !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
allocate (BL(NELEC))               !�P�[�u���̒���
allocate (CAL(NELEC))              !�P�[�u���̖��Ђ��ݒ���
allocate (STRSC(NELEC))            !�P�[�u���̒���
allocate (STRC1(NELEC))            !���݂̏������s���O�̃P�[�u���̒���
allocate (NEA(NELEC))              !�P�[�u���̂��݂̔����󋵁i=0�F�������Ă���A=1�F���Ȃ��j
allocate (PROPC(NVARC,3))          !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj

allocate (BETA(NELEF))             !�Ȃ����ނ̉�]�p
allocate (CALFR(NELEF))            !���ނ̒���
allocate (VALK(NELEF,2))           !�ޒ[1,2�̐ڍ������i��1.0E+20�F���ځA��0.0�F�s���ځA���̑��F�����ځj
allocate (NODEF(NELEF,3))          !�Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
allocate (NPIN(NELEF))             !���ޒ[���̐ڍ������i=0�F�����A=1�F�s�����A=2�F���s���A=3�F�s���s���j
allocate (PROPF(NVARF,10))         !�Ȃ��v�f�̍ޗ������i��ޔԍ��A��������(1�`10)�j
allocate (STRSF(NELEF,NDOFN*2))    !�Ȃ��v�fIELEF�̕��ޗ́i����W�n�ɂ����镔�ޒ[�׏d)
allocate (STRF1(NELEF,NDOFN*2))    !�Ȃ��v�fIELEF�̕��ޗ́i���ލ��W�n�ɂ����镔�ޒ[�׏d)

allocate (EE(3,3))                 !�e���萔�}�g���N�X�c
allocate (NORDR(NPOIN*NDOFN))      !�S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ����ĕ������鎞�̐V�����S�̎��R�x�ԍ�
allocate (NWDTH(NPOIN*NDOFN))      !�S�̍����}�g���N�X�̑�ITOTV �s�̑�P��[���v�f����Ίp�v�f�܂ł̌�
allocate (NWSUM(0:NPOIN*NDOFN))    !�S�̍����}�g���N�X�̑�ITOTV �s�̑�P��[���v�f����Ίp�v�f�܂ł̌�

! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>
allocate (ELOAD(NPOIN*NDOFN))      !�c����
allocate (TLOAD(NPOIN*NDOFN))      !���׏d
allocate (ADIS1(NPOIN*NDOFN))      !���t'[k,j1]�ATLOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
allocate (ADIS2(NPOIN*NDOFN))      !���t'[k,j2]�AELOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
allocate (ASDIS(NPOIN*NDOFN))      !�t[k,(j)]�A��IINCS�񑝕��ɂ��ψʗ�
allocate (RPRVS(NPOIN*NDOFN+1))    !��(IITER-1)�񔽕���̂��x�N�g��
allocate (D(NPOIN*NDOFN))          !SUB SKYLNE�̍�Ɨp�z��
allocate (T(NPOIN*NDOFN))          !SUB SKYLNE�̍�Ɨp�z��
allocate (NWK(NPOIN*NDOFN))        !SUB SKYLNE�̍�Ɨp�z��
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
allocate (CTEMP1(NTIMES))          !CTEMP1(LOADCREM)�F�P�[�u���ɗ^���鉷�x����[��]�AMCTEMP=1�i�P�[�u���̉��x���͉�͂��s���j�̏ꍇ
allocate (CTEMP2(2,NELEC))         !CTEMP2(1,IELEC)�F���x���͂�^����(=1.0)�A�^���Ȃ�(=0.0)�ACTEMP2(2,IELEC)�F���x�Ђ���

! *** ���C���v�b�g�f�[�^�t�@�C������ߓ_���W�A���v�f�A�׏d�Ɋւ���f�[�^���t�@�C�����V����Ǎ��ށ� ***

CALL INPUT(ANG,COORD,DDV,DISMAX,ERR,FCOE,IFPRE,LODFIX, &
           MAXCYL,MEIGN,MSTRE,NANG,NANGLE,NCONC,NDOFN,NELEM,NFPOIN,&
           NLOAD,NODEM,NOFIX,NONISO,NPOIN,NPRET,NPRINT,NPROB,NREPT,  &
           NSTEP,NTIMES,NVFIX,NVARM,PMGNF,PPRES,PRESC,   &
           PRETEN,PROPM,SCF,SSNOW,STRSM,WCF,WIND,WWIND,MINCR)
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>
      IF(MINCR.EQ.2)THEN
        NTIMES=500
        MAXCYL=1
        ARCLG=PMGNF(1,1)    !    ARCLG=0.25D00
        IF(ARCLG.EQ.0.0) STOP ' ARCLG=0 '
        TFACT=0.0D00
        CALL ZEROR2(PMGNF,NDOFN,NTIMES,NDOFN,NTIMES)
        CALL ZEROR1(PPRES,NTIMES,NTIMES)
        CALL ZEROR1(WWIND,NTIMES,NTIMES)
        CALL ZEROR1(SSNOW,NTIMES,NTIMES)
        CALL ZEROR1(ASDIS,MTOTV,NTOTV)
        CALL ZEROR1(TLOAD,MTOTV,NTOTV)
        WRITE(13,*) ' >>>>> ARCH-LENGTH METHOD >>>>> '
        WRITE(13,*) ' (NTIMES,MAXCYL,ARCLG) = (',NTIMES,MAXCYL,ARCLG,' )'
      ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! *** ���P�[�u���v�f�Ɋւ���f�[�^���t�@�C�����V����Ǎ��ށ� ***
!
IF(NELEC.NE.0)THEN
  CALL CABLE1(BL,CAL,COORD,DV,NODEC,NPOIN,NELEC,NPRSTR,NVARC,&
                  STRSC,STRC1,PROPC,MCTEMP,NTIMES,CTEMP1,CTEMP2)
ENDIF
!
! *** ���Ȃ��v�f�Ɋւ���f�[�^���t�@�C�����V����Ǎ��ށ� ***
!
IF(NELEF.NE.0)THEN
  CALL FRAME(1,CALFR,COORD,VALK,NPOIN,NDOFN,NELEF,NODEF,NPIN,&
             NPFRM,NVARF,PROPF,STRSF,STRF1,BETA,IFFIX,&
             MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
ENDIF
!
! *** �������`��̃f�[�^���t�@�C�����X�ɕۑ����遄 ***
!
 CALL WRTCOD(1,COORD,NODEC,NELEM,NELEC,NODEM,NPOIN,NELEF,NODEF)
!
! *** ���v�f�����}�� VECTOR SCRIPT �`���̃t�@�C���ŏo�͂��遄
!
LOADCREM=0
CALL WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                   NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,NPROB,NTIMES,&
                   STRSF,STRSM,STRSC,TTDIS)
!
! *** ���e�ϐ��̏����l�̐ݒ聄 ***
!
MTOTV=NPOIN*NDOFN
NTOTV=NPOIN*NDOFN
CALL ZEROR1(TTDIS,MTOTV,NTOTV)
CALL ZEROR1(F,MTOTV,NTOTV)
CALL ZEROR1(FF,MTOTV,NTOTV)
NFREE=0
!
! *** ��CMEM�^�̍S���ߓ_�����ASKYLNE�^�ɕύX���遄 ***
!
CALL FIXNO(IFFIX,IFPRE,NODEM,NODEC,MTOTV,NDOFN,NELEM,NELEC,NFPOIN,  &
           NFREE,NOFIX,NORDR,NPOIN,NTOTV,NVFIX,NELEF,NODEF)
!
! *** ���S�̍����}�g���N�X�̐����ƂP�����z�񉻂��ꂽ�X�J�C���C���}�g���N�X�̑Ή��������w�W�z����쐬���遄 ***
!
CALL SIZING(NODEM,NODEC,MSIZE,MTOTV,NELEM,NELEC,NDOFN,NFREE,NORDR,     &
            NPOIN,NTOTV,NWDTH,NWSUM,NELEF,NODEF)
!
! *** ���S�̍����}�g���N�X�i�X�J�C���C���}�g���N�X�j�̔z�񊄂蓖�ā� ***
!
allocate (TSTIF(MSIZE))            !MSIZE�F�X�J�C���C���}�g���N�X�̐����l
!
! ------------------------------------------------------------------------------------------------
! *** ���A���S���Y���̐���ϐ��� ***
!     LOADCREM : �׏d�����񐔁i�P�`NTIMES)
!     NNN      : �c���͏����̂��߂̌J�Ԃ��v�Z�񐔁i1�`NREPT+1)
!     INCREM   : �׏d�����񐔁i�P�`NSTEP)
!     NCYCLE   : �j���[�g������t�\���@�ɂ������v�Z�̌J�Ԃ��񐔁i1�`MAXCYL)
!     NCHECK   : =0 �������Ă��Ȃ���
!                =1 ����������
!     NFINAL   : =0 �v�Z�r��
!                =1 �v�Z�I�����O��,�ό`��̍��W���Q�Ƃ��ĕs�ލ��͂��v�Z���鎞�B
!     NUNBAL   : =0 �j���[�g������t�\�������v�Z��
!                =1 �j���[�g���B���t�\�������v�Z�I����C�����͂����߂鎞�B�������͂O�ŏ���̖����͂����߂鎞�B
!                =2 �j���[�g���B���t�\�������v�Z�I����C�����͂����߂鎞�B
! ------------------------------------------------------------------------------------------------
!
DO 998 LOADCREM=1,NTIMES  ! *** ������i�׏d����1�`NTIMES�j�̃��[�v�j ***
  IF(LOADCREM.GT.1)THEN
    NPRET=1
    NSTEP=1
  ENDIF
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>������
  IF(MINCR.EQ.2)CALL ZEROR1(ASDIS,MTOTV,NTOTV)
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
! *** ���׏d�̐ݒ聄 ***
!
  CALL SETLOAD(LOADCREM,NPOIN,F,PRES,WIND,SNOW,LODFIX,NTIMES,NDOFN,FCOE,PMGNF,PPRES,WWIND,SSNOW)
!
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>������
! *** �����x�ω��ɂ��P�[�u���̂Ђ��ݑ����̐ݒ聄 ***
!
  IF(MCTEMP.EQ.1)THEN   ! MCTEMP=1�F�P�[�u���̉��x���͕ω����l������ꍇ�A��������ꍇ��=0
    DO IELEC=1,NELEC
      IF(CTEMP2(1,IELEC).NE.0.0)THEN
        CTEMP2(2,IELEC)=CTEMP1(LOADCREM)*12.0E-06  ! CTEMP1(LOADCREM)�F���x����[��]�A12.0E-06�F�M�c���W��[/ ��]
      ENDIF
    ENDDO
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
! *** ����������e��VOL0�̌v�Z�i�����ϓ��i�{�C�����j���l������ꍇ�ɕK�v�Ȍv�Z�j�� ***
!
   PRES0=PRES
   APRESOUT=101325.0/9.81   ! ��C����10130kg/m2�Ɖ��肷��
   IF(LODFIX.EQ.-1)THEN
     CALL VOLUME(VOL0,COORD,NODEC,NODEM,NPOIN,NELEC,NELEM)
   ENDIF
!
! *** ���׏d�p�����[�^�̕\��> ***
!
   WRITE(6,178)LOADCREM,NTIMES,(PMGNF(IDOFN,LOADCREM),IDOFN=1,NDOFN),PRES,WIND,SNOW
   WRITE(8,178)LOADCREM,NTIMES,(PMGNF(IDOFN,LOADCREM),IDOFN=1,NDOFN),PRES,WIND,SNOW
   IF(MINCR.EQ.1)THEN
     WRITE(11,178)LOADCREM,NTIMES,(PMGNF(IDOFN,LOADCREM),IDOFN=1,NDOFN),PRES,WIND,SNOW
   ENDIF
  178 FORMAT(/,'   LOADCREM (',I3,'/',I3,' )',                          &
               '  Px -dir:',F9.2,'  Py -dir:',F9.2,'  Pz -dir:',F9.2,   &
               '  Pxx-dir:',F9.2,'  Pyy-dir:',F9.2,'  Pzz-dir:',F9.2,   &
         /,22X,'  Pres   :',F9.2,'  Wind   :',F9.2,'  Snow   :',F9.2)
!
! *** ���w�x�y�����̑��׏d�̌v�Z�� ***
!
  CALL ZEROR1(PSUM,3,3)
  WRITE(8,176)
  176 FORMAT(//' (NODAL LOADINGS BY DIRECT INPUT) [N]',/,'  NOD',9X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD', &
                8X,'XX-LOAD',8X,'YY-LOAD',8X,'ZZ-LOAD')
  DO IPOIN=1,NPOIN
    NZERO=0
    ITOTV1=(IPOIN-1)*NDOFN
    DO IDOFN=1,NDOFN
      ITOTV=ITOTV1+IDOFN
      PSUM(IDOFN)=PSUM(IDOFN)+F(ITOTV)
      IF(F(ITOTV).NE.0.0)NZERO=1
    ENDDO
    IF(NZERO.EQ.1)WRITE(8,131)IPOIN,(F(IITTV),IITTV=ITOTV1+1,ITOTV)
  ENDDO
  WRITE(8,177) (PSUM(IDOFN),IDOFN=1,NDOFN)
  131 FORMAT(I5,6F15.4)
  177 FORMAT(1X,95('-'),/,' SUM ',6F15.4)
!
! *** ���e�ϐ��̏����l�̐ݒ聄 ***
!
  IF(NPRET.EQ.0)THEN
    CALL ZEROR2(STRSM,NELEM,MSTRE,NELEM,MSTRE)
  ELSE
    PRETEN=0.0
  ENDIF
!
  DO IELEC=1,NELEC
    NEA(IELEC)=1
  ENDDO
!
  DO IELEM=1,NELEM
    MEME1(IELEM)=1
    MEME2(IELEM)=1
    DO ISTRE=1,3
      STRM1(IELEM,ISTRE)=STRSM(IELEM,ISTRE)
    ENDDO
  ENDDO
!
  DO ITOTV=1,NTOTV
    FINIT(ITOTV)= F(ITOTV)
  ENDDO
!
  DO IPOIN=1,NPOIN
    DO IDIR=1,3
      COOR0(IPOIN,IDIR)=COORD(IPOIN,IDIR)
    ENDDO
  ENDDO
!
  IF(LOADCREM.EQ.1)THEN
    DO NE=1,NELEM
      I=NODEM(NE,1)
      J=NODEM(NE,2)
      K=NODEM(NE,3)
      A=(COORD(J,2)-COORD(I,2))*(COORD(K,3)-COORD(I,3))-(COORD(K,2)-COORD(I,2))*(COORD(J,3)-COORD(I,3))
      B=(COORD(J,3)-COORD(I,3))*(COORD(K,1)-COORD(I,1))-(COORD(K,3)-COORD(I,3))*(COORD(J,1)-COORD(I,1))
      C=(COORD(J,1)-COORD(I,1))*(COORD(K,2)-COORD(I,2))-(COORD(K,1)-COORD(I,1))*(COORD(J,2)-COORD(I,2))
      AREA(NE)=SQRT(A*A+B*B+C*C)
    ENDDO
  ENDIF
!
  NFINAL = 0
  IF(NREPT.EQ.0) NFINAL=1
!
  DO 3000 NNN=1,NREPT+2  ! *** ������i�s�ލ��͏����̂��߂̌J�Ԃ��v�Z�P�`NPREPT+2�̃��[�v�j***
!   NNN=NREPT+1�F�����Ɏ���Ȃ������ꍇ�ɍŏI�s�ލ��͂̌v�Z���s��
!   NNN=NREPT+2�F�����Ɏ���Ȃ������ꍇ�ɃP�[�u�����͂̌v�Z���s��
!
    CALL ZEROR1(F,MTOTV,NTOTV)
!
    IF(LODFIX.NE.2)THEN
!
! >>>>>>>>>> �|���f�B���O��͂��s���ꍇ�̐ݒ�>>>>>>������
! LOADCREM=2�`NTIMES�FWWIND(LOADCREM)�𒙐����AWIND�𐅈�[N/m2]�AWCF(IELEM)���O�p�`�v�f�̏d�S�̂y���W�ɉ����������W���Ɛݒ肵�Čv�Z����B
    IF(MPOND.EQ.1)THEN
      IF(LOADCREM.GE.2)THEN
        ZMAX=COORD(1,3)              ! ZMAX�F�y���W�̍ő�l
        ZMIN=COORD(1,3)              ! ZMIN�F�y���W�̍ŏ��l
        DO IPOIN=2,NPOIN
          IF(COORD(IPOIN,3).GT.ZMAX)ZMAX=COORD(IPOIN,3)
          IF(COORD(IPOIN,3).LT.ZMIN)ZMIN=COORD(IPOIN,3)
        ENDDO
        ZHGHT=ZMAX-ZMIN             ! ZHGHT�F�y���W�̍ō��_����Œ�_�܂ł̍���
        ZSURF=ZMIN+ZHGHT*WWIND(LOADCREM)  ! ZSURF�F���ʂ̂y���W
        WIND=98000.0          ! WIND�F���ɂ�鈳��98000[N/m2]
        DO IELEM=1,NELEM
          NFLAG=0
          VALZ=0.0            ! VALZ�F�O�p�`�v�f�̏d�S�̂y���W���v�Z����B
          DO INODE=1,3
            VALZ=VALZ+COORD(NODEM(IELEM,INODE),3)
            IF(COORD(NODEM(IELEM,INODE),3).GT.ZSURF)NFLAG=NFLAG+1
          ENDDO
          WCF(IELEM)=(ZSURF-VALZ/3.0)                     ! WCF(IELEM)�F�����W���i�����[�j
          IF(WCF(IELEM).LT.0.0)WCF(IELEM)=0.0             ! �d�S�̂y���W�����ʂ���ɂ���ꍇ
          IF(NFLAG.EQ.1)WCF(IELEM)=WCF(IELEM)*2.0/3.0     ! �R�ߓ_�̂����P�_�����ʂ���ɂ���ꍇ
          IF(NFLAG.EQ.2)WCF(IELEM)=WCF(IELEM)*1.0/3.0     ! �R�ߓ_�̂����Q�_�����ʂ���ɂ���ꍇ
          IF(NFLAG.EQ.3)WCF(IELEM)=0.0                    ! �R�ߓ_���ׂĂ����ʂ���ɂ���ꍇ
        ENDDO
        IF(NNN.EQ.1)THEN
          WRITE(8,198)WIND,ZSURF,ZMIN,(IELEM,WCF(IELEM),IELEM=1,NELEM)
        ELSE
          WRITE(8,198)WIND,ZSURF,ZMIN
        ENDIF
        198 FORMAT(//,' (WATER PRESSURE CO-EFFICIENTS AT ELEMENTS)',/,'   PRESSURE = ',F10.4,' [N/m2/m]'&
                   '  SUREFACE LEVEL =',F10.4,'  BOTTOM LEVEL =',F10.4,/,(5(I5,F8.4,7X)))
      ENDIF
    ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
! *** ���ό`��̓���e��VOL1�̌v�Z�i�����ϓ��i�{�C�����j���l������ꍇ�ɕK�v�Ȍv�Z�j�� ***
!
      IF(NNN.GE.2.AND.LODFIX.EQ.-1.AND. &
        (WIND.NE.0.0.OR.SNOW.NE.0.0.OR.PSUM(1).NE.0.0.OR.PSUM(2).NE.0.0.OR.PSUM(3).NE.0.0))THEN
        CALL VOLUME(VOL1,COORD,NODEC,NODEM,NPOIN,NELEC,NELEM)
        IF(VOL1.EQ.0.0)STOP ' VOL1=0.0 (STOP AT LINE 327) '
!       PRES=PRES0*VOL0/VOL1
        VEXTRA=VOL0*200.0   ! ��͑Ώۂ̍\������VOL0��200�{�̗e�ς̍��̂ɘA������Ă���Ɖ���B
        PRES=(PRES0+APRESOUT)*(VOL0+VEXTRA)/(VOL1+VEXTRA)-APRESOUT
        WRITE(6,6001)PRES,VOL0,VOL1,VOL0/VOL1
        WRITE(11,6001)PRES,VOL0,VOL1,VOL0/VOL1
        6001 FORMAT(' PRES=',G10.3,' VOL0=',G10.3,' VOL1=',G10.3,' VOL0/VOL1=',G10.3)
        IF(ABS(VOL0/VOL1-1).GT.PRES0/APRESOUT)THEN
!         PRES=10.0
        ELSEIF(ABS(PRES).GT.100.0)THEN
!         PRES=PRES/10.0
!         PRES=(PRES0+APRESOUT)*VOL0/VOL1*1.0E-3-APRESOUT
        ENDIF
      ENDIF
!
! *** ���׏d�x�N�g���̌v�Z�� ***
!
      CALL LOAD(CAL,COORD,F,FF,FINIT,NODEC,NVARC,NDOFN,NELEM,NELEC,NFINAL,&
                NNN,NODEM,NPOIN,NPRINT,PRES,SCF,SNOW,WCF,WIND,NELEF,NODEF,&
                CALFR,NVARF,MPOND,PROPC,PROPF,PROPM,NVARM,LOADCREM,WWIND,NTIMES)
!
      IF(LODFIX.EQ.1)THEN
        LODFIX=2  ! LODFIX=1�i�׏d�x�N�g���Œ�j�̏ꍇ�ɂ́A������LODFIX=2�ƂȂ莟��ȍ~�׏d�x�N�g�����X�V����Ȃ��Ȃ�
        DO ITOTV=1,NTOTV
          FINIT(ITOTV)=FINIT(ITOTV)+F(ITOTV)
          F(ITOTV)=0.0
          WCF(ITOTV)=FF(ITOTV)
        ENDDO
      ENDIF
    ENDIF
!
    IF(LODFIX.EQ.1.OR.LODFIX.EQ.2)THEN
      DO ITOTV=1,NTOTV
        FF(ITOTV)=WCF(ITOTV)
      ENDDO
    ENDIF
!
    IF(NNN.EQ.1) WRITE(8,955)
    955 FORMAT(///,' *************** PROCESS OF COMPUTATION **************',//)
!
    IF(NNN.GE.2) NSTEP=1  ! �c���͎����̂��߂̔�����NNN=2��ڈȍ~�́A�׏d�������s��Ȃ��B
!
    DO ITOTV=1,NTOTV
      FFF(ITOTV)=F(ITOTV)+FINIT(ITOTV)
      IF(NLOAD.EQ.0)FFF(ITOTV)=(FFF(ITOTV)+FF(ITOTV))/FLOAT(NSTEP)
      IF(NLOAD.GE.1)FFF(ITOTV)=FFF(ITOTV)/FLOAT(NSTEP)
    ENDDO
!
    DO 2000 INCREM=1,NSTEP  ! *** ������i�׏d�����P�`NSTEP�̃��[�v�j***
!
      CALL ZEROR1(DIS,MTOTV,NTOTV)
!
      IF(NFINAL.EQ.0)THEN  ! NFINAL=0�F�v�Z�r���̏ꍇ
        DO IPOIN=1,NPOIN
          DO IDIR=1,3
            COOR1(IPOIN,IDIR)=COORD(IPOIN,IDIR)
          ENDDO
        ENDDO
      ENDIF
!
! *** ���׏d�����̐ݒ聄 ***
!
      FINC=FLOAT(INCREM-1)
      IF(NNN.GE.2.OR.INCREM.EQ.1)FINC=1.0
      DO ITOTV=1,NTOTV
        IF(NLOAD.EQ.0)THEN
          FFF(ITOTV) = FFF(ITOTV)*FLOAT(INCREM)/FINC
        ELSEIF(NLOAD.GE.1)THEN
          IF(INCREM.GE.2)FFF(ITOTV)=FFF(ITOTV)-FF(ITOTV)
          FFF(ITOTV)=FFF(ITOTV)*FLOAT(INCREM)/FINC+FF(ITOTV)
        ENDIF
      ENDDO
!
! *** ���P�[�u�����͂̌v�Z�� ***
!
      IF(NELEC.NE.0.AND.(NNN.GT.1.OR.INCREM.GT.1))THEN
        CALL CABLE3(0,BL,CAL,NELEC,NVARC,NODEC,STRSC,STRC1,PROPC,MCTEMP,CTEMP2)
      ENDIF
!
      NUNBAL=0
      IF(NFINAL.EQ.1) NUNBAL=2
!
      DO 1000 NCYCLE=1,MAXCYL  ! ***������i�j���[�g���@�̔����v�Z�̃��[�v�j***
!
        CALL ZEROR1(Q,MTOTV,NTOTV)
        CALL ZEROR1(TSTIF,MSIZE,MSIZE)

        IF(NSTEP.GT.1)THEN
          IF(LOADCREM.EQ.1.AND.NNN.EQ.1.AND.NCYCLE.EQ.1)THEN
            IF(INCREM.EQ.1)PRESC0=PRESC/NSTEP
            PRESC=PRESC0
            !READ(7,908)(NOFIX(IVFIX),IFPRE(IVFIX),(PRESC(IVFIX,IDOFN),IDOFN=1,NDOFN),IVFIX=1,NVFIX)
          ENDIF
        ENDIF
!
! **** �����v�f�̍����}�g���N�X���쐬���遄 ***
!
        460 CONTINUE
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>������
        IF(MINCR.EQ.2)THEN
          DO ITOTV=1,NTOTV
            FFF(ITOTV)=TLOAD(ITOTV)
          ENDDO
        ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
        IF(NELEM.GT.0)THEN
          CALL MEMBRN(ANG,AREA,COOR1,COORD,DDV,EE,IFFIX,INCREM,MEME1,&
                      MEME2,NPOIN,MSIZE,MSTRE,MTOTV,NANG,NANGLE,NCYCLE,&
                      NDOFN,NELEM,NFINAL,NNN,NODEM,NONISO,NORDR,NPROB,&
                      NUNBAL,NVARM,NWSUM,PRETEN,PROPM,Q,STRM1,STRSM,TSTIF)
        ENDIF
!
! *** ���P�[�u���v�f�̍����}�g���N�X���쐬���遄 ***
!
        IF(NELEC.NE.0)THEN
          CALL CABLE2(BL,CAL,COORD,COOR1,DV,IFFIX,INCREM,NODEC,NVARC,&
                      MSIZE,MTOTV,NCYCLE,NDOFN,NEA,NELEC,NFPOIN,NNN,NORDR,&
                      NPOIN,NPRET,NPRSTR,NTOTV,NUNBAL,NWSUM,Q,&
                      TSTIF,STRC1,PROPC,MCTEMP,CTEMP2)
        ENDIF
!
! *** ���Ȃ��v�f�̍����}�g���N�X���쐬���遄 ***
!
        IF(NELEF.GT.0)THEN
          CALL FRAME(2,CALFR,COORD,VALK,NPOIN,NDOFN,NELEF,NODEF,&
                      NPIN,NPFRM,NVARF,PROPF,STRSF,STRF1,BETA,&
                      IFFIX,MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
        ENDIF
!
!
        IF(NUNBAL.NE.0)GOTO 730  ! ***�iNUNBAL��0�Ȃ�΁C�j���[�g���@�̔����v�Z�̃��[�v���甲����j***
!
! *** ���S�̍����}�g���N�XTSTIF�̌ŗL�l���v�Z���遄 ***
!
        IF(MEIGN.EQ.1)THEN
          CALL EIGEN(EIGMN,MSIZE,MTOTV,NFREE,NORDR,NPOIN,NWDTH,NWSUM,FFF,TSTIF)
        ENDIF
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>������
        IF(MINCR.EQ.2)THEN            
          IDOFN=3
          ITOTV=(NPROB-1)*NDOFN+IDOFN
          WRITE(8,8002)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,Q(ITOTV),FFF(ITOTV)
 8002     FORMAT(  'LOADCREM      NNN  NCYCLE   NPROB   IDOFN    TFACTT       Q(ITOTV)     FFF(ITOTV)',/,5I8,4G15.4)
        ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
        DO ITOTV=1,NTOTV
          Q(ITOTV)=FFF(ITOTV)-Q(ITOTV)
          F(ITOTV)=Q(ITOTV)
        ENDDO
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>������
            IF(MINCR.EQ.2)THEN
              DO ITOTV=1,NTOTV
                ELOAD(ITOTV)=Q(ITOTV)
              ENDDO
            ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
        CALL ZEROR1(DIS,MTOTV,NTOTV)
!
! *** ����������C������NCHECK=1�ɂ�NFINAL=1�Ƃ��čŌ�̌J�Ԃ��v�Z�֐i�ށ� ***
!
        IF(LOADCREM.EQ.101)NPAUSE=1
        IF(NCYCLE.EQ.1.AND.NNN.GT.1)THEN
          CALL CONVER(ERR,IFFIX,NPOIN,NCHECK,NDOFN,NFPOIN,Q)
          IF(NCHECK.EQ.1)THEN
            NFINAL=1
            GOTO 3000
          ENDIF
        ENDIF
!
!
        IF(MINCR.EQ.1)THEN  ! ������ �׏d�����@�̏ꍇ�̌v�Z ������
!
! *** ���W���}�g���b�N�X����щE�Ӄx�N�g���ɑ΂��āA�􉽊w�I�ȏ������������遄 ***
!
          CALL BCONUU(F,IFFIX,MSIZE,MTOTV,NDOFN,NFREE,NOFIX,NORDR,NTOTV,NVFIX,NWSUM,PRESC,TSTIF)
!
! *** ���X�J�C���C���@�ɂ��A���ꎟ�������̉������߂遄 ***
!
          CALL SKYLNE(1,TSTIF,F,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
!
! *** �������ψʃx�N�g��DIS�̎Z�聄 ***
!
          CALL RARNG(DIS,F,MTOTV,NORDR,NTOTV)
!
        ELSEIF(MINCR.EQ.2)THEN ! ������ �ʒ������@�̏ꍇ�̌v�Z ������  
!
! *** �ʒ������@(ARCH LENGTH METHOD)�̓K�p(INDEX=1�F�A���ꎟ�������̉�@(SKYLNE)�̑O�̏���)
!
          CALL ARCLM(1,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,IFFIX,NNN,LOADCREM,&
                     MTOTV,NDOFN,NOFIX,NORDR,NTOTV,NVFIX,PRESC,RAMDA,FCOE,RPRVS,RSCLE,&
                     TTDIS,TFACT,TLOAD)
!
! *** �S�̍����}�g���N�XTSTIF�̌ŗL�l���v�Z����
!
          CALL EIGEN(EIGN1,MSIZE,MTOTV,NFREE,NORDR,NPOIN,NWDTH,NWSUM,ADIS1,TSTIF)
!
! *** �X�J�C���C���@�ɂ��A���ꎟ�������̉������߂�
!
          CALL SKYLNE(1,TSTIF,ADIS1,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
          IF(NNN.GT.1)THEN
            CALL SKYLNE(2,TSTIF,ADIS2,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
          ENDIF
!
! *** �ʒ������@(ARCH LENGTH METHOD)�̓K�p(INDEX=2�F�A���ꎟ�������̉�@(SKYLNE)�̌�̏���)
!
          CALL ARCLM(2,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,IFFIX,&
                     NNN,LOADCREM,MTOTV,NDOFN,NOFIX,NORDR,NTOTV,NVFIX,&
                     PRESC,RAMDA,FCOE,RPRVS,RSCLE,TTDIS,TFACT,TLOAD)
!
          DO ITOTV=1,NTOTV
            DIS(ITOTV)=ADIS1(ITOTV)
          ENDDO
!
        ENDIF  ! ������ �����܂ŁiMINCR�ɂ���������j������ 
!
!
! *** ���ߓ_���W�̍X�V�� ***
!
        DO IPOIN=1,NFPOIN
          DO IDIR=1,3
            ITOTV=(IPOIN-1)*NDOFN+IDIR
            COORD(IPOIN,IDIR)=COORD(IPOIN,IDIR)+DIS(ITOTV)
          ENDDO
        ENDDO
!
! *** �������`��ɑ΂���ψʃx�N�g��TTDIS�̎Z�聄 ***
!
        DO ITOTV=1,NTOTV
          TTDIS(ITOTV)=TTDIS(ITOTV)+DIS(ITOTV)
        ENDDO
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>������
        IF(MINCR.EQ.2)THEN
          IDOFN=3
          ITOTV=(NPROB-1)*NDOFN+IDOFN
          WRITE(8,8000)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,DRAMD,TLOAD(ITOTV),ELOAD(ITOTV),   &
                       DIS(ITOTV),TTDIS(ITOTV),ASDIS(ITOTV),Qmax,EIGMN
          IF(LOADCREM.EQ.1.AND.NNN.EQ.1)WRITE(11,8000)
          IF(NNN.EQ.1)  &
          WRITE(11,8000)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,DRAMD,TLOAD(ITOTV),ELOAD(ITOTV),  &
                        DIS(ITOTV),TTDIS(ITOTV),ASDIS(ITOTV),Qmax,EIGMN
 8000     FORMAT(/ ,'LOADCREM      NNN  NCYCLE   NPROB   IDOFN     TFACT     DRAMD        TLOAD(ITOTV)',  &
                 '   ELOAD(ITOTV)     DIS(ITOTV)   TTDIS(ITOTV)   ASDIS(ITOTV)            Qmax          EIGMN', &
                 /,5I8,9G15.4)
        ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
! *** ��NCYCLE=INCREM=NNN=1�̎��C�����ψʗʂ��[���ɂ��遄 ***
!
        IF(NCYCLE.EQ.1.AND.NNN.EQ.1)THEN
          CALL ZEROR2(PRESC,NVFIX,NDOFN,NVFIX,NDOFN)
        ENDIF
!
! *** ���v�Z�ߒ����t�@�C�����W�ɕۑ����遄 ***
!
        CALL PROC(0,DIS,EIGMN,IFFIX,INCREM,LOADCREM,NPOIN,NCYCLE,NDOFN,NFPOIN,&
                  NNN,NPRINT,NPROB,NTIMES,NTOTV,PRES,Q,COORD,QMAX,MINCR)
!
! *** �������ψʗʂɂ��������聄 ***
!
        NDCHK=0
        DO ITOTV=1,NTOTV
          IF(IFFIX(ITOTV).EQ.0.AND.ABS(DIS(ITOTV)).GT.DISMAX)THEN
            NDCHK=1
            GOTO 1000
          ENDIF
        ENDDO
!
! *** ��NDCHK=0(����)�܂���NCYCLE��MAXCYL�̏ꍇ�CNUNBAL=1or2�Ƃ��čs�ԍ�460��DO���[�v�֖߂遄 ***
!
! *** ���v�Z�ߒ��̏o�́�
!
        IF(NDCHK.EQ.0.OR.NCYCLE.GE.MAXCYL)THEN
          CALL PROC(1,DIS,EIGMN,IFFIX,INCREM,LOADCREM,        &
                    NPOIN,NCYCLE,NDOFN,NFPOIN,NNN,NPRINT,NPROB,     &
                    NTIMES,NTOTV,PRES,Q,COORD,QMAX,MINCR)
          PRETEN=0.0
          CALL ZEROR1(Q,MTOTV,NTOTV)
          NUNBAL=1
          IF(NPRET.NE.0.OR.INCREM.GE.2.0.OR.NNN.GE.2) NUNBAL=2
          GO TO 460
        ENDIF
!
      1000 CONTINUE ! *** �����܂Łi�j���[�g���@�̔����v�Z�̃��[�v�j***
!
      730 IF(NFINAL.EQ.1) GO TO 3100  ! ***�iNFINAL=1�Ȃ�΁C�v�Z���I�������ʂ��o�͂���j***
!
! *** ���c���͂p�̌v�Z��
!
      DO ITOTV=1,NTOTV
        Q(ITOTV)=Q(ITOTV)-FFF(ITOTV)
      ENDDO
!
      IF(NPRINT.LE.1)THEN
        WRITE(8,750)
        750 FORMAT(' NODE',3X,'X',12X,'Y',12X,'Z',15X,'DX',11X,'DY',    &
                   11X,'DZ',11X,'DXX',10X,'DYY',10X,'DZZ',13X,'QX',     &
                   11X,'QY',11X,'QZ',11X,'MX',11X,'MY',11X,'MZ')
        DO IPOIN=1,NPOIN
          ITOTV1=(IPOIN-1)*NDOFN+1
          ITOTV2=IPOIN*NDOFN
          WRITE(8,780)IPOIN,(COORD(IPOIN,IDIR),IDIR=1,3),(DIS(ITOTV),ITOTV=ITOTV1,ITOTV2),(Q(ITOTV),ITOTV=ITOTV1,ITOTV2)
          780 FORMAT(I5,3E13.5,2(3X,6E13.5))
        ENDDO
      ENDIF
!
    2000 CONTINUE ! *** �����܂Łi�׏d�����̃��[�v�j***
!
    ! �j���[�g���@�̔����v�Z�Ŏ������Ȃ������ꍇ�ɂ́ALASTCYL=NCYCLE,NFINAL=1�Ƃ���B 
    IF(NFINAL.NE.1.AND.NNN.EQ.NREPT+1)THEN
      LASTCYL=NCYCLE
      NFINAL=1
    ENDIF
!
  3000 CONTINUE ! *** �����܂Łi�c���͏����̂��߂̌J�Ԃ��v�Z�̃��[�v�j***
!
  3100 CONTINUE
!
  DO NF = 9, 10
    WRITE(NF,178)LOADCREM,NTIMES,(PMGNF(IDOFN,LOADCREM),IDOFN=1,NDOFN),PRES,WIND,SNOW
  ENDDO
!
! >>>>>>>>>> �ʒ������@��p�̐ݒ�>>>>>>>>>>>>>>>>>>��������
  IF(MINCR.EQ.2)THEN
    IDOFN=3
    ITOTV=(NPROB-1)*NDOFN+IDOFN
    IF(LOADCREM.EQ.1)WRITE(13,8000)
    WRITE(13,8000)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,DRAMD,TLOAD(ITOTV),ELOAD(ITOTV),DIS(ITOTV),TTDIS(ITOTV),ASDIS(ITOTV),Qmax,EIGMN
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<�����܂�
!
  IF(LODFIX.EQ.-1)WRITE(8,8100) PRES0,PRES,PRES-PRES0,VOL0,VOL1,VOL1-VOL0
  8100 FORMAT(//,15X,'     INITIAL     PRESENT     DIFFER.',/,5X,' PRESSURE:',3F12.3,' [N/m2]',/,5X,' VOLUME:  ',3F12.3,' [m*m*m] ')
!
! *** ���v�Z���ʂ̏o�́�
!
  CALL RESLT(BL,CAL,COOR0,COORD,FF,FFF,IFFIX,MSTRE,NDOFN,NELEM,&
             NELEC,NFINAL,NFPOIN,NPOIN,NVARC,NVARM,NODEC,NODEM,Q,STRSM,STRSC,TTDIS,&
             STRC1,PROPC,PROPM,NELEF,STRSF,MCTEMP,CTEMP2,NPROB,NVARF,NODEF)
!
! *** ���ό`��̍��W�ŉ׏d���v�Z��
!
      NNN=1
  CALL LOAD(CAL,COORD,F,FF,FINIT,NODEC,NVARC,NDOFN,NELEM,NELEC,NFINAL,&
            NNN,NODEM,NPOIN,NPRINT,PRES,SCF,SNOW,WCF,&
            WIND,NELEF,NODEF,CALFR,NVARF,MPOND,PROPC,PROPF,PROPM,NVARM,LOADCREM,WWIND,NTIMES)
!
! *** ���v�f�����}�Ɖ��͕��z�}�� VECTOR SCRIPT �`���̃t�@�C���ŏo�͂��遄
!
  IF(MINCR.NE.2.AND.LASTCYL.GE.2)NTIMES=LOADCREM    ! �j���[�g���@�̔����v�Z�Ŏ������Ȃ������ꍇ�iLASTCYL>1)
! IF(LOADCREM.EQ.1.OR.LOADCREM.EQ.NTIMES)THEN
    CALL WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                       NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,ABS(NPROB),NTIMES,&
                       STRSF,STRSM,STRSC,TTDIS)
    CALL CAL_PRES(NTIMES, NPOIN, NELEM, COORD, NODEM, PPRES, V_INIT, ATMOS_P, Z_REF)
! ENDIF
!
! *** ���Ȃ����ނ̌��聄
!
  CALL ALLOWANCE(COORD,NPOIN,NVARF,NDOFN,NELEF,NODEF,PROPF,STRSF,VALK)
!
  ! �j���[�g���@�̔����v�Z�Ŏ������Ȃ������ꍇ�iLASTCYL>1)�̏ꍇ�͌v�Z���I������B
  IF(MINCR.NE.2.AND.LASTCYL.GE.2)THEN
    WRITE(8,*)' LASTCYL > 1'
    STOP ' LASTCYL > 1'
  ENDIF
!
998 CONTINUE ! *** �����܂Łi�׏d����NTIMES�̃��[�v�j

close(8)
close(9)
close(10)
close(11)
close(12)

999 STOP
END
!
!
!
! =====================================================================
      SUBROUTINE ANGLE(E,ANGNE,ET1,ET2,POA1,POA2,GXY)
! =====================================================================
! �ޗ��̎厲�̕����ɍ��킹�Ēe���萔�}�g���N�X��ϊ�����
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION E(3,3)
      SITA=ANGNE*ATAN(1.0)/45.0
      POA12=1.0-POA1*POA2
      A11=ET1/POA12
      A12=POA2*A11
      A22=ET2/POA12
      A33=GXY
      C=COS(SITA)
      S=SIN(SITA)
      E(1,1)=C**4*A11++S**4*A22+2.0*C**2*S**2*(A12+2.0*A33)
      E(1,2)=C**2*S**2*(A11+A22-4.0*A33)+(C**4+S**4)*A12
      E(1,3)=-C*S*(C**2*A11-S**2*A22)+C*S*(C**2-S**2)*(A12+2.0*A33)
      E(2,1)=E(1,2)
      E(2,2)=S**4*A11+C**4*A22+2.0*C**2*S**2*(A12+2.0*A33)
      E(2,3)=-C*S*(S**2*A11-C**2*A22)-C*S*(C**2-S**2)*(A12+2.0*A33)
      E(3,1)=E(1,3)
      E(3,2)=E(2,3)
      E(3,3)=C**2*S**2*(A11-2.0*A12+A22)+(C**2-S**2)**2*A33
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE CABLE2(BL,CAL,COORD,COOR1,DV,IFFIX,INCREM,             &
                        NODEC,NVARC,MSIZE,MTOTV,              &
                        NCYCLE,NDOFN,NEA,NELEC,NFPOIN,                   &
                        NNN,NORDR,NPOIN,NPRET,NPRSTR,NTOTV,NUNBAL,      &
                        NWSUM,Q,TSTIF,STRC1,PROPC,            &
                        MCTEMP,CTEMP2)
! =====================================================================
! �P�[�u���v�f�̍����}�g���N�X���쐬����
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NORDR(MTOTV)
      DIMENSION NWSUM(0:MTOTV)
      DIMENSION TSTIF(MSIZE)
      DIMENSION IFFIX(MTOTV)
      DIMENSION COORD(NPOIN,3)
      DIMENSION COOR1(NPOIN,3)
      DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION BL(NELEC)
      DIMENSION CAL(NELEC)
      DIMENSION STRC1(NELEC)
      DIMENSION NEA(NELEC)
      dimension PROPC(NVARC,3)    !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj
      DIMENSION Q(MTOTV)
      DIMENSION PP(MTOTV)
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
!     MCTEMP                        ! MCTEMP�F�P�[�u���̉��x���͉�͂��s��Ȃ��i=0�j�A�s���i=1�j
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)�F���x���͂�^����(=1.0)�A�^���Ȃ�(=0.0)�ACTEMP2(2,IELEC)�F���x�Ђ���
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DIMENSION SM(6,6),IND(6),TT(6,6),U(6),DISE(6),BMT(6),GM(6,6)
      DIMENSION BB(6,6),BBG(6,6),BEBT(6,6),P(6)
!
      CALL ZEROR1(PP,MTOTV,NTOTV)  ! �P�[�u���܂��̓g���X�ɂ�铙���ߓ_�̓x�N�g��
!
      DO 210 NE=1,NELEC
        IV=NODEC(NE,3)
        EA=PROPC(IV,1)*ABS(PROPC(IV,2))  !EA[N]=�����O��[N/mm2]*�f�ʐ�[mm2]
        CALL ZEROR2(TT,6,6,6,6)
        I=NODEC(NE,1)
        J=NODEC(NE,2)
        CL=SQRT((COOR1(J,1)-COOR1(I,1))**2+(COOR1(J,2)-COOR1(I,2))**2+(COOR1(J,3)-COOR1(I,3))**2)
        TT(1,1)=(COOR1(J,1)-COOR1(I,1))/CL
        TT(1,2)=(COOR1(J,2)-COOR1(I,2))/CL
        TT(1,3)=(COOR1(J,3)-COOR1(I,3))/CL
        TT(3,3)=DSQRT(TT(1,1)*TT(1,1)+TT(1,2)*TT(1,2))
        IF(TT(3,3).EQ.0.0) TT(3,3)=1.0
        TT(2,1)=-TT(1,2)/TT(3,3)
        TT(2,2)= TT(1,1)/TT(3,3)
        TT(3,1)=-TT(1,1)*TT(1,3)/TT(3,3)
        TT(3,2)=-TT(1,2)*TT(1,3)/TT(3,3)
        DO II=1,3
          DO JJ=1,3
            TT(II+3,JJ+3)=TT(II,JJ)
          ENDDO
        ENDDO
!
        DO IDIR=1,3
          U(IDIR  )=COORD(I,IDIR)-COOR1(I,IDIR)
          U(IDIR+3)=COORD(J,IDIR)-COOR1(J,IDIR)
        ENDDO
!
        CALL ZEROR1(DISE,6,6)
        DO II=1,6
          DO KK=1,6
            DISE(II)=DISE(II)+TT(II,KK)*U(KK)
          ENDDO
        ENDDO
!
        BMT(1)=-1.0/CL
        BMT(2)=-(DISE(5)-DISE(2))/(CL*CL)
        BMT(3)=-(DISE(6)-DISE(3))/(CL*CL)
        BMT(4)=-BMT(1)
        BMT(5)=-BMT(2)
        BMT(6)=-BMT(3)
!
        BL(NE)=DSQRT((COORD(J,1)-COORD(I,1))**2+(COORD(J,2)-COORD(I,2))**2+(COORD(J,3)-COORD(I,3))**2)
        IF(NPRET.EQ.0.AND.NPRSTR.EQ.0.AND.NCYCLE.EQ.1.AND.INCREM.EQ.1.AND.NNN.EQ.1)THEN
          BL(NE)=1.00001*CAL(NE)
        ENDIF
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@ GRANPA DOME �̃����O���i�g�����j�A���f���j�̌v�Z @@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        EEA=EA
        IF(PROPC(IV,2).GT.0.0.AND.NEA(NE).EQ.0)THEN     !PROPC(IV,2)�f�ʐ�
          EEA=EEA/DV
        ENDIF
!
        IF(MCTEMP.EQ.0)THEN     ! MCTEMP=0�F�P�[�u���̉��x���͉�͂��s��Ȃ��ꍇ
          SSION=(BL(NE)-CAL(NE))*EA/CAL(NE)
          SION =BL(NE)*STRC1(NE)/CL+(BL(NE)-CL)*EEA/CL
        ELSE                    ! MCTEMP=1�F�P�[�u���̉��x���͉�͂��s���ꍇ                    
          CAL1=CAL(NE)*(1.0+CTEMP2(2,NE))
          SSION=((BL(NE)-CAL1)*EA/CAL1)
          SION =((BL(NE)-CAL1)*EA/CAL1)
        ENDIF
!
        IF(NUNBAL.NE.0)THEN
          NEA(NE)=1
          IF(SSION.LT.0.0.AND.PROPC((IV),2).GT.0.0)THEN     !PROPC(IV,2)�f�ʐ�
            SION=0.0
            NEA(NE)=0
          ENDIF
        ENDIF
!
! *** ��NUNBAL=0�F�j���[�g���@�ɂ������v�Z����
!
        IF(NUNBAL.EQ.0)THEN
!
          CALL ZEROR2(GM,6,6,6,6)
          TL=SION/CL
          GM(2,2)=TL
          GM(3,3)=TL
          GM(5,5)=TL
          GM(6,6)=TL
          GM(2,5)=-TL
          GM(3,6)=-TL
          GM(5,2)=-TL
          GM(6,3)=-TL
          CALEA=CL*EEA
!
          DO II=1,6
            DO JJ=1,6
              BB(II,JJ)=BMT(II)*BMT(JJ)*CALEA
            ENDDO
          ENDDO
!
          DO II=1,6
            DO JJ=1,6
              BBG(II,JJ)=BB(II,JJ)+GM(II,JJ)
            ENDDO
          ENDDO
!
          CALL ZEROR2(BEBT,6,6,6,6)
          DO II=1,6
            DO JJ=1,6
              DO KK=1,6
                BEBT(II,JJ)=BEBT(II,JJ)+BBG(II,KK)*TT(KK,JJ)
              ENDDO
            ENDDO
          ENDDO
!
          CALL ZEROR2(SM,6,6,6,6)
          DO II=1,6
            DO JJ=1,6
              DO KK=1,6
                SM(II,JJ)=SM(II,JJ)+TT(KK,II)*BEBT(KK,JJ)
              ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
        CALL ZEROR1(P,6,6)
        DO II=1,6
          DO JJ=1,6
            P(II)=P(II)+TT(JJ,II)*BMT(JJ)*CL*SION
          ENDDO
        ENDDO
!
! *** ���S�̍����}�g���N�X�ւ̏d�ˍ��킹��
!
        IF(NUNBAL.EQ.0)THEN
          CALL MATRIX(SM,NE,NODEC,NELEC,MSIZE,NDOFN,3,2,NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
        ENDIF
!
! *** ���P�[�u�����͂̓����ߓ_�̓x�N�g���o�o�̌v�Z�� 
!
        DO IDIR=1,3
          ITOTV=(I-1)*NDOFN+IDIR
          JTOTV=(J-1)*NDOFN+IDIR
          PP(ITOTV)=PP(ITOTV)+P(IDIR)
          PP(JTOTV)=PP(JTOTV)+P(IDIR+3)
        ENDDO
!
  210 CONTINUE
!
! *** ���P�[�u�����͂̓����ߓ_�̓x�N�g���o�o��S�̂̐ߓ_�̓x�N�g���p�ɉ����遄 
!
      DO ITOTV=1,NTOTV
        Q(ITOTV)=Q(ITOTV)+PP(ITOTV)
      ENDDO
!
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE CONVER(ERR,IFFIX,NPOIN,NCHECK,NDOFN,NFPOIN,Q)
! =====================================================================
! *** �c���͂̎��������lERR�ɑ΂��锻��i�������Ă����NCHECK=1�j
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IFFIX(NPOIN*NDOFN)
      DIMENSION Q(NPOIN*NDOFN)
!
      DO IPOIN=1,NFPOIN
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          IF(IFFIX(ITOTV).EQ.0.AND.ABS(Q(ITOTV)).GT.ERR)THEN
            NCHECK=0
            RETURN
          ENDIF
        ENDDO
      ENDDO
      NCHECK=1
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE CABLE3(INDEX,BL,CAL,NELEC,NVARC,NODEC,STRSC,STRC1,PROPC,MCTEMP,CTEMP2)
! =====================================================================
! �P�[�u���v�f�̊ɂ݂ƒ��͂��v�Z���A���ʂ��t�@�C���ɏo�͂���
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BL(NELEC)
      DIMENSION CAL(NELEC)
      DIMENSION STRSC(NELEC)
      DIMENSION STRC1(NELEC)
      DIMENSION NODEC(NELEC,3)       !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      dimension PROPC(NVARC,3)    !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj
      DIMENSION VALMM(NVARC,2)      ! VALMM(ICVR, )�F����ICVR�̍ő剞�́A�ŏ�����
      DIMENSION NVALMM(NVARC,2)     ! NVALMM(ICVR, )�F����ICVR�̍ő剞�́A�ŏ����͂ɑΉ�����v�f�ԍ�
!      DIMENSION NFAIL(5000)         ! NFAIL�F���݂𐶂��Ă���P�[�u���̗v�f�ԍ�
integer, allocatable :: NFAIL(:)    ! ���݂𐶂��Ă���P�[�u���̗v�f�ԍ�
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
!     MCTEMP                        ! MCTEMP�F�P�[�u���̉��x���͉�͂��s��Ȃ��i=0�j�A�s���i=1�j
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)�F���x���͂�^����(=1.0)�A�^���Ȃ�(=0.0)�ACTEMP2(2,IELEC)�F���x�Ђ���
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IFAIL=0   ! IFAIL�F���݂������Ă���P�[�u���̖{��
!
! *** ���P�[�u�����ށi�f�ʐς̕������v���X�̕��ށj�Ɉ��k�͂������Ă���ꍇ�ɂ́A�[���ɂ��遄
      DO NE=1,NELEC
        IV=NODEC(NE,3)
        EA=PROPC(IV,1)*ABS(PROPC(IV,2))   !EA[N]=�����O��[N/mm2]*�f�ʐ�[mm2]
        IF(MCTEMP.EQ.0)THEN     ! MCTEMP=0�F�P�[�u���̉��x���͉�͂��s��Ȃ��ꍇ
          STRC1(NE)=(BL(NE)-CAL(NE))*EA/CAL(NE)
        ELSE                    ! MCTEMP=1�F�P�[�u���̉��x���͉�͂��s���ꍇ
          STRC1(NE)=((BL(NE)-CAL(NE))*EA/CAL(NE))-CTEMP2(2,NE)*EA
        ENDIF
        STRSC(NE)=STRC1(NE)
        IF((STRSC(NE).LT.0.0).AND.(PROPC(IV,2).GT.0.0))THEN        !PROPC(IV,2)�f�ʐ�
          STRSC(NE)=0.0
          IFAIL=IFAIL+1
        ENDIF
      ENDDO
!
allocate (NFAIL(IFAIL))           ! ���݂𐶂��Ă���P�[�u���̗v�f�ԍ�
i=0
do NE=1,NELEC
  if(STRSC(NE).EQ.0.0)then
    i=i+1
    NFAIL(i)=NE
  ENDIF
enddo
!
      IF(INDEX.EQ.1)THEN
        WRITE(8,120)
        WRITE(8,130)(IELEC,STRSC(IELEC),IELEC=1,NELEC)
      ENDIF
  120 FORMAT(//,' (FORCES OF CABLE MEMBERS) [N]',/,5(' ELEM',5X,'FORCE',10X))
  130 FORMAT(5(I5,E15.6,5X))
!
! *** ���P�[�u�����̉��͂̍ő�lVALMM(ISTRE,1)�ƍŏ��lVALMM(ISTRE,2)���t�@�C�����W�ɕۑ�����iINDEX=1�̏ꍇ�j��
!
      IF(INDEX.EQ.1)THEN
        CALL ZEROI2(NVALMM,NVARC,2,NVARC,2)
        DO IELEC=1,NELEC
          ICVR=NODEC(IELEC,3)
          IF(NVALMM(ICVR,1).EQ.0)THEN
            DO NUM=1,2
              VALMM(ICVR,NUM)=STRSC(IELEC)
              NVALMM(ICVR,NUM)=IELEC
            ENDDO
          ELSE
            IF(STRSC(IELEC).GT.VALMM(ICVR,1))THEN
              VALMM(ICVR,1)=STRSC(IELEC)
              NVALMM(ICVR,1)=IELEC
            ELSEIF(STRSC(IELEC).LT.VALMM(ICVR,2))THEN
              VALMM(ICVR,2)=STRSC(IELEC)
              NVALMM(ICVR,2)=IELEC
            ENDIF
          ENDIF
        ENDDO
        WRITE(8,816)(j,(NVALMM(j,i),VALMM(j,i),i=1,2),j=1,NVARC)
      ENDIF
  816 FORMAT(/,' (MAXIMUM AND MINIMUM OF CABLE STRESS) [N]',/,'     NVARC      ELEM     MAX-STRESS      ELEM     MIN-STRESS',/,(I10,I10,E15.5,I10,E15.5))
!
! *** �����݂��������P�[�u���̌�IFAIL�Ƃ��̕��ޔԍ�NFAIL(IFAIL)���t�@�C���ɏ������ށ�
!
      IF(IFAIL.GE.1)WRITE(8,150)IFAIL,(NFAIL(JFAIL),JFAIL=1,IFAIL)
  150 FORMAT(' **** NUMBER OF SLACKED CABLE ELEMENTS (',I5,' ) ****',/,(10X,20I5))
!
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE MEMBRN(ANG,AREA,COOR1,COORD,DDV,EE,IFFIX,  &
                        INCREM,MEME1,MEME2,NPOIN,MSIZE,MSTRE,     &
                        MTOTV,NANG,NANGLE,NCYCLE,NDOFN,NELEM,NFINAL,    &
                        NNN,NODEM,NONISO,NORDR,NPROB,NUNBAL,NVARM,NWSUM,   &
                        PRETEN,PROPM,Q,STRM1,STRSM,TSTIF)
! =====================================================================
! **** �����v�f�̍����}�g���N�X���쐬���遄 ***
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION ANG(NELEM)
      DIMENSION AREA(NELEM)
      DIMENSION COOR1(NPOIN,3)
      DIMENSION COORD(NPOIN,3)
      DIMENSION EE(3,3)
      DIMENSION IFFIX(NPOIN*NDOFN)
      DIMENSION MEME1(NELEM)
      DIMENSION MEME2(NELEM)
      DIMENSION NANG(NELEM)
      DIMENSION NODEM(NELEM,4)          ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NORDR(NPOIN*NDOFN)
      DIMENSION NWSUM(0:NPOIN*NDOFN)
      DIMENSION PROPM(NVARM,7,3)        !���v�f�̍ޗ������i��ޔԍ��A��������(1�`7)�A�i�KTri-Linear(1�`3)�j
      DIMENSION Q(NPOIN*NDOFN)
      DIMENSION STRM1(NELEM,MSTRE)
      DIMENSION STRSM(NELEM,MSTRE)
      DIMENSION TSTIF(MSIZE)
      DIMENSION E(3,3)
      DIMENSION U(9)
      DIMENSION SM(9,9)
      DIMENSION P(9)
!
      DO 520 NE=1,NELEM
        II=NODEM(NE,1)
        JJ=NODEM(NE,2)
        KK=NODEM(NE,3)
IVARM=NODEM(NE,4)   !NODEM(NE,4)�F���ޔԍ�
        XXI=COOR1(II,1)                                                 ! XX(II)
        XXJ=COOR1(JJ,1)                                                 ! XX(JJ)
        XXK=COOR1(KK,1)                                                 ! XX(KK)
        YYI=COOR1(II,2)                                                 ! YY(II)
        YYJ=COOR1(JJ,2)                                                 ! YY(JJ)
        YYK=COOR1(KK,2)                                                 ! YY(KK)
        ZZI=COOR1(II,3)                                                 ! ZZ(II)
        ZZJ=COOR1(JJ,3)                                                 ! ZZ(JJ)
        ZZK=COOR1(KK,3)                                                 ! ZZ(KK)
        PX=STRSM(NE,1)+PRETEN
        PY=STRSM(NE,2)+PRETEN
        PXY=STRSM(NE,3)
        PPX=STRM1(NE,1)+PRETEN
        PPY=STRM1(NE,2)+PRETEN
        PPXY=STRM1(NE,3)

!
! *** ��Tri-Linear���f���ɂ��ޗ������̐ݒ聄 ***
!
  CNV=1.0E+00               ! CNV�F�P�ʂ̊��Z�W���i������mm�P�ʂ̏ꍇ1.0�A������m�P�ʂ̏ꍇ1000.0�j
  TH=PROPM(IVARM,6,1)       ! PROPM(i,6,1)�F����[mm]
  F=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)    ! F�F��������[N/m]
  F=F/CNV/TH             ! F�F�P�ʂ̕ϊ��i���W��mm�P�ʂ̏ꍇN/mm2�A���W��m�P�ʂ̏ꍇN/m��N/mm2�j
  if(F.lt.PROPM(IVARM,5,1))then              ! PROPM(IVRAM,5,1)�F��P�~���_���̓�y
    j=1
  elseif(F.lt.PROPM(IVARM,5,2))then          ! PROPM(IVRAM,5,2)�F��Q�~���_���̓�y
    j=2
  else
    j=3                                 ! j=1,2,3�i��P�`��R���z�j
  endif
  ET1=PROPM(IVARM,1,j)*TH*CNV     !��������ET=�����O��E*����t
  ET2=PROPM(IVARM,2,j)*TH*CNV
  POA1=PROPM(IVARM,3,j)
  POA2=POA1*ET2/ET1           ! �����̒藝�𖞑�����悤��POA2��ݒ肷��
  GXY=PROPM(IVARM,4,j)*TH*CNV

!  ET1=PROPM(IVARM,1,1)
!  ET2=PROPM(IVARM,2,1)
!  POA1=PROPM(IVARM,3,1)
!  POA2=POA1*ET2/ET1           ! �����̒藝�𖞑�����悤��POA2��ݒ肷��
!  GXY=PROPM(IVARM,4,1)

        ETT1=ET1
        ETT2=ET2
        EXY=GXY
        POA11=POA1
        POA22=POA2
        IF(INCREM.GT.1.OR.NNN.GT.1)THEN
          IF(MEME1(NE).EQ.0)THEN
            ETT1=ETT1/DDV
            POA11=POA11/DDV
            EXY=EXY/DDV
          ENDIF
          IF(MEME2(NE).EQ.0)THEN
            ETT2=ETT2/DDV
            POA22=POA22/DDV
            EXY=EXY/DDV
          ENDIF
        ENDIF
!
! ���� EE��MAIN PROGRAM���Ɉړ������v���O���������邪�A�ޗ�����`�ŉ�͂��s���ꍇ�ɂ́A�����ōX�V����K�v������B
        EE(1,1)=ET1/(1.0-POA1*POA2)
        EE(1,2)=EE(1,1)*POA2
        EE(1,3)=0.0
        EE(2,1)=EE(1,2)
        EE(2,2)=ET2/(1.0-POA1*POA2)
        EE(2,3)=0.0
        EE(3,1)=0.0
        EE(3,2)=0.0
        EE(3,3)=GXY
!
        E(1,1)=ETT1/(1.0-POA11*POA22)
        E(1,2)=E(1,1)*POA22
        E(1,3)=0.0
        E(2,1)=ETT2*POA11/(1.0-POA11*POA22)
        E(2,2)=ETT2/(1.0-POA11*POA22)
        E(2,3)=0.0
        E(3,1)=0.0
        E(3,2)=0.0
        E(3,3)=EXY
        U(1)=COORD(II,1)-COOR1(II,1)                                    ! X(II)-XX(II)
        U(2)=COORD(II,2)-COOR1(II,2)                                    ! Y(II)-YY(II)
        U(3)=COORD(II,3)-COOR1(II,3)                                    ! Z(II)-ZZ(II)
        U(4)=COORD(JJ,1)-COOR1(JJ,1)                                    ! X(JJ)-XX(JJ)
        U(5)=COORD(JJ,2)-COOR1(JJ,2)                                    ! Y(JJ)-YY(JJ)
        U(6)=COORD(JJ,3)-COOR1(JJ,3)                                    ! Z(JJ)-ZZ(JJ)
        U(7)=COORD(KK,1)-COOR1(KK,1)                                    ! X(KK)-XX(KK)
        U(8)=COORD(KK,2)-COOR1(KK,2)                                    ! Y(KK)-YY(KK)
        U(9)=COORD(KK,3)-COOR1(KK,3)                                    ! Z(KK)-ZZ(KK)
        SURF=AREA(NE)
!
! *** �����v�fIJK��IJ�����ƍޗ��̎厲�������X���Ă���ꍇ�̒e�������}�g���N�X�c�̕ϊ���
!
        IF(NONISO.EQ.2)THEN
          DO I=1,NANGLE
            IF(NE.EQ.NANG(I))THEN
                   IF(I.EQ.1)CALL WMTR2(2,E,3,3,3,3,'E_ORG','NE   ',NE)
              CALL ANGLE(E,ANG(I),ETT1,ETT2,POA11,POA22,EXY)
              CALL ANGLE(EE,ANG(I),ETT1,ETT2,POA11,POA22,EXY)
                   IF(I.EQ.1)CALL WMTR2(2,E,3,3,3,3,'E_CNV','NE   ',NE)
              EXIT
            ENDIF
          ENDDO
        ENDIF
! *** ���v�f�����}�g���N�X�̍쐬�� ***
!
        CALL NEWTON(SM,PX,PY,PXY,P,XXI,YYI,ZZI,XXJ,YYJ,ZZJ,XXK,YYK,ZZK,U,NUNBAL,E,PPX,PPY,PPXY,EE,NFINAL,NE,SURF)
!
! *** �������N�����O�̔���Ə����� ***
!
        IF(NUNBAL.NE.0.AND.(NCYCLE.GT.1.OR.INCREM.GT.1.OR.NNN.GT.1))THEN
          STRM1(NE,1)=PPX
          STRM1(NE,2)=PPY
          STRM1(NE,3)=PPXY
          STRSM(NE,1)=PX
          STRSM(NE,2)=PY
          STRSM(NE,3)=PXY
          MEME2(NE)=1
          MEME1(NE)=1
          SQ=SQRT((PPX-PPY)*(PPX-PPY)*0.25+PPXY*PPXY)
          PP1=(PPX+PPY)*0.5+SQ
          PP2=(PPX+PPY)*0.5-SQ
          IF(PP1.LE.0.0)THEN
            MEME1(NE)=0
            MEME2(NE)=0
          ELSE
            IF(PP2.LE.0.0)THEN
              IF(PPX.LE.PPY.AND.PPX.LT.0.0) MEME1(NE)=0
              IF(PPX.GT.PPY.AND.PPY.LT.0.0) MEME2(NE)=0
            ENDIF
          ENDIF
        ENDIF
!
! *** �������͂̓����ߓ_�͂�S�̂̐ߓ_�̓x�N�g���p�ɉ����遄 ***
!
        DO IDOFN=1,3
          ITOTV=(II-1)*NDOFN+IDOFN
          Q(ITOTV)=Q(ITOTV)+P(IDOFN)
          JTOTV=(JJ-1)*NDOFN+IDOFN
          Q(JTOTV)=Q(JTOTV)+P(IDOFN+3)
          KTOTV=(KK-1)*NDOFN+IDOFN
          Q(KTOTV)=Q(KTOTV)+P(IDOFN+6)
        ENDDO
!
! *** ���v�f�����}�g���b�N�X���S�̍����}�g���b�N�X�� ***
!
        CALL MATRIX(SM,NE,NODEM,NELEM,MSIZE,NDOFN,3,3,NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
  520 CONTINUE
!
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE NEWTON(SM,PX,PY,PXY,P,XXI,YYI,ZZI,XXJ,YYJ,ZZJ,XXK,YYK, &
                        ZZK,U,NUNBAL,E,PPX,PPY,PPXY,EE,NFINAL,NE,SURF)
! =====================================================================
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TT(9,9),B0MT(3,9),BLMT(3,9),SE(9),DISE(9),GM(9,9),      &
                E(3,3),BEB(9,9),BEBT(9,9),SM(9,9),U(9),P(9),EB0(3,9),   &
                EBL(3,9),B0EB0(9,9),B0EBL(9,9),BLEB0(9,9),BLEBL(9,9),   &
                TEN(3),EE(3,3),EEB0(3,9),EEBL(3,9),TTEN(3)
!
      XXJK=XXJ-XXK
      YYJK=YYJ-YYK
      ZZJK=ZZJ-ZZK
      XXKI=XXK-XXI
      YYKI=YYK-YYI
      ZZKI=ZZK-ZZI
      XXIJ=XXI-XXJ
      YYIJ=YYI-YYJ
      ZZIJ=ZZI-ZZJ
      AA=-YYIJ*ZZKI+ZZIJ*YYKI
      BB=-ZZIJ*XXKI+XXIJ*ZZKI
      CC=-XXIJ*YYKI+YYIJ*XXKI
      S2=DSQRT(AA*AA+BB*BB+CC*CC)
      AALIJ=DSQRT(XXIJ*XXIJ+YYIJ*YYIJ+ZZIJ*ZZIJ)
      SS2=SURF
!
      CALL ZEROR2(TT,9,9,9,9)
      TT(1,1)=-XXIJ/AALIJ
      TT(1,2)=-YYIJ/AALIJ
      TT(1,3)=-ZZIJ/AALIJ
      AALIJ2=AALIJ*S2
      TT(2,1)=(CC*YYIJ-BB*ZZIJ)/AALIJ2
      TT(2,2)=(AA*ZZIJ-CC*XXIJ)/AALIJ2
      TT(2,3)=(BB*XXIJ-AA*YYIJ)/AALIJ2
      TT(3,1)=AA/S2
      TT(3,2)=BB/S2
      TT(3,3)=CC/S2
      DO I=1,3
        DO J=1,3
          TT(I+3,J+3)=TT(I,J)
          TT(I+6,J+6)=TT(I,J)
        ENDDO
      ENDDO
!
      YYEJK=TT(2,1)*XXJK+TT(2,2)*YYJK+TT(2,3)*ZZJK
      XXEJK=TT(1,1)*XXJK+TT(1,2)*YYJK+TT(1,3)*ZZJK
      YYEKI=TT(2,1)*XXKI+TT(2,2)*YYKI+TT(2,3)*ZZKI
      XXEKI=TT(1,1)*XXKI+TT(1,2)*YYKI+TT(1,3)*ZZKI
      YYEIJ=TT(2,1)*XXIJ+TT(2,2)*YYIJ+TT(2,3)*ZZIJ
      XXEIJ=TT(1,1)*XXIJ+TT(1,2)*YYIJ+TT(1,3)*ZZIJ
!
      CALL ZEROR1(DISE,9,9)
      DO I=1,9
        DO K=1,9
          DISE(I)=DISE(I)+TT(I,K)*U(K)
        ENDDO
      ENDDO
!
      ALFA2=( YYEJK*DISE(1)+YYEKI*DISE(4)+YYEIJ*DISE(7))/SS2
      ALFA3=(-XXEJK*DISE(1)-XXEKI*DISE(4)-XXEIJ*DISE(7))/SS2
      ALFA5=( YYEJK*DISE(2)+YYEKI*DISE(5)+YYEIJ*DISE(8))/SS2
      ALFA6=(-XXEJK*DISE(2)-XXEKI*DISE(5)-XXEIJ*DISE(8))/SS2
      ALFA8=( YYEJK*DISE(3)+YYEKI*DISE(6)+YYEIJ*DISE(9))/SS2
      ALFA9=(-XXEJK*DISE(3)-XXEKI*DISE(6)-XXEIJ*DISE(9))/SS2
!
!     <B0>,<BL> MATRIX
      CALL ZEROR2(B0MT,3,9,3,9)
      B0MT(1,1)= YYEJK/SS2
      B0MT(1,4)= YYEKI/SS2
      B0MT(1,7)= YYEIJ/SS2
      B0MT(2,2)=-XXEJK/SS2
      B0MT(2,5)=-XXEKI/SS2
      B0MT(2,8)=-XXEIJ/SS2
      B0MT(3,1)=B0MT(2,2)
      B0MT(3,2)=B0MT(1,1)
      B0MT(3,4)=B0MT(2,5)
      B0MT(3,5)=B0MT(1,4)
      B0MT(3,7)=B0MT(2,8)
      B0MT(3,8)=B0MT(1,7)
      SS4=SS2*2.0
      BLMT(1,1)= YYEJK*ALFA2/SS2
      BLMT(1,2)= YYEJK*ALFA5/SS2
      BLMT(1,3)= YYEJK*ALFA8/SS2
      BLMT(1,4)= YYEKI*ALFA2/SS2
      BLMT(1,5)= YYEKI*ALFA5/SS2
      BLMT(1,6)= YYEKI*ALFA8/SS2
      BLMT(1,7)= YYEIJ*ALFA2/SS2
      BLMT(1,8)= YYEIJ*ALFA5/SS2
      BLMT(1,9)= YYEIJ*ALFA8/SS2
      BLMT(2,1)=-XXEJK*ALFA3/SS2
      BLMT(2,2)=-XXEJK*ALFA6/SS2
      BLMT(2,3)=-XXEJK*ALFA9/SS2
      BLMT(2,4)=-XXEKI*ALFA3/SS2
      BLMT(2,5)=-XXEKI*ALFA6/SS2
      BLMT(2,6)=-XXEKI*ALFA9/SS2
      BLMT(2,7)=-XXEIJ*ALFA3/SS2
      BLMT(2,8)=-XXEIJ*ALFA6/SS2
      BLMT(2,9)=-XXEIJ*ALFA9/SS2
      BLMT(3,1)=(YYEJK*ALFA3-XXEJK*ALFA2)/SS2
      BLMT(3,2)=(YYEJK*ALFA6-XXEJK*ALFA5)/SS2
      BLMT(3,3)=(YYEJK*ALFA9-XXEJK*ALFA8)/SS2
      BLMT(3,4)=(YYEKI*ALFA3-XXEKI*ALFA2)/SS2
      BLMT(3,5)=(YYEKI*ALFA6-XXEKI*ALFA5)/SS2
      BLMT(3,6)=(YYEKI*ALFA9-XXEKI*ALFA8)/SS2
      BLMT(3,7)=(YYEIJ*ALFA3-XXEIJ*ALFA2)/SS2
      BLMT(3,8)=(YYEIJ*ALFA6-XXEIJ*ALFA5)/SS2
      BLMT(3,9)=(YYEIJ*ALFA9-XXEIJ*ALFA8)/SS2
      IF(NFINAL.EQ.1) GO TO 160
!
      CALL ZEROR2(EB0,3,9,3,9)
      CALL ZEROR2(EBL,3,9,3,9)
      CALL ZEROR2(EEB0,3,9,3,9)
      CALL ZEROR2(EEBL,3,9,3,9)
      DO I=1,3
        DO J=1,9
          DO K=1,3
            EB0(I,J)=EB0(I,J)+E(I,K)*B0MT(K,J)
            EBL(I,J)=EBL(I,J)+E(I,K)*BLMT(K,J)
            EEB0(I,J)=EEB0(I,J)+EE(I,K)*B0MT(K,J)
            EEBL(I,J)=EEBL(I,J)+EE(I,K)*BLMT(K,J)
          ENDDO
        ENDDO
      ENDDO
!
      EPRX=ALFA2+0.5*(ALFA2*ALFA2+ALFA5*ALFA5+ALFA8*ALFA8)
      EPRY=ALFA6+0.5*(ALFA3*ALFA3+ALFA6*ALFA6+ALFA9*ALFA9)
      GAMMA=ALFA3+ALFA5+ALFA2*ALFA3+ALFA5*ALFA6+ALFA8*ALFA9
      DO I=1,3
        TEN(I)=E(I,1)*EPRX+E(I,2)*EPRY+E(I,3)*GAMMA
        TTEN(I)=EE(I,1)*EPRX+EE(I,2)*EPRY+EE(I,3)*GAMMA
      ENDDO
!
      PPX=TTEN(1)+PPX
      PPY=TTEN(2)+PPY
      PPXY=TTEN(3)+PPXY
      PX=TEN(1)+PX
      PY=TEN(2)+PY
      PXY=TEN(3)+PXY
      IF(NUNBAL.NE.1) GO TO 160
      PPX=TTEN(1)
      PPY=TTEN(2)
      PPXY=TTEN(3)
      PX=TEN(1)
      PY=TEN(2)
      PXY=TEN(3)
  160 CONTINUE
      IF(NUNBAL.EQ.0) GO TO 165
      SQ=DSQRT(0.25*(PPX-PPY)**2+PPXY**2)
      PP1=0.5*(PPX+PPY)+SQ
      PP2=0.5*(PPX+PPY)-SQ
      IF(PPX.EQ.PPY) TNSITA=0.1D+09
      IF(PPX.EQ.PPY) GO TO 161
      TNSITA=2.0*PPXY/(PPX-PPY)
  161 PP3=28.66*DATAN(TNSITA)
      IF(PP2.LT.0.0) GO TO 162
      PX=PPX
      PY=PPY
      PXY=PPXY
      GO TO 165
  162 IF(PP1.LT.0.0) GO TO 164
      IF(PPX.LT.PPY) GO TO 163
      PX=0.5*PP1*(1.0+1.0/DSQRT(1.0+TNSITA*TNSITA))
      PY=PP1-PX
      PXY=TNSITA*(PX-PY)*0.5
      GO TO 165
  163 PY=0.5*PP1*(1.0+1.0/DSQRT(1.0+TNSITA*TNSITA))
      PX=PP1-PY
      PXY=TNSITA*(PX-PY)*0.5
      GO TO 165
  164 PX=0.0
      PY=0.0
      PXY=0.0
  165 CONTINUE
!
!     SE(I) : NODAL FORCES EQUIVALENT TO MEMBRANE FORCES IN LOCAL.
      DO I=1,9
        SE(I)=0.5*SS2*((B0MT(1,I)+BLMT(1,I))*PX+(B0MT(2,I)+BLMT(2,I))*PY+(B0MT(3,I)+BLMT(3,I))*PXY)
      ENDDO
!
      IF(NUNBAL.NE.0) GO TO 280
!
      CALL ZEROR2(GM,9,9,9,9)
      W=PX
      QQ=PY
      R=PXY
      GM(1,1)=YYEJK*YYEJK*W+XXEJK*XXEJK*QQ-2.0*XXEJK*YYEJK*R
      GM(1,4)=YYEJK*YYEKI*W+XXEJK*XXEKI*QQ-(YYEJK*XXEKI+XXEJK*YYEKI)*R
      GM(1,7)=YYEJK*YYEIJ*W+XXEJK*XXEIJ*QQ-(YYEJK*XXEIJ+YYEIJ*XXEJK)*R
      GM(4,4)=YYEKI*YYEKI*W+XXEKI*XXEKI*QQ-2.0*XXEKI*YYEKI*R
      GM(4,7)=YYEKI*YYEIJ*W+XXEKI*XXEIJ*QQ-(XXEKI*YYEIJ+YYEKI*XXEIJ)*R
      GM(7,7)=YYEIJ*YYEIJ*W+XXEIJ*XXEIJ*QQ-2.0*XXEIJ*YYEIJ*R
!
      DO I=1,2
        GM(I+1,I+1)=GM(I,I)
        GM(I+4,I+4)=GM(I+3,I+3)
        GM(I+7,I+7)=GM(I+6,I+6)
        GM(I+1,I+4)=GM(I,I+3)
        GM(I+1,I+7)=GM(I,I+6)
        GM(I+4,I+7)=GM(I+3,I+6)
      ENDDO
      DO I=1,6
        GM(I+3,I)=GM(I,I+3)
      ENDDO
      DO I=1,3
        GM(I+6,I)=GM(I,I+6)
      ENDDO
      DO I=1,9
        DO J=1,9
          GM(I,J)=GM(I,J)/SS4
        ENDDO
      ENDDO
!
      CALL ZEROR2(B0EB0,9,9,9,9)
      CALL ZEROR2(B0EBL,9,9,9,9)
      CALL ZEROR2(BLEB0,9,9,9,9)
      CALL ZEROR2(BLEBL,9,9,9,9)
      DO I=1,9
        DO J=1,9
          DO K=1,3
            B0EB0(I,J)=B0EB0(I,J)+B0MT(K,I)*EB0(K,J)
            B0EBL(I,J)=B0EBL(I,J)+B0MT(K,I)*EBL(K,J)
            BLEB0(I,J)=BLEB0(I,J)+BLMT(K,I)*EB0(K,J)
            BLEBL(I,J)=BLEBL(I,J)+BLMT(K,I)*EBL(K,J)
          ENDDO
        ENDDO
      ENDDO
      DO I=1,9
        DO J=1,9
          BEB(I,J)=B0EB0(I,J)+B0EBL(I,J)+BLEB0(I,J)+BLEBL(I,J)
          BEB(I,J)=BEB(I,J)*SS2*0.5
        ENDDO
      ENDDO
      DO I=1,9
        DO J=1,9
          BEB(I,J)=BEB(I,J)+GM(I,J)
        ENDDO
      ENDDO
!
      CALL ZEROR2(BEBT,9,9,9,9)
      DO I=1,9
        DO J=1,9
          DO K=1,9
            BEBT(I,J)=BEBT(I,J)+BEB(I,K)*TT(K,J)
          ENDDO
        ENDDO
      ENDDO
!
      CALL ZEROR2(SM,9,9,9,9)
      DO I=1,9
        DO J=1,9
          DO K=1,9
            SM(I,J)=SM(I,J)+TT(K,I)*BEBT(K,J)
          ENDDO
        ENDDO
      ENDDO
  280 CONTINUE
!
      CALL ZEROR1(P,9,9)
      DO I=1,9
        DO J=1,9
          P(I)=P(I)+TT(J,I)*SE(J)
        ENDDO
      ENDDO
!
      RETURN
      END
