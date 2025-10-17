C ==================================================================================================
C �v���O�����̍\��
C SKYLNE2M.FOR�F�X�J�C���C���@�ɂ��A���ꎟ�������̉�@�Ɋւ���v���O�����A�iCMEM7.FOR�d�l�j
C   FIXNO  NET�n�̍S���ߓ_�����ASKYLNE�^�ɕύX����
C   SIZING �S�̍����}�g���N�X�̐����ƂP�����z�񉻂��ꂽ�X�J�C���C���}�g���N�X�̑Ή��������w�W�z����쐬����
C   MATRIX �v�f�����}�g���N�X�̐�����S�̍����}�g���N�X�i�P�����z�񉻂��ꂽ�X�J�C���C���}�g���N�X�j�֑������
C   BCONUU �W���}�g���b�N�X����щE�Ӄx�N�g���ɑ΂��āA�􉽊w�I�ȏ�������������
C   RARNG  �e��ψʃx�N�g���̏��������ɖ߂�
C   SKYLNE �X�J�C���C���@�ɂ��A���ꎟ�������̉������߂�
C ======================================================================
      SUBROUTINE FIXNO(IFFIX,IFPRE,NODEM,NODEC,
     &                 MTOTV,NDOFN,NELEM,NELEC,NFPOIN,
     &                 NFREE,NOFIX,NORDR,NPOIN,NTOTV,NVFIX,
     &                 NELEF,NODEF)
C ======================================================================
C *** CMEM�^�̍S���ߓ_�����ASKYLNE�^�ɕύX����
C  * CMEM�^ 
C     NFPOIN            �F�����R�ߓ_��
C     NPOIN             �F���ߓ_��
C  * SKYLNE�^
C     LNODS(IELEM,INODE)�F�v�fIELEM�̗v�f���ߓ_INODE�̑S�̐ߓ_�ԍ�
C     NELEM             �F�v�f��
C     IFFIX(ITOTV)      �F�S���Ɋւ�����i0:���R�A1:�S��)
C     IFPRE(IVFIX)      �F�S����Ԃ�\���ϐ��i��@101:x,z�����S���Ay�������R�j
C     NDOFN             �F�P�ߓ_����̎��R�x��
C     NFREE             �F�����R�x���i�S��������)
C     NOFIX(IVFIX)      �F�S���ߓ_�̔ԍ�
C     NTOTV             �F�����R�x���i�S�����܂�), ��MTOTV
C     NVFIX             �F�S���ߓ_�̑���
C     
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IFFIX(MTOTV)
      DIMENSION NOFIX(NVFIX)
      DIMENSION NORDR(MTOTV)
	DIMENSION IFPRE(NVFIX)
      DIMENSION NODEM(NELEM,4)       ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NODEC(NELEC,3)     !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION NODEF(NELEF,3)       ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** FIXNO **********'
C
C *** �S�����������z��IFFIX�̐ݒ�A���ނɐڑ�����Ă��Ȃ��ߓ_�͍S���ߓ_IIFIX=1�Ƃ��ď�������
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
C *** �ߓ_�ԍ� NFPOIN+1 �` NPOIN �̎��R�x�����ׂčS������
C
      DO 400 IPOIN=NFPOIN+1,NPOIN
      DO 400 IDOFN=1,NDOFN
        ITOTV=(IPOIN-1)*NDOFN+IDOFN
        IFFIX(ITOTV)=1
  400 CONTINUE
C
C *** �����R�x��(�S��������)NFREE�̌v�Z
      NFREE=NTOTV
      DO 500 ITOTV=1,NTOTV
        NFREE=NFREE-IFFIX(ITOTV)
  500 CONTINUE
C
C *** NORDR(ITOTV):�S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ�����
C                  �������鎞�̐V�����S�̎��R�x�ԍ��̔z��̍쐬
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
C *** �S�̍����}�g���N�X�̐����ƂP�����z�񉻂��ꂽ�X�J�C���C���}�g���N�X�̑Ή��������w�W�z����쐬����B
C     NODEM(IELEM,INODE) :���v�fIELEM�̗v�f���ߓ_INODE�̑S�̐ߓ_�ԍ�
C     NODEC(IELEC,INODE):�P�[�u���v�fIELEC�̗v�f���ߓ_INODE�̑S�̐ߓ_�ԍ�
C     NTOTV       : �����R�x��(=NPOIN*NDOFN<MTOTV)
C     NELEM       : ���v�f��
C     NELEC       : �P�[�u���v�f��
C     NDOFN       : �P�ߓ_������̎��R�x��
C     NPOIN       : �ߓ_��
C     NORDR(ITOTV): �S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ����ĕ������鎞�̐V�����S�̎��R�x�ԍ�
C     NWDTH(ITOTV): �S�̍����}�g���N�X�̑�ITOTV�s�̑�P��[���v�f����Ίp�v�f�܂ł̌�
C     NWSUM(ITOTV): �S�̍����}�g���N�X�̑�P�`ITOTV�s�̑�P��[���v�f����Ίp�v�f�܂ł̗v�f���̗ݘa
C     LEFTS(ITOTV): �S�̍����}�g���N�X�̑�IPOIN�s�ڂ̑�P��[���v�f�̗�ԍ�
C
      PARAMETER(MCHK=0,MMTOTV=20000)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NODEM(NELEM,4)       ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NODEC(NELEC,3)       !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION NORDR(MTOTV)
      DIMENSION NWDTH(MTOTV)
      DIMENSION NWSUM(0:MTOTV)
      DIMENSION LEFTS(MMTOTV)
      DIMENSION NODEF(NELEF,3)       ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
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
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:���׊�����̍s�ԍ�
              DO 200 JDOFN=1,3                                          JTOTV:���׊�����̗�ԍ�
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
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:���׊�����̍s�ԍ�
              DO 250 JDOFN=1,3                                          JTOTV:���׊�����̗�ԍ�
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
              ITOTV=NORDR((IPOIN-1)*NDOFN+IDOFN)                        ITOTV:���׊�����̍s�ԍ�
              DO 270 JDOFN=1,6                                          JTOTV:���׊�����̗�ԍ�
                JTOTV=NORDR((JPOIN-1)*NDOFN+JDOFN)
                IF(LEFTS(ITOTV).GT.JTOTV)LEFTS(ITOTV)=JTOTV
  270 CONTINUE
C
      DO 300 IPOIN=1,NPOIN
      DO 300 IDOFN=1,NDOFN
        ITOTV=(IPOIN-1)*NDOFN+IDOFN                                     ITOTV:���׊����O�̍s�ԍ�
        JTOTV=NORDR(ITOTV)                                              JTOTV:���׊�����̍s�ԍ�
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
      MSIZE=NFSIZE    !���I�z��̃T�C�Y�ATSTIF(MSIZE)
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
C *** �v�f�����}�g���N�X�̐�����S�̍����}�g���N�X�i�P�����z�񉻂��ꂽ�X�J�C���C���}�g���N�X�j�֑������B
C     ESTIF(IEVAB,JEVAB) : �v�f�����}�g���N�X(IEVAB,JEVAB=1�`NDOF1*NNODE)
C     TSTIF(JPOSI)       : �S�̍����}�g���N�X�i�X�J�C���C���}�g���N�X�j
C     LNODS(IELEM,INODE) : �v�fIELEM�̗v�f���ߓ_INODE�̑S�̐ߓ_�ԍ�
C     IELEM   : �v�f�ԍ�
C     NDOFN   : �P�ߓ_������̎��R�x�� PARAMETER(NDOFN=6)
C     NDOF1   : �e�v�f�����}�g���N�X�ɂ����Ē�`����Ă���P�ߓ_������̎��R�x��
C     NNODE   : �e�v�f�����}�g���N�X�ɂ����Ē�`����Ă���P�v�f������̐ߓ_��
C     NPOIN   : �ߓ_��
C     NSIZE   : �X�J�C���C���}�g���N�X�̃T�C�Y
C     NWSUM(ITOTV)�F�S�̍����}�g���N�X�A��P�`ITOTV�s�̑�P��[���v�f����Ίp�v�f�܂ł̗v�f���̗ݘa
C     IFFIX(ITOTV)�F�S���Ɋւ�����i0:���R�A1:�S��)�AITOTV�F���בւ��O�̑S�̎��R�x�ԍ�
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
C *** IEVAB:�v�f�����}�g���N�X�̍s�ԍ�
C     JEVAB:�v�f�����}�g���N�X�̗�ԍ�
C     ITOTV:���׊�����̑S�̍����}�g���N�X�̍s�ԍ�
C     JTOTV:���׊�����̑S�̍����}�g���N�X�̗�ԍ�
C     IITTV:���׊����O�̑S�̍����}�g���N�X�̍s�ԍ�
C     JJTTV:���׊����O�̑S�̍����}�g���N�X�̗�ԍ�
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
C *** �W���}�g���b�N�X����щE�Ӄx�N�g���ɑ΂��āA�􉽊w�I�ȏ�������������
C     ELOAD(ITOTV)�F�c���̓x�N�g��
C     IFFIX(ITOTV)�F�S���Ɋւ�����i0:���R�A1:�S��)
C     MSIZE       �F�X�J�C���C����}�g���N�X�̃T�C�Y�ő�l
C     NDOFN       �F�P�ߓ_����̎��R�x��
C     NFREE       �F�����R�x��
C     NOFIX(IVFIX)�F�S���ߓ_�̔ԍ�
C     NTOTV       �F=NPOIN*NDOFN, ��MTOTV
C     NVFIX       �F�S���ߓ_�̑���
C     NWSUM(ITOTV):�S�̍����}�g���N�X�A��P�`ITOTV�s�̑�P��[���v�f����Ίp�v�f�܂ł̗v�f���̗ݘa
C     PRESC(IVFIX,IDOFN):IVFIX�Ԗڂ̍S���ߓ_��IDOFN�Ԗڂ̎��R�x�̍S����
C     TSTIF(JPOSI): �S�̍����}�g���N�X�i�X�J�C���C���}�g���N�X�j
C *L  ET(IMTOTV)  �F��Ɨp�z��
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
C *** ITOTV:���׊����O�̍s�ԍ��AJTOTV:���׊�����̍s�ԍ�
C
	DO 100 ITOTV=1,NTOTV
	  JTOTV=NORDR(ITOTV)
        ET(JTOTV)=ELOAD(ITOTV)                                          �c���̓x�N�g���̕����ɒ��ӁI
  100 CONTINUE
C
C *** ELOAD(NFREE+1�`NTOTVC)�ɋ����ψʗʂ����Ă����B
C     SUBROUTINE SKYLINE�ɂ�ELOAD(1�`NFREE)�Ɋe���R�x�̕ψʂ�����B
C
      DO 200 IVFIX=1,NVFIX
	  IPOIN=NOFIX(IVFIX)
	  DO 200 IDOFN=1,NDOFN
	    ITOTV=(IPOIN-1)*NDOFN+IDOFN
	    JTOTV=NORDR(ITOTV)
	    IF(IFFIX(ITOTV).NE.0)ET(JTOTV)=PRESC(IVFIX,IDOFN)
  200 CONTINUE
C
C *** ITOTV:���׊�����̍s�ԍ����ő����׏dELOAD���쐬����B
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
C *** �e��ψʃx�N�g���̏��������ɖ߂�
C     ASDIS(ITOTV)�F�����ψ�
C     ELOAD(ITOTV)�F�c���̓x�N�g��
C     NORDR(ITOTV): �S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ����ĕ������鎞�̐V�����S�̎��R�x�ԍ�
C     NTOTV       �F=NPOIN*NDOFN, ��MTOTV
C  *L ASD(IMTOTV) �F��Ɨp�z��
C
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASDIS(MTOTV)
      DIMENSION ELOAD(MTOTV)
      DIMENSION NORDR(MTOTV)
      IF(MCHK.EQ.1)WRITE(8,*) ' *********** RARNG **********'
C
C *** ITOTV:���׊����O�̍s�ԍ��AJTOTV:���׊�����̍s�ԍ�
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
C *** �X�J�C���C���@�ɂ��A���ꎟ�������̉������߂�
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
! *** INDEX=1�F�t�}�g���N�X�����߂�
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
