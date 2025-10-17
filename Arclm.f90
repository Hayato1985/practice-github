! ==================================================================================================
! �v���O�����̍\��
! ARCLM.F90  �F�ʒ������@(ARCH LENGTH METHOD)�̓K�p
!   ARCLM �ʒ������@(ARCH LENGTH METHOD)�̓K�p
!   SCPDT ���ς̌v�Z
! ======================================================================
      SUBROUTINE ARCLM(INDEX,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,&
                       IFFIX,IITER,IINCS,MTOTV,NDOFN,NOFIX,NORDR,&
                       NTOTV,NVFIX,PRESC,RAMDA,RLOAD,RPRVS,RSCLE,TDISP,&
                       TFACT,TLOAD)
! ======================================================================
! *** �ʒ������@(ARCH LENGTH METHOD)�̓K�p
! INDEX=1�F�A���ꎟ�������̉�@(SKYLNE)�̑O�̏���
!      =2�F�A���ꎟ�������̉�@(SKYLNE)�̌�̏���
!       TFACT       �F��[k,j]�A���׏d�W��
!       RAMDA       �F��[k,(j)]�A��IINCS(=k)�񑝕��ɂ����׏d�W���̕ω���
!       FACTO       �F=DRAMD�A1��̔����v�Z�ɂ����׏d�W���̑���
!       DSCLE       �F�ts�A�ψʂ̃X�P�[�����O�W��
!       RSCLE       �F��s�A�ɂ̃X�P�[�����O�W��
!  *L   DRAMD       �F����[k,j]�A��IITER(=j)�񔽕��ɂ����׏d�W���̕ω���
!
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ADIS1(MTOTV)				! ADIS1(ITOTV)�F���t'[k,j1]�ATLOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
      DIMENSION ADIS2(MTOTV)				! ADIS2(ITOTV)�F���t'[k,j2]�AELOAD�ɑ΂���X�P�[�����O���ꂽ�ψʗ�
      DIMENSION ASDIS(MTOTV)				! ASDIS(ITOTV)�F�t[k,(j)]�A��IINCS�񑝕��ɂ��ψʗ�
      DIMENSION ELOAD(MTOTV)				! ELOAD(ITOTV)�F�c���̓x�N�g��
      DIMENSION IFFIX(MTOTV)
      DIMENSION NOFIX(NVFIX)
      DIMENSION NORDR(MTOTV)				! NORDR(ITOTV)�F���̑S�̎��R�x�ԍ�ITOTV�����בւ���̑S�̎��R�x�ԍ�NORDR(ITOTV)
      DIMENSION PRESC(NVFIX,NDOFN)
      DIMENSION RPRVS(MTOTV+1)			    ! RPRVS(ITOTV+1)�F��(IITER-1)�񔽕���̂��x�N�g��
      DIMENSION RLOAD(MTOTV)
      DIMENSION TDISP(MTOTV)				! TDISP(ITOTV)�F���ψ�
      DIMENSION TLOAD(MTOTV)				! TLOAD(ITOTV)�F�׏d���[�h�x�N�g��
      DOUBLEPRECISION, ALLOCATABLE :: ADIS3(:)
      DOUBLEPRECISION, ALLOCATABLE :: DDISP(:)
      DOUBLEPRECISION, ALLOCATABLE :: RPRSN(:)
      ALLOCATE (ADIS3(MTOTV))				! ADIS3(ITOTV)�F�t'[k,(j)]�A�X�P�[�����O���ꂽ�ψʗ�ASDIS
      ALLOCATE (DDISP(MTOTV))				! DDISP(ITOTV)�F���t[k,j]�A��IITER�񔽕��ɂ��ψʗʁ�ADIS1�ɕۑ������
      ALLOCATE (RPRSN(MTOTV+1))			    ! RPRSN(ITOTV+1)�F��(IITER  )�񔽕���̂��x�N�g��
      IF(MCHK.EQ.1)WRITE(8,*)' *********** ARCLM ( INDEX=',INDEX,') **********'
!
      IF(IITER.EQ.1.AND.IINCS.EQ.1)THEN
        CALL ZEROR1(ADIS2,MTOTV,NTOTV)
        CALL ZEROR1(RPRVS,MTOTV+1,NTOTV+1)
      ENDIF
      CALL ZEROR1(RPRSN,MTOTV+1,NTOTV+1)
!
! *** INDEX=1�̏ꍇ�A�E�Ӄx�N�g��ADIS1,ADIS2�ɕ��׊�����̏���JTOTV�ŉ׏d��������
!
    IF(INDEX.EQ.1)THEN
!
      DO ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        ADIS1(JTOTV)=RLOAD(ITOTV)
        ADIS2(JTOTV)=ELOAD(ITOTV)
      ENDDO
!
! *** �E�Ӄx�N�g��ADIS1,ADIS2�̋��E(NFREE+1�`NTOTVC)�ɋ����ψʗʂ����Ă����B
!     SUBROUTINE SKYLINE�ł̌v�Z�ɂ��E�Ӄx�N�g��ADIS1,ADIS2�Ɋe���R�x���Ƃ̕ψʂ�����B
!
      DO IVFIX=1,NVFIX
        IPOIN=NOFIX(IVFIX)
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          JTOTV=NORDR(ITOTV)
          IF(IFFIX(ITOTV).NE.0)THEN
            ADIS1(JTOTV)=PRESC(IVFIX,IDOFN)
            ADIS2(JTOTV)=PRESC(IVFIX,IDOFN)
          ENDIF
        ENDDO
      ENDDO
!
!      IF(MCHK.EQ.1)THEN
!        CALL WRTLD1('RLOAD','ADIS1',RLOAD,MTOTV,NTOTV/2,1,NDOFN)
!        CALL WRTLD1('ELOAD','ADIS2',ELOAD,MTOTV,NTOTV/2,1,NDOFN)
!        CALL WRTLD1('ADIS1','LOAD1',ADIS1,MTOTV,NTOTV,1,NDOFN)
!        CALL WRTLD1('ADIS2','LOAD2',ADIS2,MTOTV,NTOTV,1,NDOFN)
!      ENDIF
!
      RETURN
!
! *** INDEX=2�̏ꍇ�A�׏d���萔����DRAMDA���v�Z���A�����ψʂ����߂�B
!
    ELSEIF(INDEX.EQ.2)THEN
!
! *** �X�P�[�����O�W���̐ݒ�
!
      IF(IINCS.EQ.1.AND.IITER.EQ.1)THEN
        SCPDT1=SCPDT(ADIS1,ADIS1,NTOTV,MTOTV)
        DSCLE=SQRT(SCPDT1)
!        DSCLE=1.0
        RSCLE=1.0
        IF(DSCLE.EQ.0.0)THEN
          WRITE(8,*)' *** ERROR AT SUB ARCLM (DSCLE=0.0) ***'
          STOP ' *** ERROR AT SUB ARCLM (DSCLE=0.0) ***'
        ENDIF
      ENDIF
!
! *** �ψʂ̃X�P�[�����O
!
      DO ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        ADIS1(JTOTV)=ADIS1(JTOTV)/DSCLE
        ADIS2(JTOTV)=ADIS2(JTOTV)/DSCLE
        ADIS3(JTOTV)=ASDIS(ITOTV)/DSCLE
      ENDDO
!
      RAMDAPRV=RPRVS(NTOTV+1)                       ! RAMDAPRV�F�O��̔����v�Z��(IITER-1)�̃�
!
! *** ��P�񔽕��iIITER=1�j�ɂ����郢�t[k,1]�A����[k,1]�̌v�Z
!
      IF(IITER.EQ.1)THEN
        SCPDT1=SCPDT(ADIS1,ADIS1,NTOTV,MTOTV)
        DRAMD=ARCLG/SQRT(SCPDT1+1.0/RSCLE**2)
        DO ITOTV=1,NTOTV
          RPRSN(ITOTV)=DRAMD*ADIS1(ITOTV)
        ENDDO
        RPRSN(NTOTV+1)=DRAMD
!
        IF(IINCS.GE.2)THEN
          SCPDT3=SCPDT(RPRVS,RPRSN,NTOTV+1,MTOTV+1)   ! RPRVS��RPRSN�̓��ς����i�Ȃ��p���s�p�j
          IF(SCPDT3.LT.0.0)THEN                       ! �ɂȂ�悤�Ƀɂ̕�����I������
            DO ITOTV=1,NTOTV+1
              RPRSN(ITOTV)=-RPRSN(ITOTV)
            ENDDO
          ENDIF
        ENDIF
!
! *** ��Q��ȍ~�����iIITER��2�j�ɂ����郢�t[k,j]�A����[k,j]�̌v�Z
!
      ELSE
        SCPDT1=SCPDT(ADIS3,ADIS1,NTOTV,MTOTV)
        SCPDT2=SCPDT(ADIS3,ADIS2,NTOTV,MTOTV)
        DRAMD=-SCPDT2/(SCPDT1+RAMDA/RSCLE)
!
        DO ITOTV=1,NTOTV
          RPRSN(ITOTV)=RPRVS(ITOTV)+DRAMD*ADIS1(ITOTV)+ADIS2(ITOTV)
        ENDDO
        RPRSN(NTOTV+1)=RPRVS(NTOTV+1)+DRAMD
!
      ENDIF
!
! *** RADIAL RETURN�ɂ�钴���ʏ�ւ̈����߂�
!
      SCPDT1=SCPDT(RPRSN,RPRSN,NTOTV+1,MTOTV+1)
      RNORM=SQRT(SCPDT1)
      IF(RNORM.EQ.0.0)STOP ' RNORM=0 '
      DO ITOTV=1,NTOTV+1
        RPRSN(ITOTV)=ARCLG*RPRSN(ITOTV)/RNORM
      ENDDO
!
      IF(IITER.EQ.1)THEN
        DRAMD=RPRSN(NTOTV+1)
      ELSE
        DRAMD=RPRSN(NTOTV+1)-RAMDAPRV
      ENDIF
!
      DO ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        DDISP(ITOTV)=RPRSN(JTOTV)-ASDIS(ITOTV)
        ASDIS(ITOTV)=RPRSN(JTOTV)
        RPRVS(JTOTV)=RPRSN(JTOTV)                       ! ����̌ʒ������̂��߂̏���
      ENDDO
      RPRVS(NTOTV+1)=RPRSN(NTOTV+1)                     ! ����̌ʒ������̂��߂̏���
!
! *** ���׏d�W��TFACT�A���ψ�TDISP�A���׏dTLOAD�̐ݒ�B
!
      TFACT=TFACT+DRAMD
      DO ITOTV=1,NTOTV
! >>>>>>>CMEM�̏ꍇ          TDISP(ITOTV)=TDISP(ITOTV)+DDISP(ITOTV)
        TLOAD(ITOTV)=RLOAD(ITOTV)*TFACT
        ADIS1(ITOTV)=DDISP(ITOTV)                       ! ����̔����i��IITER��ځj�ɂ��ψʗʂ�ADIS1�ɓ����
      ENDDO
!
      IF(MCHK.EQ.1)THEN
        WRITE(8,9001)IITER,RAMDA,DRAMD,ARCLG,RNORM,TFACT,SCPDT3
 9001   FORMAT(' *** IITER',6X,'RAMDA',6X,'DRAMD',6X,'ARCLG',6X,'RNORM',6X,'TFACT',6X,'SCPDT3',/,7X,I3,6E11.4)
!        CALL WMTR1(NDOFN,ADIS1,MTOTV,0,NTOTV,0,'ADIS1','ARCLM ',0)
!        CALL WMTR1(NDOFN,ADIS2,MTOTV,0,NTOTV,0,'ADIS2','ARCLM ',0)
!        CALL WMTR1(NDOFN,ADIS3,MTOTV,0,NTOTV,0,'ADIS3','ARCLM ',0)
!        CALL WMTR1(NDOFN,TDISP,MTOTV,0,NTOTV,0,'TDISP','ARCLM ',0)
!        CALL WMTR1(NDOFN,DDISP,MTOTV,0,NTOTV/2,0,'DDISP','ARCLM ',0)
!        CALL WMTR1(NDOFN,RPRVS,MTOTV+1,0,NTOTV+1,0,'RPRVS','ARCLM ',0)
!        CALL WMTR1(NDOFN,RPRSN,MTOTV+1,0,NTOTV+1,0,'RPRSN','ARCLM ',0)
!        CALL WRTLD1('TLOAD','ARCLM',TLOAD,MTOTV,NTOTV,2,NDOFN)
      ENDIF
    ENDIF   ! �i�����܂�INDEX=2�̏ꍇ�j
!
	DEALLOCATE(ADIS3,DDISP,RPRSN)	! ���I�z����J������
!
    RETURN
    END
!
!
!
! ======================================================================
      FUNCTION SCPDT(A,B,N,M)
! ======================================================================
! *** ���ς̌v�Z
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(M),B(M)
      SCPDT=0.0
      DO 100 I=1,N
        SCPDT=SCPDT+A(I)*B(I)
  100 CONTINUE
      RETURN
      END