! =====================================================================
      SUBROUTINE FRAME(INDEX,CALFR,COORD,VALK,&
                       NPOIN,NDOFN,NELEF,NODEF,NPIN,NPFRM,  &
                       NVARF,PROPF,STRSF,STRF1,                   &
                       BETA,IFFIX,MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
! =====================================================================
! *** ���Ȃ��v�f�̓��o�́A�����}�g���N�X�쐬�Ɋւ���T�u���[�`���� ***
!     NELEF�F�Ȃ��v�f��
!     NPFRM�F�������́i=0�F�Ȃ��A=1�F����j
!     NVARF�F�Ȃ��v�f�̎�ސ�
!
      PARAMETER( MFORMAT=1 )    ! MFORMAT�F�Œ�t�H�[�}�b�g�i=0�j�A�z��f�[�^�̈ꕔ�̓t���[�t�H�[�}�b�g�i=1�j
!
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION CALFR(NELEF)              ! CALFR(IELEF)�F���ނ̒���
      DIMENSION COORD(NPOIN,3)            ! COORD(IPOIN,IDIR)�F�ߓ_���W
      DIMENSION VALK(NELEF,2)             ! VALK(IELEF, )�F�ޒ[1,2�̐ڍ������i��1.0E+20�F���ځA��0.0�F�s���ځA���̑��F�����ځj
      DIMENSION NODEF(NELEF,3)            ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION NPIN(NELEF)               ! NPIN(IELEF)�F���ޒ[���̐ڍ������i=0�F�����A=1�F�s�����A=2�F���s���A=3�F�s���s���j
      DIMENSION STRSF(NELEF,NDOFN*2)      ! STRSF(IELEF,IDOFN)�F�Ȃ��v�f�̉��́i����W�n�ɂ����镔�ޒ[�׏d)
      DIMENSION STRF1(NELEF,NDOFN*2)      ! STRSF(IELEF,IDOFN)�F�Ȃ��v�f�̉��́i���ލ��W�n�ɂ����镔�ޒ[�׏d)
      DIMENSION BETA(NELEF)               ! BETA(IELEF)�F�Ȃ����ނ̉�]�p
      DIMENSION IFFIX(MTOTV)              ! IFFIX(ITOTV)�F�S���Ɋւ�����i0:���R�A1:�S��)�AITOTV�F���בւ��O�̑S�̎��R�x�ԍ�
      DIMENSION NORDR(MTOTV)              ! NORDR(ITOTV): �S�̍����}�g���N�X[K]�̑S�̎��R�x�ԍ�ITOTV��FREE��FIX�ɕ����ĕ������鎞�̐V�����S�̎��R�x�ԍ�
      DIMENSION NWSUM(0:MTOTV)            ! NWSUM(ITOTV):�S�̍����}�g���N�X�A��P�`ITOTV�s�̑�P��[���v�f����Ίp�v�f�܂ł̗v�f���̗ݘa
      DIMENSION PROPF(NVARF,10)           ! �Ȃ��v�f�̍ޗ������i��ޔԍ��A��������(1�`10)�j
      DIMENSION TSTIF(MSIZE)              ! TSTIF(JPOSI): �S�̍����}�g���N�X�i�X�J�C���C���}�g���N�X�j
      DIMENSION DIS(NPOIN*NDOFN)          ! DIS(IPOIN*IDOFN)�F�ψʑ���
      DIMENSION Q(NPOIN*NDOFN)            ! Q(ITOTV)�F�c���̓x�N�g��
      DIMENSION C(12,6),CSK(12,6),CT(6,12),SK(6,6),SKP(12,12)
      DIMENSION C0(12,6),CSK0(12,6),SKP0(12,12)
      DIMENSION SKG(12,12),SKGT(12,12),TSKGT(12,12)
      DIMENSION TR(12,12),TRT(12,12),SKK(12,12)
      DIMENSION DM(12,1),PM(12,1)
      DIMENSION VECLY(3)                                                ! VECLY(3)�F�Ǐ����W�n�̂���������\�킷�x�N�g��
character*170 text

      NEVAB=NDOFN*2
!
! ---------------------------------------------------------------------
      IF(INDEX.EQ.1)THEN ! *** ���Ȃ��v�f�Ɋւ���f�[�^���t�@�C�����V����Ǎ��ށ� ***
! ---------------------------------------------------------------------

! *** ���Ȃ��v�f�̏d�ʁA�����A����f�␳�W���A���ގ厲�̖@���x�N�g����

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

! *** ���Ȃ��v�f�̐ߓ_�ԍ��A�ޒ[�o�l�̌W����

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
!!!!!!�P�ʌn�ύX  WEIGHT=PROPF(IV,7)*CALFR(IE)  !�d��[kgf]���P�ʏd��[kg/m]*����[m]
  WEIGHT=PROPF(IV,7)*CALFR(IE)*9.8  !�d��[N]���P�ʏd��[kg/m]*����[m]*9.8
  WSUM=WSUM+WEIGHT
  write(8,8004)IE,NODEF(IE,1),NODEF(IE,2),VALK(IE,1),VALK(IE,2),NODEF(IE,3),CALFR(IE),WEIGHT
enddo
8003 format(//&
'�Ȃ��v�f',/&
'       �v�f�ԍ�       �ߓ_�ԍ�                �ޒ[�̃o�l�W��                  ���ގ�ޔԍ�       ���ޒ���           �d��',/&
'                           I�[            J�[            I�[            J�[                           [m]          [kgf]',/&
'     I=1�`NELEF     NODEF(I,1)     NODEF(I,2)      VALK(I,1)      VALK(I,2)     NODEF(I,3)')
8004 format(3I15,2G15.4,I15,2F15.4)
WRITE(8,8005) WSUM
8005 FORMAT(4X,120('-'),/,85X,'  �Ȃ����ނ̑��d�� (',F15.4,' [N])')

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
      ELSEIF(INDEX.EQ.2)THEN ! *** ���Ȃ��v�f�̍����}�g���N�X���쐬���遄 ***
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
!  *** �����ލ��W�n�ɂ���������ψʍ����}�g���N�X�̍쐬�� ***
!
          CALL STIFF(SK,EA,EIz,EIy,GIp,GA,CAL1,VALKAPPA,VALK(NE,1),VALK(NE,2))   
!         CALL WMTR2(2,SK,6,6,6,6,'SK   ','NE   ',NE)
!
!  *** ���}�g���N�X�b0�C�bt�C���W�ϊ��}�g���N�X�s�q�̍쐬�� ***
!
          DO IDIR=1,3
            VECLY(IDIR)=PROPF(IV,IDIR+7)    ! �厲�@���x�N�g��
          ENDDO
          CALL TRANST(TRT,CAL1,CAL2,DELX,DELY,DELZ,BETA(NE),VECLY)
          CALL TRANS(12,12,TRT,TR)
          CALL CNECT1(C0,CT,CAL1,CAL2,DELX,DELY,DELZ,BETA(NE),TRT)
          CALL TRANS(6,12,CT,C)
!
!  *** ���e�������}�g���N�X�r�j�o�̍쐬�� ***
!
          CALL MLTPLY(12,6,6,C,SK,CSK)
          CALL MLTPLY(12,6,12,CSK,CT,SKP)
!
!  *** �����ލ��W�n�ɂ����镔�ޒ[�׏d�Z��̂��߂̃}�g���N�X�r�j�o�O�̍쐬�� ***
!
          CALL MLTPLY(12,6,6,C0,SK,CSK0)
          CALL MLTPLY(12,6,12,CSK0,CT,SKP0)
!
! *** ���Ȃ��v�f�̕��މ��͂��v�Z���遄 ***
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
!  *** �����ލ��W�n�ɂ�����􉽍����}�g���N�X�̍쐬�� ***
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
! *** ���v�f�����}�g���b�N�X���S�̍����}�g���b�N�X�� ***
!
          CALL MATRIX(SKK,NE,NODEF,NELEF,MSIZE,NDOFN,6,2,NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
!
! *** ���Ȃ����ނ̕��ޒ[�׏d�x�N�g���𕔍ލ��W�n�r�s�q�e�q�������W�n�r�s�q�e�P�ɕϊ����遄 ***
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
! *** ���Ȃ����މ��͂̓����ߓ_�̓x�N�g��DM��S�̂̐ߓ_�̓x�N�g���p�ɉ����遄 ***
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
      ELSEIF(INDEX.EQ.3)THEN ! *** ���Ȃ��v�f�̕��މ��͂��t�@�C�����W�ɏo�͂��遄 ***
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
      DIMENSION H(6,6)  ! �}�g���N�X�g
      DIMENSION F(6,6)  ! �}�g���N�X��
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
! *** �Ȃ��v�f�̊􉽍����}�g���N�X���쐬����T�u���[�`���B
!     �e�����́u��������́v�i��䒉�F���E�|���فjPP.91�`96���Q�l�ɂ��Ă���B
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION SKG(12,12)                            ! SKG(12,12)�F�􉽍����}�g���N�X
      DIMENSION SKGF(4,4)                             ! SKMNF(4,4)�F��Gf��(�`-1)t�E��Gf��(�`-1)t�E��Gf�E�`-1
      DIMENSION AINV(4,4)                             ! AINV(4,4) �F�`-1
      DIMENSION FK(0:6)                               ! FK(6)     �F��k����f(z)z^kdz
      DIMENSION CL(0:6)                               ! CL(6)     �FCL(n)=CAL1**n
      DIMENSION STRSM(NELEM,NDOFN*2)
      INTEGER NORDR(12)                               ! NORDR(12) �FSKGF(4,4)��SKG(12,12) �ɑ�����鎞�̑����
      DATA NORDR /2,6,8,12,3,5,9,11,4,0,10,0/
!
!     NDOFN �F�P�ߓ_����̎��R�x��(6)
!     CAL1  �F���ނ̒���
!     FZ    �F���ޒ[�׏d�iPx,Py,Pz,Myz,Mzx,Mxy�j
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
! *** �b0�C�b�s�}�g���N�X�̍쐬
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION C0(12,6)
      DIMENSION CT(6,12)
      DIMENSION TRT(12,12)    ! TRT(12,12)�F���W�ϊ��}�g���N�X�s�q�̓]�u
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
! *** ���W�ϊ��}�g���N�X�̍쐬
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TRT(12,12)    ! TRT(12,12)�F���W�ϊ��}�g���N�X�s�̓]�u
      DIMENSION VECLX(3)      ! VECLX(3)�F�Ǐ����W�n�̂���������\�킷�x�N�g��
      DIMENSION VECLY(3)      ! VECLY(3)�F�Ǐ����W�n�̂���������\�킷�x�N�g��
      DIMENSION VECLZ(3)      ! VECLZ(3)�F�Ǐ����W�n�̂���������\�킷�x�N�g��
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
! *** �x�N�g���`�A�a�̊O�ςb�����߂�B
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3),B(3),C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      RETURN
      END
