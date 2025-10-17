! =====================================================================
SUBROUTINE INPUT(ANG,COORD,DDV,DISMAX,ERR,FCOE,IFPRE,LODFIX, &
                 MAXCYL,MEIGN,MSTRE,NANG,NANGLE,NCONC,NDOFN,NELEM,NFPOIN,&
                 NLOAD,NODEM,NOFIX,NONISO,NPOIN,NPRET,NPRINT,NPROB,NREPT,  &
                 NSTEP,NTIMES,NVFIX,NVARM,PMGNF,PPRES,PRESC,   &
                 PRETEN,PROPM,SCF,SSNOW,STRSM,WCF,WIND,WWIND,MINCR)
! =====================================================================
      PARAMETER( MFORMAT=1 )    ! MFORMAT�F�Œ�t�H�[�}�b�g�i=0�j�A�z��f�[�^�̈ꕔ�̓t���[�t�H�[�}�b�g�i=1�j
! 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ANG(NELEM)
      DIMENSION COORD(NPOIN,3)
      DIMENSION FCOE(NPOIN*NDOFN)
      DIMENSION IFPRE(NVFIX)
      DIMENSION NANG(NELEM)
      DIMENSION NODEM(NELEM,4)          ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NOFIX(NVFIX)
      DIMENSION PMGNF(NDOFN,NTIMES)
      DIMENSION PPRES(NTIMES)
      DIMENSION PRESC(NVFIX,NDOFN)
      DIMENSION SCF(NELEM)
      DIMENSION SSNOW(NTIMES)
      DIMENSION STRSM(NELEM,MSTRE)
      DIMENSION WWIND(NTIMES)
      DIMENSION WCF(NELEM)

DIMENSION PROPM(NVARM,7,3)        !���v�f�̍ޗ������i��ޔԍ��A��������(1�`7)�A�i�KTri-Linear(1�`3)�j
character*150 text
!
! *** ���׏d�{���A�����A���x���A��׏d�̐ݒ聄 ***
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

! *** ���ߓ_���W�A���ڐߓ_�׏d�̐ݒ聄 ***

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

! *** ���S�������̐ݒ聄 ***

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

! *** ���O�p�`���v�f�̕��ޓ����̐ݒ聄 ***

if(NVARM.NE.0)then
  do i=1,6
    read(7,'(A)')text
    write(8,'(A)')text
  enddo
  do i = 1,NVARM
!PROPM(NVARM,7,3)   !���v�f�̍ޗ������i��ޔԍ�IVARM�A��������(1�`7)�A�i�KTri-Linear(1�`3)�j
    read(7,*,err=999)text,(PROPM(i,j,1),j=1,7)
    do k=2,3
      read(7,*,err=999)text,(PROPM(i,j,k),j=1,5)
    enddo
    write(8,168)i,(PROPM(i,j,1),j=1,7),((PROPM(i,j,k),j=1,5),k=2,3)
  enddo
endif
168 format(&
I15,'     ��P���z  ',7F15.4,/&
15X,'     ��Q���z  ',5F15.4,/&
15X,'     ��R���z  ',5F15.4)

! *** ���O�p�`���v�f�̐ߓ_�ԍ��A�������́A�׏d�W���̐ݒ聄 ***

!NSURF=0   !NSURF=0�F����׏d�Ȃ��A=1:���A=2�F��A=3:�G���[
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
    901 format(///,'error�@���׏d�Ɛ�׏d�𓯎��Ɍv�Z���邱�Ƃ͂ł��܂���',/&
                  '       WWIND��SSNOW�̂ǂ��炩����͑S�ă[���ɂ��Ă�������')
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
      !@@@@@@ GRANPA DOME �̕��͌W���̌v�Z @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if(NSURF.eq.1)then
        NWCOE=0  ! =1 �ȉ��̃v���O�����ŕ��͌W����ݒ肷�� =0 �������Ȃ��i���̓f�[�^�̕��͌W�����g���Čv�Z����j
        if(NWCOE.eq.1)then
          center=0.0
          span=13.6*2.0
          do j=1,3
            center=center+COORD(NODEM(i,j),1)/3.0
          enddo
          if(center.LE.-span/8.0)then  !���� �[���`1/8
            WCF(i)=-0.25
          elseif(center.LE.-span/4.0)then  !����1/8�`1/4
            WCF(i)=-0.25
          elseif(center.LE.0.0)then   !����1/4�`����
            WCF(i)=-0.61
          elseif(center.LE.span/4.0)then  !�����`����1/4
            WCF(i)=-0.61
          elseif(center.LE.span/8.0)then  !����1/4�`1/8
            WCF(i)=-0.01
          else  !���� 1/8�`�[��
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
'���v�f                                                                            ��������                                     �׏d�W��',/&
'          �ԍ�               I              J              K       ��ޔԍ�     IJ�������� IJ�����������       ����f��     �i��or��j',/&
'                                                                                     [N/m]          [N/m]          [N/m]',/&
'     I=1�`NELEM     NODEM(I,1)     NODEM(I,2)     NODEM(I,3)                    STRSM(I,1)     STRSM(I,2)     STRSM(I,3)')
136 format(5I15,3F15.4)
137 format(5I15,4F15.4)

!!
!! *** �����v�f��IJ�����ƍޗ��̎厲�����̂Ȃ��p�x�̐ݒ聄 ***
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
! *** ���P�[�u���܂��̓g���X�Ɋւ���f�[�^���t�@�C�����V����Ǎ��ށ� ***
!     NELEC�F�P�[�u���܂��̓g���X�̗v�f��
!     NVARC�F���ސ�
!     NPRSTR�F�P�[�u���܂��̓g���X�̏������́i=0�F�Ȃ��A=1�F����j
!     DV�F�P�[�u�����ɂ񂾏ꍇ�ɁA������EA/DV�Ƃ���B�ʏ��1.0�Ƃ���B
!
      PARAMETER( MFORMAT=1 )    ! MFORMAT�F�Œ�t�H�[�}�b�g�i=0�j�A�z��f�[�^�̈ꕔ�̓t���[�t�H�[�}�b�g�i=1�j
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION BL(NELEC)        !   BL(IELEC)�F���݂̒���
      DIMENSION CAL(NELEC)       !   CAL(IELEC)�F���Ђ��ݎ��̒���
      DIMENSION STRSC(NELEC)     !   STRSC( )�F��������
      DIMENSION STRC1(NELEC)
      dimension PROPC(NVARC,3)   !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
      DIMENSION CTEMP1(NTIMES)   ! CTEMP1(LOADCREM)�F�P�[�u���ɗ^���鉷�x����[��]�AMCTEMP=1�i�P�[�u���̉��x���͉�͂��s���j�̏ꍇ
      DIMENSION CTEMP2(2,NELEC)  ! CTEMP2(1,IELEC)�F���x���͂�^����(=1)�A�^���Ȃ�(=0)�ACTEMP2(2,IELEC)�F���x�Ђ���
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
character*150 text

! *** ���P�[�u���܂��̓g���X�̕��ޔԍ��A�f�ʐρA���x�A�����O����

do i=1,6
  read(7,'(A)')text
  write(8,'(A)')text
enddo
read(7,*,err=999) (num,PROPC(i,1),PROPC(i,2),PROPC(i,3),i=1,NVARC)
write(8,101) (i,PROPC(i,1),PROPC(i,2),PROPC(i,3),i=1,NVARC)
101 format(I15,2F15.4,G15.4)

! *** ���P�[�u���܂��̓g���X�v�f�̐ߓ_�ԍ��A���ޔԍ��A�������́�

read(7,'(A)')(text,i=1,6)
read(7,*,err=999)(num,(NODEC(i,j),j=1,3),STRSC(i),i=1,NELEC)

! *** ���P�[�u���܂��̓g���X�̗v�f�����i�����A�d�ʁA�����j�̌v�Z�Əo�́�

call ZEROR1(STRC1,NELEC,NELEC)
WSUM=0.0
write(8,102)
DO NE=1,NELEC
  STRC1(NE)=STRSC(NE)
  I=NODEC(NE,1)
  J=NODEC(NE,2)
  IV=NODEC(NE,3)
  BL(NE)=SQRT((COORD(I,1)-COORD(J,1))**2+(COORD(I,2)-COORD(J,2))**2+(COORD(I,3)-COORD(J,3))**2)
  EA=PROPC(IV,1)*DABS(PROPC(IV,2))   !EA[N]=�����O��[N/mm2]*�f�ʐ�[mm2]
  IF((STRSC(NE).EQ.0.0).OR.(EA.EQ.0.0))THEN
    CAL(NE)=BL(NE)
  ELSE
    CAL(NE)=BL(NE)/(STRSC(NE)/EA+1.0)
  ENDIF
!!!!!�P�ʌn�ύX  W=PROPC(IV,3)*CAL(NE)   !�d��[kgf]���P�ʏd��[kg/m]*����[m]
  W=PROPC(IV,3)*CAL(NE)*9.8   !�d��[N]���P�ʏd��[kg/mm2/m]*����[m]
  WSUM=WSUM+W
  write(8,103)NE,(NODEC(NE,i),i=1,3),STRSC(NE),PROPC(IV,2),BL(NE),W
ENDDO
102 format(//&
'�P�[�u���v�f',/&
'       �v�f�ԍ�       �ߓ_�ԍ�                      ��ޔԍ�       ��������         �f�ʐ�       ���ޒ���           �d��',/&
'                           I�[            J�[                           [N]          [cm2]            [m]            [N]',/&
'      I=1�`NELEC     NODEC(I,1)     NODEC(I,2)     NODEC(I,3)         STRSC(I)     PROPC(I,2)                              ')
103 format(4I15,4F15.4)
WRITE(8,337) WSUM
337 FORMAT(4X,120('-'),/,85X,'  �P�[�u���̑��d�� (',F15.4,' [N]')
!
! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
      IF(MCTEMP.EQ.1)THEN
        READ(7,*,err=999) (I,CTEMP1(ITIMES),ITIMES=1,NTIMES)    ! CTEMP1(LOADCREM)�F�P�[�u���ɗ^���鉷�x����[��]�AMCTEMP=1�i�P�[�u���̉��x���͉�͂��s���j�̏ꍇ
        WRITE(8,410) (ITIMES,CTEMP1(ITIMES),ITIMES=1,NTIMES)
  410   FORMAT(/,5X,' TEMPERATURE INCREMENT OF CABLE ELEMENTS',/,'  LOADCREM TEMPERATURE [DEGREE]',/,(I10,F10.4))
        READ(7,*,err=999) NCTEMP    ! NCTEMP�F�P�[�u���̉��x���͌v�Z�̑ΏۂƂ���v�f�̐�
        IF(NCTEMP.LT.NELEC)THEN
          READ(7,*,err=999) (IELEC,CTEMP2(1,IELEC),ICTEMP=1,NCTEMP)   ! CTEMP2(1,IELEC)�F���x���͂�^����(=1)�A�^���Ȃ�(=0)�ACTEMP2(2,IELEC)�F���x�Ђ���
        ELSE
          DO IELEC=1,NELEC    ! NCTEMP=NELEC�Ƃ����ꍇ�ɂ́ACTEMP2(1,IELEC)����͂��Ȃ��Ă������I��1.0�ɐݒ肷��B
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
! *** �v�Z�ߒ�����ʁA�t�@�C���ɏo�͂���
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
! *** �ő�ψ�DD�ƍő�c����QQ�̌v�Z
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
DIMENSION STRSF(NELEF,NDOFN*2)                                    ! STRSF(IELEF,IEVAB)�F�Ȃ��v�f�̉��́i���ލ��W�n�ɂ����镔�ޒ[�׏d)
DIMENSION NODEF(NELEF,3)            ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
DIMENSION FFF(NPOIN*NDOFN)
DIMENSION BL(NELEC)
DIMENSION CAL(NELEC)
DIMENSION STRSC(NELEC)
DIMENSION STRC1(NELEC)
DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
dimension PROPC(NVARC,3)    !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj
DIMENSION COORD(NPOIN,3)
DIMENSION COOR0(NPOIN,3)
DIMENSION FF(NPOIN*NDOFN)
DIMENSION Q(NPOIN*NDOFN)

DIMENSION NODEM(NELEM,4)          ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
DIMENSION PROPM(NVARM,7,3)        !���v�f�̍ޗ������i��ޔԍ��A��������(1�`7)�A�i�KTri-Linear(1�`3)�j

! >>>>>>>>>> �P�[�u���̉��x���͉�͐�p�̐ݒ�>>>>>>
!     MCTEMP                        ! MCTEMP�F�P�[�u���̉��x���͉�͂��s��Ȃ��i=0�j�A�s���i=1�j
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)�F���x���͂�^����(=1.0)�A�^���Ȃ�(=0.0)�ACTEMP2(2,IELEC)�F���x�Ђ���
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
DIMENSION FALL(6),QALL(6),UNB(6)
DIMENSION PST(3),VALMM(7,2),NVALMM(7,2)
!
! *** �����v�f�̒��͂��t�@�C�����W�ɏo�͂��遄 ***
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
! *** �������͂̍ő�lVALMM(ISTRE,1)�ƍŏ��lVALMM(ISTRE,2)���t�@�C�����W�ɕۑ����遄
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
! ETFE�t�B�����̑������͂��v�Z���C���̍ő�lVALMM(1,1)�ƍŏ��lVALMM(1,2)���t�@�C�����W�ɕۑ�����
CNV=1.0E+00               ! CNV�F�P�ʂ̊��Z�W���i������mm�P�ʂ̏ꍇ1.0�A������m�P�ʂ̏ꍇ1000.0�j
DO IELEM=1,NELEM
  IV=NODEM(IELEM,4)          ! NODEM(IELEM,4)�F���ޔԍ�
  THICK=PROPM(IV,6,1)        ! PROPM(IV,6,1)�F�t�B�����̌���[mm]
  VALUE=1.0/CNV/THICK        ! VALUE�F�P�ʂ�ϊ����邽�߂̌W��
  PX=STRSM(IELEM,1)*VALUE
  PY=STRSM(IELEM,2)*VALUE
  PXY=STRSM(IELEM,3)*VALUE
  F=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)    ! F�F��������[N/mm2]
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
  WRITE(10,128) IELEM,PX,PY,PXY     !�������͂��t�@�C�����P�O�ɏo�͂���
ENDDO
if(NELEM.gt.0)WRITE(8,817)(NVALMM(1,I),VALMM(1,I),I=1,2)
817   FORMAT(/,' (MAXIMUM AND MINIMUM OF EQUIVARENT MEMBRANE STRESS) [N/mm2]',/,'           ELEM         STRESS',/,'  MAX',I10,F15.5,/,'  MIN',I10,F15.5)
!
! *** ���ό`��̉��͂��t�@�C�����P�O�ɕۑ����遄 ***
!
!IF(NELEM.GE.1)WRITE(10,128) (K,(STRSM(K,ISTRE),ISTRE=1,3),K=1,NELEM)
128 FORMAT(I5,3F15.4)
!
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  �ʂ̌v�Z���ʏo��
!      IP=111
!      WRITE(11,8170)IP,COORD(IP,3),NVALMM(1,1),VALMM(1,1),NVALMM(1,2),VALMM(1,2)
! 8170 FORMAT('  Z(',I4,')=',G12.5,'   ST_MAX=(',I4,',',G12.5,')','   ST_MIN=(',I4,',',G12.5,')')
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! *** ���P�[�u���v�f�̕��މ��͂��t�@�C�����W�ɏo�͂��遄 ***
!
      IF(NELEC.NE.0)CALL CABLE3(1,BL,CAL,NELEC,NVARC,NODEC,STRSC,STRC1,PROPC,MCTEMP,CTEMP2)
!
! *** ���Ȃ��v�f�̕��މ��͂��t�@�C�����W�ɏo�͂��遄 ***
!
      IF(NELEF.GT.0)WRITE(8,8040) (NE,(STRSF(NE,IEVAB),IEVAB=1,NDOFN*2),NE=1,NELEF)
 8040 FORMAT(//,' (MEMBER FORCES OF FRAME ELEMENTS IN LOCAL COORDINATES)',/,' ELEM',    &
            12X,'PX1',12X,'PY1',12X,'PZ1',11X,'MYZ1',11X,'MZX1',11X,'MXY1',/,           &
            17X,'PX2',12X,'PY2',12X,'PZ2',11X,'MYZ2',11X,'MZX2',11X,'MXY2',/,           &
            (I5,6E15.4,/,5X,6E15.4))
!
! *** ���Ȃ��v�f�̕��ގ��IV���Ƃɉ��͂̐�Βl���ő�ƂȂ镔�ނ̔ԍ�NVALMM(ISTRE,1)�Ɖ��͂̒lVALMM(ISTRE,1)���t�@�C�����W�ɕۑ����遄
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
!  �ʂ̌v�Z���ʏo��
!		  JE=NVALMM(5,1)
!          WRITE(11,8053)IV,JE,(STRSF(JE,ISTRE),ISTRE=1,6)
! 8053 FORMAT(17X,10X,'PX',10X,'PY',10X,'PZ',9X,'MYZ',9X,'MZX',9X,'MXY',/,&
!             'MAX-MYY(IV',I1,':',I4,'):',6E12.4)
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	  ENDDO
    ENDIF
!
! *** ���ό`��̍��W���t�@�C�����X�ɕۑ����遄 ***
!
      CALL WRTCOD(2,COORD,NODEC,NELEM,NELEC,NODEM,NPOIN,NELEF,NODEF)
!
! *** ���ό`��̃P�[�u�����́A�Ȃ����ޗ͂��t�@�C�����P�O�ɕۑ����遄 ***
!
      IF(NELEC.GE.1)WRITE(10,990) (K,STRSC(K),K=1,NELEC)
  990 FORMAT(4(I5,F15.4,5X))
!
!     WRITE(10,8041) (NE,(STRSF(NE,IEVAB),IEVAB=1,NDOFN*2),NE=1,NELEF)
      WRITE(10,8042) (NE,(STRSF(NE,IEVAB),IEVAB=7,7),NE=1,NELEF)
 8041   FORMAT(I5,12E15.4)
 8042   FORMAT(4(I5,E15.4))
!
! *** ���ό`��̍��W�Ə����`��ɑ΂���ψʂ��t�@�C�����W�ɕۑ����遄 ***
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
! *** ���ߓ_�ψʂ̍ő�lVALMM(IDOFN,1)�ƍŏ��lVALMM(IDOFN,2)���t�@�C�����W�ɕۑ����遄
!
      DO IPOIN=1,NPOIN
        VLENG=0.0   ! VLENG�F�ψʃx�N�g���̒���
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
! *** ���c���́i�s�ލ��́j�A�ߓ_�O�́A�m�������t�@�C�����W�ɕۑ����遄
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
! *** �c���́i�s�ލ��́j�A�ߓ_�O�́A�m�������t�@�C�����W�ɕۑ�����
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FF(NPOIN*NDOFN)
      DIMENSION FFF(NPOIN*NDOFN)
      DIMENSION IFFIX(NPOIN*NDOFN)
      DIMENSION Q(NPOIN*NDOFN)
      DIMENSION FALL(6,2),RALL(6,2)
      CHARACTER*2 NB(24)
!
! *** �c���͂ƌŒ�_����
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
! *** �ό`��̐ߓ_�O��
      WRITE(8,420)(IPOIN,(FF((IPOIN-1)*NDOFN+IDOFN),IDOFN=1,NDOFN),IPOIN=1,NPOIN)
  420 FORMAT(///,' (TOTAL LOAD AT EACH NODAL POINT AFTER DEFORMATION) [N]',//, &
            '  NODE',4X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD',8X,'XX-LOAD',8X,'YY-LOAD',8X,'ZZ-LOAD',/,  &
            (I5,6(F14.4,1X)))
!
! *** �����󋵂̌���
      CALL ZEROR2(FALL,NDOFN,2,NDOFN,2) ! FALL�F�O�͂̑��a
      CALL ZEROR2(RALL,NDOFN,2,NDOFN,2) ! RALL�F�c���͂̑��a
      FNORM=0.0
      RNORM=0.0
	  DO IDOFN=1,NDOFN
        DO IPOIN=1,NPOIN
	      ITOTV=(IPOIN-1)*NDOFN+IDOFN
	      IF(IFFIX(ITOTV).EQ.0)THEN
            FALL(IDOFN,1)=FALL(IDOFN,1)+FFF(ITOTV)  ! FALL(IDOFN,1)�F�O��(���R�ߓ_�j�̑��a
            RALL(IDOFN,1)=RALL(IDOFN,1)+Q(ITOTV)    ! RALL(IDOFN,1)�F�c����(���R�ߓ_�j�̑��a
          ELSE
            FALL(IDOFN,2)=FALL(IDOFN,2)+FFF(ITOTV)  ! FALL(IDOFN,2)�F�O��(�S���ߓ_�j�̑��a
            RALL(IDOFN,2)=RALL(IDOFN,2)+Q(ITOTV)    ! RALL(IDOFN,2)�F�c����(�ߓ_�j�̑��a
          ENDIF
        ENDDO
        FNORM=FNORM+FALL(IDOFN,1)**2
        RNORM=RNORM+RALL(IDOFN,1)**2
      ENDDO
!
      FNORM=SQRT(FNORM)                 ! FNORM�F�O�͂̃m����
      RNORM=SQRT(RNORM)                 ! RNORM�F�c���͂̃m����
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
! *** ���`��̃f�[�^���t�@�C�����X�ɕۑ����遄 ***
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION NODEM(NELEM,4)   ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION NODEF(NELEF,3)   ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
!
      IF(INDEX.EQ.1)THEN    ! INDEX=1�iSUBROUTINE INPUT �����CALL�j�̏ꍇ�ɂ͗v�f�����L�^����
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
! *** �C���v�b�g�f�[�^�ɃG���[���������ꍇ�́A�c��̃f�[�^���o�͂��ďI������B
character*150 text
write(6,*)'***** �C���v�b�g�f�[�^�ɃG���[�����邽�ߎc��̃C���v�b�g�f�[�^���o�͂��ďI�����܂� *****'
write(8,100)
100 format(/////,'***** �c��̃C���v�b�g�f�[�^���ȉ��ɏo�͂��܂� *****',//)
do i=1,10000
  read(7,'(A)',end=999)text
  write(6,'(A)')text
  write(8,'(A)')text
enddo
999 stop
end