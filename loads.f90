! =====================================================================
      SUBROUTINE SETLOAD(NFLAG,NPOIN,FF, PRES, WIND, SNOW,LODFIX,NTIMES,NDOFN,FCOE,PMGNF,PPRES,WWIND,SSNOW)
! =====================================================================
! *** �׏d�i���ڐߓ_���͉׏dFF�A����PRES�A���x��WIND�A�ϐ�׏dSNOW�j�̐ݒ�
!     NFLAG�F��LOADCREM�i�׏d�����v�Z�̏ꍇ�j
      IMPLICIT DOUBLEPRECISION(A-H,O-Z)
      DIMENSION FF(NPOIN*NDOFN)
      DIMENSION FCOE(NPOIN*NDOFN)
      DIMENSION PMGNF(NDOFN,NTIMES)
      DIMENSION PPRES(NTIMES)
      DIMENSION WWIND(NTIMES)
      DIMENSION SSNOW(NTIMES)
!
! *** ���ڐߓ_���͉׏d�̐ݒ�iFF���׏d�W��FCOE���׏d�{��PMGNF�j
      DO IPOIN=1,NPOIN
        DO IDOFN=1,NDOFN
          ITOTV=(IPOIN-1)*NDOFN+IDOFN
          FF(ITOTV)=FCOE(ITOTV)*PMGNF(IDOFN,NFLAG)
        ENDDO
      ENDDO
!
      IF(LODFIX.EQ.-1.AND.NFLAG.GE.2.AND.PPRES(NFLAG).EQ.PPRES(NFLAG-1))THEN
        PRES=PRES
      ELSE
        PRES=PPRES(NFLAG)
      ENDIF
      WIND=WWIND(NFLAG)
      SNOW=SSNOW(NFLAG)
!
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE VOLUME(VTOTAL,COORD,NODEC,NODEM,NPOIN,NELEC,NELEM)
! =====================================================================
! *** ���ʂɓ������Ԃ̑̐�VTOTAL�̌v�Z
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION NODEM(NELEM,4)  ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION NND(3)         ! NND(3) :�O�p�`���\������R�̐ߓ_�ԍ�
      DIMENSION COD(18)        ! COD(18):�O�p�`���������ʏ�ɓ��e���Ăł���U�ʑ̂̐ߓ_���W
!
      VTOTAL=0.0
      ZMIN=0.0
!
      DO IPOIN=1,NPOIN
        IF(COORD(IPOIN,3).LT.ZMIN)ZMIN=COORD(IPOIN,3)
      ENDDO
!
      DO NE=1,NELEM
        NND(1)=NODEM(NE,1)
        NND(2)=NODEM(NE,2)
        NND(3)=NODEM(NE,3)
        DO I=1,3
          COD( I   *3-2)= COORD(NND(I),1)-COORD(NND(1),1)                ! xx(nnd(i)) - xx(nnd(1))
          COD( I   *3-1)= COORD(NND(I),2)-COORD(NND(1),2)                ! yy(nnd(i)) - yy(nnd(1))
          COD( I   *3  )= COORD(NND(I),3)-ZMIN                           ! zz(nnd(i)) - zmin
          COD((I+3)*3-2)= cod(I*3-2)
          COD((I+3)*3-1)= cod(I*3-1)
          COD((I+3)*3  )= 0.0
        ENDDO
!
        VELEM = 0.0
!
        DO NSTEP=1,3
          IF(NSTEP.EQ.1)THEN
            NP1=1
            NP2=2
            NP3=3
          ELSEIF(NSTEP.EQ.2)THEN
            NP1=2
            NP2=5
            NP3=3
          ELSEIF(NSTEP.EQ.3)THEN
            NP1=5
            NP2=6
            NP3=3
          ENDIF
          X12=COD(NP2*3-2)-COD(NP1*3-2)
          Y12=COD(NP2*3-1)-COD(NP1*3-1)
          Z12=COD(NP2*3  )-COD(NP1*3  )
          X13=COD(NP3*3-2)-COD(NP1*3-2)
          Y13=COD(NP3*3-1)-COD(NP1*3-1)
          Z13=COD(NP3*3  )-COD(NP1*3  )
          VOUTX=Y12*Z13-Z12*Y13
          VOUTY=Z12*X13-X12*Z13
          VOUTZ=X12*Y13-Y12*X13
          V0123=(VOUTX*COD(NP1*3-2)+VOUTY*COD(NP1*3-1)+VOUTZ*COD(NP1*3))/6.0
          VELEM=VELEM+V0123
        ENDDO
!
        VTOTAL=VTOTAL+VELEM
!
      ENDDO
!
      RETURN
      END
!
!
!
! =====================================================================
      SUBROUTINE LOAD(CAL,COORD,F,FF,FINIT,NODEC,NVARC,NDOFN,NELEM,NELEC,NFINAL,&
                      NNN,NODEM,NPOIN,NPRINT,PRES,SCF,SNOW,WCF,WIND,NELEF,NODEF,&
                      CALFR,NVARF,MPOND,PROPC,PROPF,PROPM,NVARM,LOADCREM,WWIND,NTIMES)
! =====================================================================
! *** �׏d�x�N�g���̌v�Z
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COORD(NPOIN,3)
      DIMENSION F(NPOIN*NDOFN)
      DIMENSION FINIT(NPOIN*NDOFN)
      DIMENSION FF(NPOIN*NDOFN)
      DIMENSION NODEM(NELEM,4)      ! ���v�f�̃f�[�^�i1�`3�F�ߓ_�ԍ��A4�F���ޔԍ��j
      DIMENSION WCF(NELEM)
      DIMENSION SCF(NELEM)
      DIMENSION NODEC(NELEC,3)   !�P�[�u���v�f�̐ߓ_�ԍ�(1�`2)�A���ޔԍ�(3)
      DIMENSION CAL(NELEC)
      dimension PROPC(NVARC,3)    !�P�[�u���v�f�̍ޗ��f�[�^�i�����O���A�f�ʐρA�P�ʏd�ʁj
      DIMENSION PROPM(NVARM,7,3)        !���v�f�̍ޗ������i��ޔԍ��A��������(1�`7)�A�i�KTri-Linear(1�`3)�j
      DIMENSION PROPF(NVARF,10)           ! �Ȃ��v�f�̍ޗ������i��ޔԍ��A��������(1�`10)�j
      DIMENSION SUMF(6), SUMFF(6), SUMFPRES(6)
      DIMENSION CALFR(NELEF)        ! ���ނ̒���
      DIMENSION NODEF(NELEF,3)      ! �Ȃ��v�fIELEF���\������ߓ_�̔ԍ�(1�`2)�A���ޔԍ�(3)
!
      DIMENSION WWIND(NTIMES)                                             ! WWIND(ILIC)�F����
!
      CALL ZEROR1(FF,NPOIN*NDOFN,NPOIN*NDOFN)
!
      AREAXY=0.0
      AREAXZ=0.0
      AREAYZ=0.0
      QQ=WIND   !WIND�F���x��[N/m2]
      DO NE=1,NELEM
        II=NODEM(NE,1)
        JJ=NODEM(NE,2)
        KK=NODEM(NE,3)
        IVARM=NODEM(NE,4)    ! ���v�f�̕��ޔԍ�
        WEIGHT=PROPM(IVARM,7,1)*PROPM(IVARM,6,1)*9.8 ! ���v�f�̒P�ʏd��[N/m2]=[g/cm3]*10^(-3)*9.8*10^(6)*����[mm]*10^(-3)
        A=(COORD(JJ,2)-COORD(II,2))*(COORD(KK,3)-COORD(II,3))  &          
         -(COORD(JJ,3)-COORD(II,3))*(COORD(KK,2)-COORD(II,2))             ! A=(Y(JJ)-Y(II))*(Z(KK)-Z(II))-(Z(JJ)-Z(II))*(Y(KK)-Y(II))
        B=(COORD(JJ,3)-COORD(II,3))*(COORD(KK,1)-COORD(II,1))  &          
         -(COORD(JJ,1)-COORD(II,1))*(COORD(KK,3)-COORD(II,3))             ! B=(Z(JJ)-Z(II))*(X(KK)-X(II))-(X(JJ)-X(II))*(Z(KK)-Z(II))
        C=(COORD(JJ,1)-COORD(II,1))*(COORD(KK,2)-COORD(II,2))  &          
         -(COORD(JJ,2)-COORD(II,2))*(COORD(KK,1)-COORD(II,1))             ! C=(X(JJ)-X(II))*(Y(KK)-Y(II))-(Y(JJ)-Y(II))*(X(KK)-X(II))
        SMAE=SQRT(A*A+B*B+C*C)*0.5
        SIDEA2=(COORD(JJ,1)-COORD(KK,1))**2+(COORD(JJ,2)-COORD(KK,2))**2+(COORD(JJ,3)-COORD(KK,3))**2
        SIDEB2=(COORD(KK,1)-COORD(II,1))**2+(COORD(KK,2)-COORD(II,2))**2+(COORD(KK,3)-COORD(II,3))**2
        SIDEC2=(COORD(JJ,1)-COORD(II,1))**2+(COORD(JJ,2)-COORD(II,2))**2+(COORD(JJ,3)-COORD(II,3))**2
        SIDEA=SQRT(SIDEA2)
        SIDEB=SQRT(SIDEB2)
        SIDEC=SQRT(SIDEC2)
        CSA=(SIDEB2+SIDEC2-SIDEA2)/(2.0*SIDEB*SIDEC)
        CSB=(SIDEC2+SIDEA2-SIDEB2)/(2.0*SIDEA*SIDEC)
        CSC=(SIDEA2+SIDEB2-SIDEC2)/(2.0*SIDEA*SIDEB)
        IF(CSC.LT.0.0)THEN
          SNA=2.0*SMAE/(SIDEB*SIDEC)
          SNB=2.0*SMAE/(SIDEA*SIDEC)
          SA=SIDEB2*SNA/(8.0*CSA)
          SB=SIDEA2*SNB/(8.0*CSB)
          SC=SMAE-SA-SB
        ELSE
          IF(CSA.LT.0.0)THEN
            SNB=2.0*SMAE/(SIDEA*SIDEC)
            SNC=2.0*SMAE/(SIDEA*SIDEB)
            SB=SIDEC2*SNB/(8.0*CSB)
            SC=SIDEB2*SNC/(8.0*CSC)
            SA=SMAE-SB-SC
          ELSE
            IF(CSB.LT.0.0)THEN
              SNA=2.0*SMAE/(SIDEB*SIDEC)
              SNC=2.0*SMAE/(SIDEA*SIDEB)
              SA=SIDEC2*SNA/(8.0*CSA)
              SC=SIDEA2*SNC/(8.0*CSC)
              SB=SMAE-SA-SC
            ELSE
              R=SIDEA*SIDEB*SIDEC/(4.0*SMAE)
              SUISEN=R*R-0.25*SIDEC2
              IF(SUISEN.LT.0.0) SUISEN=0.0
              SIJ=0.5*SIDEC*DSQRT(SUISEN)
              SUISEN=R*R-0.25*SIDEA2
              IF(SUISEN.LT.0.0) SUISEN=0.0
              SJK=0.5*SIDEA*DSQRT(SUISEN)
              SUISEN=R*R-0.25*SIDEB2
              IF(SUISEN.LT.0.0) SUISEN=0.0
              SKI=0.5*SIDEB*DSQRT(SUISEN)
              SA=0.5*(SIJ+SKI)
              SB=0.5*(SIJ+SJK)
              SC=0.5*(SJK+SKI)
            ENDIF
          ENDIF
        ENDIF
        PA=PRES
        SMAE2=2.0*SMAE
        SI=-QQ*WCF(NE)*SA/SMAE2
        SJ=-QQ*WCF(NE)*SB/SMAE2
        SK=-QQ*WCF(NE)*SC/SMAE2
        SII=PA*SA/SMAE2
        SJJ=PA*SB/SMAE2
        SKK=PA*SC/SMAE2
        F((II-1)*NDOFN+1)=F((II-1)*NDOFN+1)+A*SI
        F((II-1)*NDOFN+2)=F((II-1)*NDOFN+2)+B*SI
        F((II-1)*NDOFN+3)=F((II-1)*NDOFN+3)+C*SI-SNOW*SCF(NE)*SA*C/SMAE2
        F((JJ-1)*NDOFN+1)=F((JJ-1)*NDOFN+1)+A*SJ
        F((JJ-1)*NDOFN+2)=F((JJ-1)*NDOFN+2)+B*SJ
        F((JJ-1)*NDOFN+3)=F((JJ-1)*NDOFN+3)+C*SJ-SNOW*SCF(NE)*SB*C/SMAE2
        F((KK-1)*NDOFN+1)=F((KK-1)*NDOFN+1)+A*SK
        F((KK-1)*NDOFN+2)=F((KK-1)*NDOFN+2)+B*SK
        F((KK-1)*NDOFN+3)=F((KK-1)*NDOFN+3)+C*SK-SNOW*SCF(NE)*SC*C/SMAE2
        FF((II-1)*NDOFN+1)=FF((II-1)*NDOFN+1)+A*SII
        FF((II-1)*NDOFN+2)=FF((II-1)*NDOFN+2)+B*SII
        FF((II-1)*NDOFN+3)=FF((II-1)*NDOFN+3)+C*SII
        FF((JJ-1)*NDOFN+1)=FF((JJ-1)*NDOFN+1)+A*SJJ
        FF((JJ-1)*NDOFN+2)=FF((JJ-1)*NDOFN+2)+B*SJJ
        FF((JJ-1)*NDOFN+3)=FF((JJ-1)*NDOFN+3)+C*SJJ
        FF((KK-1)*NDOFN+1)=FF((KK-1)*NDOFN+1)+A*SKK
        FF((KK-1)*NDOFN+2)=FF((KK-1)*NDOFN+2)+B*SKK
        FF((KK-1)*NDOFN+3)=FF((KK-1)*NDOFN+3)+C*SKK
        AREAXZ=AREAXZ+0.5*A
        AREAYZ=AREAYZ+0.5*B
        AREAXY=AREAXY+0.5*C
!
        IF(NNN.EQ.1)THEN
          FINIT((II-1)*NDOFN+3)=FINIT((II-1)*NDOFN+3)-WEIGHT*SA
          FINIT((JJ-1)*NDOFN+3)=FINIT((JJ-1)*NDOFN+3)-WEIGHT*SB
          FINIT((KK-1)*NDOFN+3)=FINIT((KK-1)*NDOFN+3)-WEIGHT*SC
        ENDIF
!
      ENDDO
!
      IF(NNN.EQ.1)THEN
        DO NE=1,NELEC
          I=NODEC(NE,1)
          J=NODEC(NE,2)
          IV=NODEC(NE,3)
          WAL=CAL(NE)*PROPC(IV,3)*9.8*0.5    !�d��[N]���P�ʏd��[kg/m]*����[m]*9.8
          FINIT((I-1)*NDOFN+3)=FINIT((I-1)*NDOFN+3)-WAL
          FINIT((J-1)*NDOFN+3)=FINIT((J-1)*NDOFN+3)-WAL
        ENDDO
      ENDIF
!
      IF(NELEF.GT.0.AND.NNN.EQ.1)THEN
        DO IE=1,NELEF
          I=NODEF(IE,1)
          J=NODEF(IE,2)
          IV=NODEF(IE,3)
          WEIGHT=PROPF(IV,7)*CALFR(IE)*9.8*0.5     !�d��[N]���P�ʏd��[kg/m]*����[m]*9.8
          FINIT((I-1)*NDOFN+3)=FINIT((I-1)*NDOFN+3)-WEIGHT
          FINIT((J-1)*NDOFN+3)=FINIT((J-1)*NDOFN+3)-WEIGHT
        ENDDO
      ENDIF
!
! *** �׏d�y�ъO�͂̓����ߓ_�͂��t�@�C���ɏo�͂���
!
      IF(NFINAL.NE.1.AND.NPRINT.LE.2.AND.(NNN.EQ.1.OR.NPRINT.EQ.0))THEN
        WRITE(8,105)
  105   FORMAT(//,' (INTERNAL AND EXTERNAL FORCES) [N]',/,2(2X,'NOD',9X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD',5X))
        WRITE(8,115)(K,(F((K-1)*NDOFN+IDOFN)+FF((K-1)*NDOFN+IDOFN)+FINIT((K-1)*NDOFN+IDOFN),IDOFN=1,3),K=1,NPOIN)
  115   FORMAT(2(I5,3F15.4,5X))
!
        CALL ZEROR1(SUMF,6,6)
        CALL ZEROR1(SUMFF,6,6)
        CALL ZEROR1(SUMFPRES,6,6)
        WRITE(8,330)
  330   FORMAT(//,' (DEAD LOAD) [N]',44X,'(INTERNAL PRESSURE LOAD) [N]',/,2(2X,'NOD',9X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD',5X))
        DO I=1,NPOIN
          WRITE(8,320) I,(FINIT((I-1)*NDOFN+IDOFN),IDOFN=1,3),I,(FF((I-1)*NDOFN+IDOFN),IDOFN=1,3)
          DO IDOFN=1,NDOFN
            SUMF(IDOFN)=SUMF(IDOFN)+ F((I-1)*NDOFN+IDOFN)+ FF((I-1)*NDOFN+IDOFN)+ FINIT((I-1)*NDOFN+IDOFN)
            SUMFF(IDOFN)=SUMFF(IDOFN)+FINIT((I-1)*NDOFN+IDOFN)
            SUMFPRES(IDOFN)=SUMFPRES(IDOFN)+FF((I-1)*NDOFN+IDOFN)
          ENDDO
        ENDDO
  320   FORMAT(2(I5,3F15.4,5X))
!
        WRITE(8,325)
        WRITE(8,326) (SUMF(IDOFN),IDOFN=1,NDOFN)
        WRITE(8,327) (SUMFF(IDOFN),IDOFN=1,NDOFN)
        WRITE(8,328) (SUMFPRES(IDOFN),IDOFN=1,NDOFN)
        IF(NFINAL.EQ.1) WRITE(8,130) AREAXY,AREAYZ,AREAXZ
  130   FORMAT(//,' (WHOLE AREA OF THE STRUCTURE PROJECTED ON EACH PLANE)','-M*M-',/' XY-PLANE=',F10.4,5X'YZ-PLANE=',F10.4,5X' XZ-PLANE=',F10.4)
  325   FORMAT(/,105('-'),/,24X,'X-LOAD',9X,'Y-LOAD',9X,'Z-LOAD',8X,'XX-LOAD',8X,'YY-LOAD',8X,'ZZ-LOAD')
  326   FORMAT(3X,' (NODAL LOAD)',3X,6(F11.4,4X),14X,3(F11.4,4X))
  327   FORMAT(3X,' (DEAD  LOAD)',3X,6(F11.4,4X),14X,3(F11.4,4X))
  328   FORMAT(3X,' (PRESS LOAD)',3X,6(F11.4,4X),14X,3(F11.4,4X))
      ENDIF
      RETURN
      END
!
!
!
