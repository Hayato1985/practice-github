! ==================================================================================================
! プログラムの構成
! ARCLM.F90  ：弧長増分法(ARCH LENGTH METHOD)の適用
!   ARCLM 弧長増分法(ARCH LENGTH METHOD)の適用
!   SCPDT 内積の計算
! ======================================================================
      SUBROUTINE ARCLM(INDEX,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,&
                       IFFIX,IITER,IINCS,MTOTV,NDOFN,NOFIX,NORDR,&
                       NTOTV,NVFIX,PRESC,RAMDA,RLOAD,RPRVS,RSCLE,TDISP,&
                       TFACT,TLOAD)
! ======================================================================
! *** 弧長増分法(ARCH LENGTH METHOD)の適用
! INDEX=1：連立一次方程式の解法(SKYLNE)の前の処理
!      =2：連立一次方程式の解法(SKYLNE)の後の処理
!       TFACT       ：λ[k,j]、比例荷重係数
!       RAMDA       ：λ[k,(j)]、第IINCS(=k)回増分による比例荷重係数の変化量
!       FACTO       ：=DRAMD、1回の反復計算による比例荷重係数の増分
!       DSCLE       ：Ｕs、変位のスケーリング係数
!       RSCLE       ：λs、λのスケーリング係数
!  *L   DRAMD       ：Δλ[k,j]、第IITER(=j)回反復による比例荷重係数の変化量
!
      PARAMETER(MCHK=0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ADIS1(MTOTV)				! ADIS1(ITOTV)：ΔＵ'[k,j1]、TLOADに対するスケーリングされた変位量
      DIMENSION ADIS2(MTOTV)				! ADIS2(ITOTV)：ΔＵ'[k,j2]、ELOADに対するスケーリングされた変位量
      DIMENSION ASDIS(MTOTV)				! ASDIS(ITOTV)：Ｕ[k,(j)]、第IINCS回増分による変位量
      DIMENSION ELOAD(MTOTV)				! ELOAD(ITOTV)：残差力ベクトル
      DIMENSION IFFIX(MTOTV)
      DIMENSION NOFIX(NVFIX)
      DIMENSION NORDR(MTOTV)				! NORDR(ITOTV)：元の全体自由度番号ITOTV→並べ替え後の全体自由度番号NORDR(ITOTV)
      DIMENSION PRESC(NVFIX,NDOFN)
      DIMENSION RPRVS(MTOTV+1)			    ! RPRVS(ITOTV+1)：第(IITER-1)回反復後のｒベクトル
      DIMENSION RLOAD(MTOTV)
      DIMENSION TDISP(MTOTV)				! TDISP(ITOTV)：総変位
      DIMENSION TLOAD(MTOTV)				! TLOAD(ITOTV)：荷重モードベクトル
      DOUBLEPRECISION, ALLOCATABLE :: ADIS3(:)
      DOUBLEPRECISION, ALLOCATABLE :: DDISP(:)
      DOUBLEPRECISION, ALLOCATABLE :: RPRSN(:)
      ALLOCATE (ADIS3(MTOTV))				! ADIS3(ITOTV)：Ｕ'[k,(j)]、スケーリングされた変位量ASDIS
      ALLOCATE (DDISP(MTOTV))				! DDISP(ITOTV)：ΔＵ[k,j]、第IITER回反復による変位量→ADIS1に保存される
      ALLOCATE (RPRSN(MTOTV+1))			    ! RPRSN(ITOTV+1)：第(IITER  )回反復後のｒベクトル
      IF(MCHK.EQ.1)WRITE(8,*)' *********** ARCLM ( INDEX=',INDEX,') **********'
!
      IF(IITER.EQ.1.AND.IINCS.EQ.1)THEN
        CALL ZEROR1(ADIS2,MTOTV,NTOTV)
        CALL ZEROR1(RPRVS,MTOTV+1,NTOTV+1)
      ENDIF
      CALL ZEROR1(RPRSN,MTOTV+1,NTOTV+1)
!
! *** INDEX=1の場合、右辺ベクトルADIS1,ADIS2に並べ換え後の順序JTOTVで荷重を代入する
!
    IF(INDEX.EQ.1)THEN
!
      DO ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        ADIS1(JTOTV)=RLOAD(ITOTV)
        ADIS2(JTOTV)=ELOAD(ITOTV)
      ENDDO
!
! *** 右辺ベクトルADIS1,ADIS2の境界(NFREE+1～NTOTVC)に強制変位量を入れておく。
!     SUBROUTINE SKYLINEでの計算により右辺ベクトルADIS1,ADIS2に各自由度ごとの変位が入る。
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
! *** INDEX=2の場合、荷重比例定数増分DRAMDAを計算し、増分変位を求める。
!
    ELSEIF(INDEX.EQ.2)THEN
!
! *** スケーリング係数の設定
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
! *** 変位のスケーリング
!
      DO ITOTV=1,NTOTV
        JTOTV=NORDR(ITOTV)
        ADIS1(JTOTV)=ADIS1(JTOTV)/DSCLE
        ADIS2(JTOTV)=ADIS2(JTOTV)/DSCLE
        ADIS3(JTOTV)=ASDIS(ITOTV)/DSCLE
      ENDDO
!
      RAMDAPRV=RPRVS(NTOTV+1)                       ! RAMDAPRV：前回の反復計算時(IITER-1)のλ
!
! *** 第１回反復（IITER=1）におけるΔＵ[k,1]、Δλ[k,1]の計算
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
          SCPDT3=SCPDT(RPRVS,RPRSN,NTOTV+1,MTOTV+1)   ! RPRVSとRPRSNの内積が正（なす角が鋭角）
          IF(SCPDT3.LT.0.0)THEN                       ! になるようにλの符号を選択する
            DO ITOTV=1,NTOTV+1
              RPRSN(ITOTV)=-RPRSN(ITOTV)
            ENDDO
          ENDIF
        ENDIF
!
! *** 第２回以降反復（IITER≧2）におけるΔＵ[k,j]、Δλ[k,j]の計算
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
! *** RADIAL RETURNによる超球面上への引き戻し
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
        RPRVS(JTOTV)=RPRSN(JTOTV)                       ! 次回の弧長増分のための準備
      ENDDO
      RPRVS(NTOTV+1)=RPRSN(NTOTV+1)                     ! 次回の弧長増分のための準備
!
! *** 比例荷重係数TFACT、総変位TDISP、総荷重TLOADの設定。
!
      TFACT=TFACT+DRAMD
      DO ITOTV=1,NTOTV
! >>>>>>>CMEMの場合          TDISP(ITOTV)=TDISP(ITOTV)+DDISP(ITOTV)
        TLOAD(ITOTV)=RLOAD(ITOTV)*TFACT
        ADIS1(ITOTV)=DDISP(ITOTV)                       ! 今回の反復（第IITER回目）による変位量をADIS1に入れる
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
    ENDIF   ! （ここまでINDEX=2の場合）
!
	DEALLOCATE(ADIS3,DDISP,RPRSN)	! 動的配列を開放する
!
    RETURN
    END
!
!
!
! ======================================================================
      FUNCTION SCPDT(A,B,N,M)
! ======================================================================
! *** 内積の計算
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(M),B(M)
      SCPDT=0.0
      DO 100 I=1,N
        SCPDT=SCPDT+A(I)*B(I)
  100 CONTINUE
      RETURN
      END