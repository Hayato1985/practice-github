   
    
    ! ======================================================================
!     NONLINEAR ANALYSIS OF CABLE REINFORCED MEMBRANE STRUCTURE    
!                                      - EXACT SOL. BY F.E.M.
!                                      - THE NEWTON RAPHSON METHOD
!                                      - SKYLINE MATRIX
!                                      - STRESS TRANSFER METHOD
! ======================================================================
! プログラムの構成
! Cmem7.f90	         メインプログラム、膜およびケーブルの剛性マトリクス作成
! Allowance.f90      鉄骨部材（鋼管）の許容応力度検定
! Arclm.f90          弧長増分法による解析のためのルーチン
! EIGRS.FOR          固有値解析ルーチン
! Frame.f90          曲げ要素の剛性マトリクス作成、曲げ要素関連データのファイル入出力
! IOcntl.f90         INPUTファイルからの入力、OUTPUTファイルへの出力
! loads.f90          自重、外力（風荷重、雪荷重、節点集中荷重）、内圧の等価節点力の計算
! skylne2M.for       スカイライン法による連立一次方程式の解法ルーチン
! TOOLS.F90          補助的なルーチン
! write_vsctipt.f90  要素分割図と応力分布図を VECTOR SCRIPT 形式のファイルで出力するためのサブルーチン
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ 121214 GRANPA DOME  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ IOcntl 風力係数の計算 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ Cmemケーブルのヤング係数PROPC @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
USE WRITE_VSCRIPT_MOD
implicit doubleprecision(A-H,O-Z)
character*15 text
!
!*** ＜動的配列の宣言＞
!
REAL*8 ATMOS_P PARAMETER (ATMOS_P=101325.0D0)
REAL*8 V_INIT
 
doubleprecision, allocatable :: CELL1(:,:)    !Excelの入力データ
doubleprecision, allocatable :: PMGNF(:,:)    !荷重倍率
doubleprecision, allocatable :: PPRES(:)      !内圧[N/m2]
doubleprecision, allocatable :: WWIND(:)      !速度圧[N/m2]
doubleprecision, allocatable :: SSNOW(:)      !雪荷重[N/m2]
doubleprecision, allocatable :: F(:)          !風荷重、雪荷重ベクトル
doubleprecision, allocatable :: FF(:)         !内圧力ベクトル
doubleprecision, allocatable :: FFF(:)        !増分荷重ベクトル
doubleprecision, allocatable :: FINIT(:)      !固定荷重ベクトル
doubleprecision, allocatable :: FCOE(:)       !直接節点荷重の荷重係数
doubleprecision, allocatable :: Q(:)          !残差力ベクトル
doubleprecision, allocatable :: PSUM(:)       !直接節点外力の各方向の合計

doubleprecision, allocatable :: COORD(:,:)    !節点座標
doubleprecision, allocatable :: COOR0(:,:)    !初期形状時の節点座標
doubleprecision, allocatable :: COOR1(:,:)    !ニュートン法による繰返し計算前の 節点座標
doubleprecision, allocatable :: DIS(:)        !変位
doubleprecision, allocatable :: TTDIS(:)      !初期形状に対する変位

integer, allocatable         :: IFFIX(:)      !拘束に関する情報（0:自由、1:拘束)
integer, allocatable         :: NOFIX(:)      !拘束節点の番号
doubleprecision, allocatable :: PRESC(:,:)    !IVFIX番目の拘束節点のIDOFN番目の自由度の拘束量
doubleprecision, allocatable :: PRESC0(:,:)   !IVFIX番目の拘束節点のIDOFN番目の自由度の拘束量（設定値を格納するための配列）
integer, allocatable         :: IFPRE(:)      !拘束状態を表す変数（例　101:x,z方向拘束、y方向自由）

integer, allocatable         :: NODEM(:,:)    !膜要素の節点番号
doubleprecision, allocatable :: PROPM(:,:,:)  !膜要素の材料特性（種類番号、特性項目、段階）
doubleprecision, allocatable :: STRSM(:,:)    !膜要素の張力
doubleprecision, allocatable :: STRM1(:,:)    !リンクリング処理を行う前の膜要素の張力
integer, allocatable         :: NANG(:)       !膜要素IJKのIJ方向と材料の主軸方向が傾いている要素の番号
doubleprecision, allocatable :: ANG(:)        !I膜要素IJKのIJ方向と材料の主軸方向のなす角度
integer, allocatable         :: MEME1(:)      !主応力方向Ｘのリンクリング発生状況（=0：発生している、=1：いない）
integer, allocatable         :: MEME2(:)      !主応力方向Ｙのリンクリング発生状況（=0：発生している、=1：いない）
doubleprecision, allocatable :: AREA(:)       !膜要素の初期形状時の面積
doubleprecision, allocatable :: WCF(:)        !風圧係数
doubleprecision, allocatable :: SCF(:)        !雪荷重係数

integer, allocatable         :: NODEC(:,:)    !ケーブル要素の節点番号(1～2)、部材番号(3)
doubleprecision, allocatable :: BL(:)         !ケーブルの長さ
doubleprecision, allocatable :: CAL(:)        !ケーブルの無ひずみ長さ
doubleprecision, allocatable :: STRSC(:)      !ケーブルの張力
doubleprecision, allocatable :: STRC1(:)      !ゆるみの処理を行う前のケーブルの張力
integer, allocatable         :: NEA(:)        !ケーブルのゆるみの発生状況（=0：発生している、=1：いない）
doubleprecision, allocatable :: PROPC(:,:)    !ケーブル要素の材料データ（ヤング率、断面積、単位重量）

doubleprecision, allocatable :: BETA(:)       !曲げ部材の回転角
doubleprecision, allocatable :: CALFR(:)      !部材の長さ
doubleprecision, allocatable :: VALK(:,:)     !材端1,2の接合条件（＝1.0E+20：剛接、＝0.0：ピン接、その他：半剛接）
integer, allocatable         :: NODEF(:,:)    !曲げ要素を構成する節点の番号(1～2)、部材番号(3)
integer, allocatable         :: NPIN(:)       !部材端部の接合条件（=0：剛剛、=1：ピン剛、=2：剛ピン、=3：ピンピン）
doubleprecision, allocatable :: PROPF(:,:)    !曲げ要素の材料特性（種類番号、特性項目）
doubleprecision, allocatable :: STRSF(:,:)    !曲げ要素IELEFの部材力（基準座標系における部材端荷重)
doubleprecision, allocatable :: STRF1(:,:)    !曲げ要素IELEFの部材力（部材座標系における部材端荷重)

doubleprecision, allocatable :: EE(:,:)       !弾性定数マトリクスＤ
doubleprecision, allocatable :: TSTIF(:)      !全体剛性マトリクス（スカイラインマトリクス）
integer, allocatable         :: NORDR(:)      !全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて並換える時の新しい全体自由度番号
integer, allocatable         :: NWDTH(:)      !全体剛性マトリクスの第ITOTV 行の第１非ゼロ要素から対角要素までの個数
integer, allocatable         :: NWSUM(:)      !全体剛性マトリクスの第ITOTV 行の第１非ゼロ要素から対角要素までの個数

! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>
doubleprecision, allocatable :: ELOAD(:)      !残差力
doubleprecision, allocatable :: TLOAD(:)      !総荷重
doubleprecision, allocatable :: ADIS1(:)      !ΔＵ'[k,j1]、TLOADに対するスケーリングされた変位量
doubleprecision, allocatable :: ADIS2(:)      !ΔＵ'[k,j2]、ELOADに対するスケーリングされた変位量
doubleprecision, allocatable :: ASDIS(:)      !Ｕ[k,(j)]、第IINCS回増分による変位量
doubleprecision, allocatable :: RPRVS(:)      !第(IITER-1)回反復後のｒベクトル
doubleprecision, allocatable :: D(:)          !SUB SKYLNEの作業用配列
doubleprecision, allocatable :: T(:)          !SUB SKYLNEの作業用配列
integer, allocatable         :: NWK(:)        !SUB SKYLNEの作業用配列
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
doubleprecision, allocatable :: CTEMP1(:)     !CTEMP1(LOADCREM)：ケーブルに与える温度増分[°]、MCTEMP=1（ケーブルの温度応力解析を行う）の場合
doubleprecision, allocatable :: CTEMP2(:,:)   !CTEMP2(1,IELEC)：温度応力を与える(=1.0)、与えない(=0.0)、CTEMP2(2,IELEC)：温度ひずみ
!
! *** ＜入出力ファイルの設定＞ ***
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
!*** ＜動的配列の割り当てに関する変数のデータをファイルから読み込む＞
!
!基本的な変数の説明
!NPOIN  節点数
!NFPOIN 自由節点数（NFPOIN+1～NPOINの節点は全拘束）
!NVFIX  拘束節点数
!NELEM  膜要素数
!NVARM  膜要素の種類数
!NPRET  膜要素の初期張力データ、=0：なし、=1：あり
!NELEC  ケーブル要素数
!NVARC  ケーブル要素の種類数
!NPRSTR ケーブルの初期張力データ、=0なし、=1あり
!NELEF  曲げ要素数
!NVARF  曲げ要素の種類数
!NPFRM  曲げ要素の初期応力データ、=0なし、=1あり
!NTIMES 荷重増分法で解析する荷重ケース数
!NCONC  荷重値を直接入力する節点の数

!LODFIX 荷重ベクトルの扱い、=0 変形追従、=1初期値で固定
!MAXCYL ニュートン法の最大繰り返し回数
!DISMAX ニュートン法の収束判定基準
!ERR    節点不釣合力の収束判定基準
!NREPT  不釣合力収束のための最大繰り返し計算回数

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
           ///,'入力データの制御変数')
502 format(9X,'NPOIN ',I10,5X,'節点数')
503 format(9X,'NFPOIN',I10,5X,'自由節点数（NFPOIN+1～NPOINの節点は全拘束）')
504 format(9X,'NVFIX ',I10,5X,'拘束節点数')
505 format(9X,'NELEM ',I10,5X,'膜要素数')
506 format(9X,'NVARM ',I10,5X,'膜要素の種類数')
507 format(9X,'NPRET ',I10,5X,'膜要素の初期張力データ、=0：なし、=1：あり')
508 format(9X,'NELEC ',I10,5X,'ケーブル要素数')
509 format(9X,'NVARC ',I10,5X,'ケーブル要素の種類数')
510 format(9X,'NPRSTR',I10,5X,'ケーブルの初期張力データ、=0なし、=1あり')
511 format(9X,'NELEF ',I10,5X,'曲げ要素数')
512 format(9X,'NVARF ',I10,5X,'曲げ要素の種類数')
513 format(9X,'NPFRM ',I10,5X,'曲げ要素の初期応力データ、=0なし、=1あり')
514 format(9X,'NTIMES',I10,5X,'荷重増分法で解析する荷重ケース数')
515 format(9X,'NCONC ',I10,5X,'直接節点荷重データ、=0なし、=1あり')
516 format(/,'プログラム内部の制御変数')
517 format(9X,'LODFIX',I10,5X,'荷重ベクトルの扱い、=0 変形追従、=1初期値で固定')
518 format(9X,'MAXCYL',I10,5X,'ニュートン法の最大繰り返し回数')
519 format(9X,'DISMAX',F10.4,5X,'ニュートン法の収束判定基準')
520 format(9X,'ERR   ',F10.4,5X,'節点不釣合力の収束判定基準')
521 format(9X,'NREPT ',I10,5X,'不釣合力収束のための最大繰り返し計算回数')
522 format(9X,'NMONT ',I10,5X,'座標値をモニタリングする（WRITE(10,*)で出力する）節点の番号')

nflag=0
if(NTIMES.eq.0)write(8,*)'NTIMESは1以上にしてください'
if(NTIMES.eq.0)nflag=1
if(nflag.eq.1)call echo

!配列サイズに関係するパラメータ
NDOFN = 6                          !節点自由度数
MSTRE = 3                          !膜要素の応力成分の数

!解析オプション
NLOAD  = 0
PRETEN = 0.0                       !NPRET=0のとき、剛性マトリクス特異を避けるための微小な初期張力
DDV    = 1.0                       !膜要素リンクリング時の剛性低減率
DV     = 1.0                       !ケーブルが緩んだ場合に剛性を(1/DV)倍する
NPROB  = NMONT                     !この番号の節点座標データをファイルに出力するWRITE(11,*)

!隠し機能（特別な解析を行うときに設定するパラメータ）
MINCR  = 1                         !荷重増分法（=1）、弧長増分法（=2）
MCTEMP = 0                         !ケーブルの温度応力解析専用の設定
MPOND  = 0                         !MPOND：ポンディング解析を行わない（=0）、行う（=1）
MEIGN  = 0                         !固有値の計算（=0：しない、=1：する）
NONISO = 1                         !膜材料の特性 =0等方性、=1異方性、=2異方性＆繊維方向と膜要素ij方向が傾斜
NANGLE = 0                         !NONISO=2の場合に、該当する膜要素の数
NPRINT = 2                         !出力の制御（=0	計算過程すべて、=1	計算過程の最大変位、最大不釣合力、収束判定、荷重を出力、=2	計算過程の最大変位、最大不釣合力を出力、=3	計算結果のみ）

!できれば削除したい変数
NSTEP  = 1

!
!*** ＜動的配列の割り当て＞

allocate (PMGNF(NDOFN,NTIMES))     !荷重倍率
allocate (PPRES(NTIMES))           !内圧[N/m2]
allocate (WWIND(NTIMES))           !速度圧[N/m2]
allocate (SSNOW(NTIMES))           !雪荷重[N/m2]
allocate (F(NPOIN*NDOFN))          !風荷重、雪荷重ベクトル
allocate (FF(NPOIN*NDOFN))         !内圧力ベクトル
allocate (FFF(NPOIN*NDOFN))        !増分荷重ベクトル
allocate (FINIT(NPOIN*NDOFN))      !固定荷重ベクトル
allocate (FCOE(NPOIN*NDOFN))       !直接節点荷重の荷重係数
allocate (Q(NPOIN*NDOFN))          !残差力ベクトル
allocate (PSUM(NDOFN))             !直接節点外力の各方向の合計

allocate (COORD(NPOIN,3))          !節点座標
allocate (COOR0(NPOIN,3))          !初期形状時の節点座標
allocate (COOR1(NPOIN,3))          !ニュートン法による繰返し計算前の 節点座標
allocate (DIS(NPOIN*NDOFN))        !変位
allocate (TTDIS(NPOIN*NDOFN))      !初期形状に対する変位

allocate (IFFIX(NPOIN*NDOFN))      !拘束に関する情報（0:自由、1:拘束)
allocate (NOFIX(NVFIX))            !拘束節点の番号
allocate (PRESC(NVFIX,NDOFN))      !IVFIX番目の拘束節点のIDOFN番目の自由度の拘束量
allocate (PRESC0(NVFIX,NDOFN))     !IVFIX番目の拘束節点のIDOFN番目の自由度の拘束量（設定値を格納するための配列）
allocate (IFPRE(NVFIX))            !拘束状態を表す変数（例　101:x,z方向拘束、y方向自由）

allocate (NODEM(NELEM,4))          !膜要素のデータ（1～3：節点番号、4：部材番号）
allocate (PROPM(NVARM,7,3))        !膜要素の材料特性（種類番号、特性項目(1～7)、段階Tri-Linear(1～3)）
allocate (STRSM(NELEM,MSTRE))      !膜要素の張力
allocate (STRM1(NELEM,MSTRE))      !リンクリング処理を行う前の膜要素の張力
allocate (NANG(NELEM))             !膜要素IJKのIJ方向と材料の主軸方向が傾いている要素の番号
allocate (ANG(NELEM))              !I膜要素IJKのIJ方向と材料の主軸方向のなす角度
allocate (MEME1(NELEM))            !主応力方向Ｘのリンクリング発生状況（=0：発生している、=1：いない）
allocate (MEME2(NELEM))            !主応力方向Ｙのリンクリング発生状況（=0：発生している、=1：いない）
allocate (AREA(NELEM))             !膜要素の初期形状時の面積
allocate (WCF(NELEM))              !風圧係数
allocate (SCF(NELEM))              !雪荷重係数

allocate (NODEC(NELEC,3))          !ケーブル要素の節点番号(1～2)、部材番号(3)
allocate (BL(NELEC))               !ケーブルの長さ
allocate (CAL(NELEC))              !ケーブルの無ひずみ長さ
allocate (STRSC(NELEC))            !ケーブルの張力
allocate (STRC1(NELEC))            !ゆるみの処理を行う前のケーブルの張力
allocate (NEA(NELEC))              !ケーブルのゆるみの発生状況（=0：発生している、=1：いない）
allocate (PROPC(NVARC,3))          !ケーブル要素の材料データ（ヤング率、断面積、単位重量）

allocate (BETA(NELEF))             !曲げ部材の回転角
allocate (CALFR(NELEF))            !部材の長さ
allocate (VALK(NELEF,2))           !材端1,2の接合条件（＝1.0E+20：剛接、＝0.0：ピン接、その他：半剛接）
allocate (NODEF(NELEF,3))          !曲げ要素IELEFを構成する節点の番号(1～2)、部材番号(3)
allocate (NPIN(NELEF))             !部材端部の接合条件（=0：剛剛、=1：ピン剛、=2：剛ピン、=3：ピンピン）
allocate (PROPF(NVARF,10))         !曲げ要素の材料特性（種類番号、特性項目(1～10)）
allocate (STRSF(NELEF,NDOFN*2))    !曲げ要素IELEFの部材力（基準座標系における部材端荷重)
allocate (STRF1(NELEF,NDOFN*2))    !曲げ要素IELEFの部材力（部材座標系における部材端荷重)

allocate (EE(3,3))                 !弾性定数マトリクスＤ
allocate (NORDR(NPOIN*NDOFN))      !全体剛性マトリクス[K]の全体自由度番号ITOTVをFREEとFIXに分けて並換える時の新しい全体自由度番号
allocate (NWDTH(NPOIN*NDOFN))      !全体剛性マトリクスの第ITOTV 行の第１非ゼロ要素から対角要素までの個数
allocate (NWSUM(0:NPOIN*NDOFN))    !全体剛性マトリクスの第ITOTV 行の第１非ゼロ要素から対角要素までの個数

! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>
allocate (ELOAD(NPOIN*NDOFN))      !残差力
allocate (TLOAD(NPOIN*NDOFN))      !総荷重
allocate (ADIS1(NPOIN*NDOFN))      !ΔＵ'[k,j1]、TLOADに対するスケーリングされた変位量
allocate (ADIS2(NPOIN*NDOFN))      !ΔＵ'[k,j2]、ELOADに対するスケーリングされた変位量
allocate (ASDIS(NPOIN*NDOFN))      !Ｕ[k,(j)]、第IINCS回増分による変位量
allocate (RPRVS(NPOIN*NDOFN+1))    !第(IITER-1)回反復後のｒベクトル
allocate (D(NPOIN*NDOFN))          !SUB SKYLNEの作業用配列
allocate (T(NPOIN*NDOFN))          !SUB SKYLNEの作業用配列
allocate (NWK(NPOIN*NDOFN))        !SUB SKYLNEの作業用配列
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
allocate (CTEMP1(NTIMES))          !CTEMP1(LOADCREM)：ケーブルに与える温度増分[°]、MCTEMP=1（ケーブルの温度応力解析を行う）の場合
allocate (CTEMP2(2,NELEC))         !CTEMP2(1,IELEC)：温度応力を与える(=1.0)、与えない(=0.0)、CTEMP2(2,IELEC)：温度ひずみ

! *** ＜インプットデータファイルから節点座標、膜要素、荷重に関するデータをファイル＃７から読込む＞ ***

CALL INPUT(ANG,COORD,DDV,DISMAX,ERR,FCOE,IFPRE,LODFIX, &
           MAXCYL,MEIGN,MSTRE,NANG,NANGLE,NCONC,NDOFN,NELEM,NFPOIN,&
           NLOAD,NODEM,NOFIX,NONISO,NPOIN,NPRET,NPRINT,NPROB,NREPT,  &
           NSTEP,NTIMES,NVFIX,NVARM,PMGNF,PPRES,PRESC,   &
           PRETEN,PROPM,SCF,SSNOW,STRSM,WCF,WIND,WWIND,MINCR)
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>
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
! *** ＜ケーブル要素に関するデータをファイル＃７から読込む＞ ***
!
IF(NELEC.NE.0)THEN
  CALL CABLE1(BL,CAL,COORD,DV,NODEC,NPOIN,NELEC,NPRSTR,NVARC,&
                  STRSC,STRC1,PROPC,MCTEMP,NTIMES,CTEMP1,CTEMP2)
ENDIF
!
! *** ＜曲げ要素に関するデータをファイル＃７から読込む＞ ***
!
IF(NELEF.NE.0)THEN
  CALL FRAME(1,CALFR,COORD,VALK,NPOIN,NDOFN,NELEF,NODEF,NPIN,&
             NPFRM,NVARF,PROPF,STRSF,STRF1,BETA,IFFIX,&
             MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
ENDIF
!
! *** ＜初期形状のデータをファイル＃９に保存する＞ ***
!
 CALL WRTCOD(1,COORD,NODEC,NELEM,NELEC,NODEM,NPOIN,NELEF,NODEF)
!
! *** ＜要素分割図を VECTOR SCRIPT 形式のファイルで出力する＞
!
LOADCREM=0
CALL WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                   NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,NPROB,NTIMES,&
                   STRSF,STRSM,STRSC,TTDIS)
!
! *** ＜各変数の初期値の設定＞ ***
!
MTOTV=NPOIN*NDOFN
NTOTV=NPOIN*NDOFN
CALL ZEROR1(TTDIS,MTOTV,NTOTV)
CALL ZEROR1(F,MTOTV,NTOTV)
CALL ZEROR1(FF,MTOTV,NTOTV)
NFREE=0
!
! *** ＜CMEM型の拘束節点情報を、SKYLNE型に変更する＞ ***
!
CALL FIXNO(IFFIX,IFPRE,NODEM,NODEC,MTOTV,NDOFN,NELEM,NELEC,NFPOIN,  &
           NFREE,NOFIX,NORDR,NPOIN,NTOTV,NVFIX,NELEF,NODEF)
!
! *** ＜全体剛性マトリクスの成分と１次元配列化されたスカイラインマトリクスの対応を示す指標配列を作成する＞ ***
!
CALL SIZING(NODEM,NODEC,MSIZE,MTOTV,NELEM,NELEC,NDOFN,NFREE,NORDR,     &
            NPOIN,NTOTV,NWDTH,NWSUM,NELEF,NODEF)
!
! *** ＜全体剛性マトリクス（スカイラインマトリクス）の配列割り当て＞ ***
!
allocate (TSTIF(MSIZE))            !MSIZE：スカイラインマトリクスの制限値
!
! ------------------------------------------------------------------------------------------------
! *** ＜アルゴリズムの制御変数＞ ***
!     LOADCREM : 荷重増分回数（１～NTIMES)
!     NNN      : 残差力処理のための繰返し計算回数（1～NREPT+1)
!     INCREM   : 荷重増分回数（１～NSTEP)
!     NCYCLE   : ニュートン･ラフソン法による収束計算の繰返し回数（1～MAXCYL)
!     NCHECK   : =0 収束していない時
!                =1 収束した時
!     NFINAL   : =0 計算途中
!                =1 計算終了直前に,変形後の座標を参照して不釣合力を計算する時。
!     NUNBAL   : =0 ニュートン･ラフソン収束計算中
!                =1 ニュートン。ラフソン収束計算終了後，膜張力を求める時。初期張力０で初回の膜張力を求める時。
!                =2 ニュートン。ラフソン収束計算終了後，膜張力を求める時。
! ------------------------------------------------------------------------------------------------
!
DO 998 LOADCREM=1,NTIMES  ! *** これより（荷重増分1～NTIMES）のループ） ***
  IF(LOADCREM.GT.1)THEN
    NPRET=1
    NSTEP=1
  ENDIF
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>これより
  IF(MINCR.EQ.2)CALL ZEROR1(ASDIS,MTOTV,NTOTV)
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
! *** ＜荷重の設定＞ ***
!
  CALL SETLOAD(LOADCREM,NPOIN,F,PRES,WIND,SNOW,LODFIX,NTIMES,NDOFN,FCOE,PMGNF,PPRES,WWIND,SSNOW)
!
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>これより
! *** ＜温度変化によるケーブルのひずみ増分の設定＞ ***
!
  IF(MCTEMP.EQ.1)THEN   ! MCTEMP=1：ケーブルの温度応力変化を考慮する場合、無視する場合は=0
    DO IELEC=1,NELEC
      IF(CTEMP2(1,IELEC).NE.0.0)THEN
        CTEMP2(2,IELEC)=CTEMP1(LOADCREM)*12.0E-06  ! CTEMP1(LOADCREM)：温度増分[°]、12.0E-06：熱膨張係数[/ °]
      ENDIF
    ENDDO
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
! *** ＜初期内包容積VOL0の計算（内圧変動（ボイル則）を考慮する場合に必要な計算）＞ ***
!
   PRES0=PRES
   APRESOUT=101325.0/9.81   ! 大気圧を10130kg/m2と仮定する
   IF(LODFIX.EQ.-1)THEN
     CALL VOLUME(VOL0,COORD,NODEC,NODEM,NPOIN,NELEC,NELEM)
   ENDIF
!
! *** ＜荷重パラメータの表示> ***
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
! *** ＜ＸＹＺ方向の総荷重の計算＞ ***
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
! *** ＜各変数の初期値の設定＞ ***
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
  DO 3000 NNN=1,NREPT+2  ! *** これより（不釣合力処理のための繰返し計算１～NPREPT+2のループ）***
!   NNN=NREPT+1：収束に至らなかった場合に最終不釣合力の計算を行う
!   NNN=NREPT+2：収束に至らなかった場合にケーブル張力の計算を行う
!
    CALL ZEROR1(F,MTOTV,NTOTV)
!
    IF(LODFIX.NE.2)THEN
!
! >>>>>>>>>> ポンディング解析を行う場合の設定>>>>>>これより
! LOADCREM=2～NTIMES：WWIND(LOADCREM)を貯水率、WINDを水圧[N/m2]、WCF(IELEM)を三角形要素の重心のＺ座標に応じた水圧係数と設定して計算する。
    IF(MPOND.EQ.1)THEN
      IF(LOADCREM.GE.2)THEN
        ZMAX=COORD(1,3)              ! ZMAX：Ｚ座標の最大値
        ZMIN=COORD(1,3)              ! ZMIN：Ｚ座標の最小値
        DO IPOIN=2,NPOIN
          IF(COORD(IPOIN,3).GT.ZMAX)ZMAX=COORD(IPOIN,3)
          IF(COORD(IPOIN,3).LT.ZMIN)ZMIN=COORD(IPOIN,3)
        ENDDO
        ZHGHT=ZMAX-ZMIN             ! ZHGHT：Ｚ座標の最高点から最低点までの高さ
        ZSURF=ZMIN+ZHGHT*WWIND(LOADCREM)  ! ZSURF：水面のＺ座標
        WIND=98000.0          ! WIND：水による圧力98000[N/m2]
        DO IELEM=1,NELEM
          NFLAG=0
          VALZ=0.0            ! VALZ：三角形要素の重心のＺ座標を計算する。
          DO INODE=1,3
            VALZ=VALZ+COORD(NODEM(IELEM,INODE),3)
            IF(COORD(NODEM(IELEM,INODE),3).GT.ZSURF)NFLAG=NFLAG+1
          ENDDO
          WCF(IELEM)=(ZSURF-VALZ/3.0)                     ! WCF(IELEM)：水圧係数（＝水深）
          IF(WCF(IELEM).LT.0.0)WCF(IELEM)=0.0             ! 重心のＺ座標が水面より上にある場合
          IF(NFLAG.EQ.1)WCF(IELEM)=WCF(IELEM)*2.0/3.0     ! ３節点のうち１点が水面より上にある場合
          IF(NFLAG.EQ.2)WCF(IELEM)=WCF(IELEM)*1.0/3.0     ! ３節点のうち２点が水面より上にある場合
          IF(NFLAG.EQ.3)WCF(IELEM)=0.0                    ! ３節点すべてが水面より上にある場合
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
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
! *** ＜変形後の内包容積VOL1の計算（内圧変動（ボイル則）を考慮する場合に必要な計算）＞ ***
!
      IF(NNN.GE.2.AND.LODFIX.EQ.-1.AND. &
        (WIND.NE.0.0.OR.SNOW.NE.0.0.OR.PSUM(1).NE.0.0.OR.PSUM(2).NE.0.0.OR.PSUM(3).NE.0.0))THEN
        CALL VOLUME(VOL1,COORD,NODEC,NODEM,NPOIN,NELEC,NELEM)
        IF(VOL1.EQ.0.0)STOP ' VOL1=0.0 (STOP AT LINE 327) '
!       PRES=PRES0*VOL0/VOL1
        VEXTRA=VOL0*200.0   ! 解析対象の構造物がVOL0の200倍の容積の剛体に連結されていると仮定。
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
! *** ＜荷重ベクトルの計算＞ ***
!
      CALL LOAD(CAL,COORD,F,FF,FINIT,NODEC,NVARC,NDOFN,NELEM,NELEC,NFINAL,&
                NNN,NODEM,NPOIN,NPRINT,PRES,SCF,SNOW,WCF,WIND,NELEF,NODEF,&
                CALFR,NVARF,MPOND,PROPC,PROPF,PROPM,NVARM,LOADCREM,WWIND,NTIMES)
!
      IF(LODFIX.EQ.1)THEN
        LODFIX=2  ! LODFIX=1（荷重ベクトル固定）の場合には、ここでLODFIX=2となり次回以降荷重ベクトルが更新されなくなる
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
    IF(NNN.GE.2) NSTEP=1  ! 残差力収束のための反復回数NNN=2回目以降は、荷重増分を行わない。
!
    DO ITOTV=1,NTOTV
      FFF(ITOTV)=F(ITOTV)+FINIT(ITOTV)
      IF(NLOAD.EQ.0)FFF(ITOTV)=(FFF(ITOTV)+FF(ITOTV))/FLOAT(NSTEP)
      IF(NLOAD.GE.1)FFF(ITOTV)=FFF(ITOTV)/FLOAT(NSTEP)
    ENDDO
!
    DO 2000 INCREM=1,NSTEP  ! *** これより（荷重分割１～NSTEPのループ）***
!
      CALL ZEROR1(DIS,MTOTV,NTOTV)
!
      IF(NFINAL.EQ.0)THEN  ! NFINAL=0：計算途中の場合
        DO IPOIN=1,NPOIN
          DO IDIR=1,3
            COOR1(IPOIN,IDIR)=COORD(IPOIN,IDIR)
          ENDDO
        ENDDO
      ENDIF
!
! *** ＜荷重増分の設定＞ ***
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
! *** ＜ケーブル張力の計算＞ ***
!
      IF(NELEC.NE.0.AND.(NNN.GT.1.OR.INCREM.GT.1))THEN
        CALL CABLE3(0,BL,CAL,NELEC,NVARC,NODEC,STRSC,STRC1,PROPC,MCTEMP,CTEMP2)
      ENDIF
!
      NUNBAL=0
      IF(NFINAL.EQ.1) NUNBAL=2
!
      DO 1000 NCYCLE=1,MAXCYL  ! ***これより（ニュートン法の反復計算のループ）***
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
! **** ＜膜要素の剛性マトリクスを作成する＞ ***
!
        460 CONTINUE
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>これより
        IF(MINCR.EQ.2)THEN
          DO ITOTV=1,NTOTV
            FFF(ITOTV)=TLOAD(ITOTV)
          ENDDO
        ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
        IF(NELEM.GT.0)THEN
          CALL MEMBRN(ANG,AREA,COOR1,COORD,DDV,EE,IFFIX,INCREM,MEME1,&
                      MEME2,NPOIN,MSIZE,MSTRE,MTOTV,NANG,NANGLE,NCYCLE,&
                      NDOFN,NELEM,NFINAL,NNN,NODEM,NONISO,NORDR,NPROB,&
                      NUNBAL,NVARM,NWSUM,PRETEN,PROPM,Q,STRM1,STRSM,TSTIF)
        ENDIF
!
! *** ＜ケーブル要素の剛性マトリクスを作成する＞ ***
!
        IF(NELEC.NE.0)THEN
          CALL CABLE2(BL,CAL,COORD,COOR1,DV,IFFIX,INCREM,NODEC,NVARC,&
                      MSIZE,MTOTV,NCYCLE,NDOFN,NEA,NELEC,NFPOIN,NNN,NORDR,&
                      NPOIN,NPRET,NPRSTR,NTOTV,NUNBAL,NWSUM,Q,&
                      TSTIF,STRC1,PROPC,MCTEMP,CTEMP2)
        ENDIF
!
! *** ＜曲げ要素の剛性マトリクスを作成する＞ ***
!
        IF(NELEF.GT.0)THEN
          CALL FRAME(2,CALFR,COORD,VALK,NPOIN,NDOFN,NELEF,NODEF,&
                      NPIN,NPFRM,NVARF,PROPF,STRSF,STRF1,BETA,&
                      IFFIX,MSIZE,MTOTV,NORDR,NWSUM,TSTIF,DIS,Q)
        ENDIF
!
!
        IF(NUNBAL.NE.0)GOTO 730  ! ***（NUNBAL≠0ならば，ニュートン法の反復計算のループから抜ける）***
!
! *** ＜全体剛性マトリクスTSTIFの固有値を計算する＞ ***
!
        IF(MEIGN.EQ.1)THEN
          CALL EIGEN(EIGMN,MSIZE,MTOTV,NFREE,NORDR,NPOIN,NWDTH,NWSUM,FFF,TSTIF)
        ENDIF
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>これより
        IF(MINCR.EQ.2)THEN            
          IDOFN=3
          ITOTV=(NPROB-1)*NDOFN+IDOFN
          WRITE(8,8002)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,Q(ITOTV),FFF(ITOTV)
 8002     FORMAT(  'LOADCREM      NNN  NCYCLE   NPROB   IDOFN    TFACTT       Q(ITOTV)     FFF(ITOTV)',/,5I8,4G15.4)
        ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
        DO ITOTV=1,NTOTV
          Q(ITOTV)=FFF(ITOTV)-Q(ITOTV)
          F(ITOTV)=Q(ITOTV)
        ENDDO
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>これより
            IF(MINCR.EQ.2)THEN
              DO ITOTV=1,NTOTV
                ELOAD(ITOTV)=Q(ITOTV)
              ENDDO
            ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
        CALL ZEROR1(DIS,MTOTV,NTOTV)
!
! *** ＜収束判定，収束時NCHECK=1にはNFINAL=1として最後の繰返し計算へ進む＞ ***
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
        IF(MINCR.EQ.1)THEN  ! ◆◆◆ 荷重増分法の場合の計算 ◆◆◆
!
! *** ＜係数マトリックスおよび右辺ベクトルに対して、幾何学的な条件を処理する＞ ***
!
          CALL BCONUU(F,IFFIX,MSIZE,MTOTV,NDOFN,NFREE,NOFIX,NORDR,NTOTV,NVFIX,NWSUM,PRESC,TSTIF)
!
! *** ＜スカイライン法により連立一次方程式の解を求める＞ ***
!
          CALL SKYLNE(1,TSTIF,F,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
!
! *** ＜増分変位ベクトルDISの算定＞ ***
!
          CALL RARNG(DIS,F,MTOTV,NORDR,NTOTV)
!
        ELSEIF(MINCR.EQ.2)THEN ! ◆◆◆ 弧長増分法の場合の計算 ◆◆◆  
!
! *** 弧長増分法(ARCH LENGTH METHOD)の適用(INDEX=1：連立一次方程式の解法(SKYLNE)の前の処理)
!
          CALL ARCLM(1,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,IFFIX,NNN,LOADCREM,&
                     MTOTV,NDOFN,NOFIX,NORDR,NTOTV,NVFIX,PRESC,RAMDA,FCOE,RPRVS,RSCLE,&
                     TTDIS,TFACT,TLOAD)
!
! *** 全体剛性マトリクスTSTIFの固有値を計算する
!
          CALL EIGEN(EIGN1,MSIZE,MTOTV,NFREE,NORDR,NPOIN,NWDTH,NWSUM,ADIS1,TSTIF)
!
! *** スカイライン法により連立一次方程式の解を求める
!
          CALL SKYLNE(1,TSTIF,ADIS1,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
          IF(NNN.GT.1)THEN
            CALL SKYLNE(2,TSTIF,ADIS2,NWSUM,MSIZE,NFREE,MTOTV,IER,D,T,NWK)
          ENDIF
!
! *** 弧長増分法(ARCH LENGTH METHOD)の適用(INDEX=2：連立一次方程式の解法(SKYLNE)の後の処理)
!
          CALL ARCLM(2,ADIS1,ADIS2,ARCLG,ASDIS,DSCLE,ELOAD,DRAMD,IFFIX,&
                     NNN,LOADCREM,MTOTV,NDOFN,NOFIX,NORDR,NTOTV,NVFIX,&
                     PRESC,RAMDA,FCOE,RPRVS,RSCLE,TTDIS,TFACT,TLOAD)
!
          DO ITOTV=1,NTOTV
            DIS(ITOTV)=ADIS1(ITOTV)
          ENDDO
!
        ENDIF  ! ◆◆◆ ここまで（MINCRによる条件分岐）◆◆◆ 
!
!
! *** ＜節点座標の更新＞ ***
!
        DO IPOIN=1,NFPOIN
          DO IDIR=1,3
            ITOTV=(IPOIN-1)*NDOFN+IDIR
            COORD(IPOIN,IDIR)=COORD(IPOIN,IDIR)+DIS(ITOTV)
          ENDDO
        ENDDO
!
! *** ＜初期形状に対する変位ベクトルTTDISの算定＞ ***
!
        DO ITOTV=1,NTOTV
          TTDIS(ITOTV)=TTDIS(ITOTV)+DIS(ITOTV)
        ENDDO
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>これより
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
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
! *** ＜NCYCLE=INCREM=NNN=1の時，強制変位量をゼロにする＞ ***
!
        IF(NCYCLE.EQ.1.AND.NNN.EQ.1)THEN
          CALL ZEROR2(PRESC,NVFIX,NDOFN,NVFIX,NDOFN)
        ENDIF
!
! *** ＜計算過程をファイル＃８に保存する＞ ***
!
        CALL PROC(0,DIS,EIGMN,IFFIX,INCREM,LOADCREM,NPOIN,NCYCLE,NDOFN,NFPOIN,&
                  NNN,NPRINT,NPROB,NTIMES,NTOTV,PRES,Q,COORD,QMAX,MINCR)
!
! *** ＜増分変位量による収束判定＞ ***
!
        NDCHK=0
        DO ITOTV=1,NTOTV
          IF(IFFIX(ITOTV).EQ.0.AND.ABS(DIS(ITOTV)).GT.DISMAX)THEN
            NDCHK=1
            GOTO 1000
          ENDIF
        ENDDO
!
! *** ＜NDCHK=0(収束)またはNCYCLE≧MAXCYLの場合，NUNBAL=1or2として行番号460のDOループへ戻る＞ ***
!
! *** ＜計算過程の出力＞
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
      1000 CONTINUE ! *** ここまで（ニュートン法の反復計算のループ）***
!
      730 IF(NFINAL.EQ.1) GO TO 3100  ! ***（NFINAL=1ならば，計算を終了し結果を出力する）***
!
! *** ＜残差力Ｑの計算＞
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
    2000 CONTINUE ! *** ここまで（荷重分割のループ）***
!
    ! ニュートン法の反復計算で収束しなかった場合には、LASTCYL=NCYCLE,NFINAL=1とする。 
    IF(NFINAL.NE.1.AND.NNN.EQ.NREPT+1)THEN
      LASTCYL=NCYCLE
      NFINAL=1
    ENDIF
!
  3000 CONTINUE ! *** ここまで（残差力処理のための繰返し計算のループ）***
!
  3100 CONTINUE
!
  DO NF = 9, 10
    WRITE(NF,178)LOADCREM,NTIMES,(PMGNF(IDOFN,LOADCREM),IDOFN=1,NDOFN),PRES,WIND,SNOW
  ENDDO
!
! >>>>>>>>>> 弧長増分法専用の設定>>>>>>>>>>>>>>>>>>ここから
  IF(MINCR.EQ.2)THEN
    IDOFN=3
    ITOTV=(NPROB-1)*NDOFN+IDOFN
    IF(LOADCREM.EQ.1)WRITE(13,8000)
    WRITE(13,8000)LOADCREM,NNN,NCYCLE,NPROB,IDOFN,TFACT,DRAMD,TLOAD(ITOTV),ELOAD(ITOTV),DIS(ITOTV),TTDIS(ITOTV),ASDIS(ITOTV),Qmax,EIGMN
  ENDIF
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<ここまで
!
  IF(LODFIX.EQ.-1)WRITE(8,8100) PRES0,PRES,PRES-PRES0,VOL0,VOL1,VOL1-VOL0
  8100 FORMAT(//,15X,'     INITIAL     PRESENT     DIFFER.',/,5X,' PRESSURE:',3F12.3,' [N/m2]',/,5X,' VOLUME:  ',3F12.3,' [m*m*m] ')
!
! *** ＜計算結果の出力＞
!
  CALL RESLT(BL,CAL,COOR0,COORD,FF,FFF,IFFIX,MSTRE,NDOFN,NELEM,&
             NELEC,NFINAL,NFPOIN,NPOIN,NVARC,NVARM,NODEC,NODEM,Q,STRSM,STRSC,TTDIS,&
             STRC1,PROPC,PROPM,NELEF,STRSF,MCTEMP,CTEMP2,NPROB,NVARF,NODEF)
!
! *** ＜変形後の座標で荷重を計算＞
!
      NNN=1
  CALL LOAD(CAL,COORD,F,FF,FINIT,NODEC,NVARC,NDOFN,NELEM,NELEC,NFINAL,&
            NNN,NODEM,NPOIN,NPRINT,PRES,SCF,SNOW,WCF,&
            WIND,NELEF,NODEF,CALFR,NVARF,MPOND,PROPC,PROPF,PROPM,NVARM,LOADCREM,WWIND,NTIMES)
!
! *** ＜要素分割図と応力分布図を VECTOR SCRIPT 形式のファイルで出力する＞
!
  IF(MINCR.NE.2.AND.LASTCYL.GE.2)NTIMES=LOADCREM    ! ニュートン法の反復計算で収束しなかった場合（LASTCYL>1)
! IF(LOADCREM.EQ.1.OR.LOADCREM.EQ.NTIMES)THEN
    CALL WRITE_VSCRIPT(COORD,COOR0,NODEC,LOADCREM,MSTRE,NVARF,NELEC,&
                       NDOFN,NELEF,NELEM,NODEM,NODEF,NPOIN,ABS(NPROB),NTIMES,&
                       STRSF,STRSM,STRSC,TTDIS)
    CALL CAL_PRES(NTIMES, NPOIN, NELEM, COORD, NODEM, PPRES, V_INIT, ATMOS_P, Z_REF)
! ENDIF
!
! *** ＜曲げ部材の検定＞
!
  CALL ALLOWANCE(COORD,NPOIN,NVARF,NDOFN,NELEF,NODEF,PROPF,STRSF,VALK)
!
  ! ニュートン法の反復計算で収束しなかった場合（LASTCYL>1)の場合は計算を終了する。
  IF(MINCR.NE.2.AND.LASTCYL.GE.2)THEN
    WRITE(8,*)' LASTCYL > 1'
    STOP ' LASTCYL > 1'
  ENDIF
!
998 CONTINUE ! *** ここまで（荷重分割NTIMESのループ）

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
! 材料の主軸の方向に合わせて弾性定数マトリクスを変換する
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
! ケーブル要素の剛性マトリクスを作成する
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NORDR(MTOTV)
      DIMENSION NWSUM(0:MTOTV)
      DIMENSION TSTIF(MSIZE)
      DIMENSION IFFIX(MTOTV)
      DIMENSION COORD(NPOIN,3)
      DIMENSION COOR1(NPOIN,3)
      DIMENSION NODEC(NELEC,3)   !ケーブル要素の節点番号(1～2)、部材番号(3)
      DIMENSION BL(NELEC)
      DIMENSION CAL(NELEC)
      DIMENSION STRC1(NELEC)
      DIMENSION NEA(NELEC)
      dimension PROPC(NVARC,3)    !ケーブル要素の材料データ（ヤング率、断面積、単位重量）
      DIMENSION Q(MTOTV)
      DIMENSION PP(MTOTV)
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
!     MCTEMP                        ! MCTEMP：ケーブルの温度応力解析を行わない（=0）、行う（=1）
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)：温度応力を与える(=1.0)、与えない(=0.0)、CTEMP2(2,IELEC)：温度ひずみ
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      DIMENSION SM(6,6),IND(6),TT(6,6),U(6),DISE(6),BMT(6),GM(6,6)
      DIMENSION BB(6,6),BBG(6,6),BEBT(6,6),P(6)
!
      CALL ZEROR1(PP,MTOTV,NTOTV)  ! ケーブルまたはトラスによる等価節点力ベクトル
!
      DO 210 NE=1,NELEC
        IV=NODEC(NE,3)
        EA=PROPC(IV,1)*ABS(PROPC(IV,2))  !EA[N]=ヤング率[N/mm2]*断面積[mm2]
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
!@@@@@@ GRANPA DOME のヤング率（トリリニアモデル）の計算 @@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        EEA=EA
        IF(PROPC(IV,2).GT.0.0.AND.NEA(NE).EQ.0)THEN     !PROPC(IV,2)断面積
          EEA=EEA/DV
        ENDIF
!
        IF(MCTEMP.EQ.0)THEN     ! MCTEMP=0：ケーブルの温度応力解析を行わない場合
          SSION=(BL(NE)-CAL(NE))*EA/CAL(NE)
          SION =BL(NE)*STRC1(NE)/CL+(BL(NE)-CL)*EEA/CL
        ELSE                    ! MCTEMP=1：ケーブルの温度応力解析を行う場合                    
          CAL1=CAL(NE)*(1.0+CTEMP2(2,NE))
          SSION=((BL(NE)-CAL1)*EA/CAL1)
          SION =((BL(NE)-CAL1)*EA/CAL1)
        ENDIF
!
        IF(NUNBAL.NE.0)THEN
          NEA(NE)=1
          IF(SSION.LT.0.0.AND.PROPC((IV),2).GT.0.0)THEN     !PROPC(IV,2)断面積
            SION=0.0
            NEA(NE)=0
          ENDIF
        ENDIF
!
! *** ＜NUNBAL=0：ニュートン法による収束計算中＞
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
! *** ＜全体剛性マトリクスへの重ね合わせ＞
!
        IF(NUNBAL.EQ.0)THEN
          CALL MATRIX(SM,NE,NODEC,NELEC,MSIZE,NDOFN,3,2,NORDR,NWSUM,MTOTV,TSTIF,IFFIX)
        ENDIF
!
! *** ＜ケーブル張力の等価節点力ベクトルＰＰの計算＞ 
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
! *** ＜ケーブル張力の等価節点力ベクトルＰＰを全体の節点力ベクトルＱに加える＞ 
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
! *** 残差力の収束判定基準値ERRに対する判定（収束していればNCHECK=1）
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
! ケーブル要素の緩みと張力を計算し、結果をファイルに出力する
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BL(NELEC)
      DIMENSION CAL(NELEC)
      DIMENSION STRSC(NELEC)
      DIMENSION STRC1(NELEC)
      DIMENSION NODEC(NELEC,3)       !ケーブル要素の節点番号(1～2)、部材番号(3)
      dimension PROPC(NVARC,3)    !ケーブル要素の材料データ（ヤング率、断面積、単位重量）
      DIMENSION VALMM(NVARC,2)      ! VALMM(ICVR, )：部材ICVRの最大応力、最小応力
      DIMENSION NVALMM(NVARC,2)     ! NVALMM(ICVR, )：部材ICVRの最大応力、最小応力に対応する要素番号
!      DIMENSION NFAIL(5000)         ! NFAIL：ゆるみを生じているケーブルの要素番号
integer, allocatable :: NFAIL(:)    ! ゆるみを生じているケーブルの要素番号
! >>>>>>>>>> ケーブルの温度応力解析専用の設定>>>>>>
!     MCTEMP                        ! MCTEMP：ケーブルの温度応力解析を行わない（=0）、行う（=1）
      DIMENSION CTEMP2(2,NELEC)      ! CTEMP2(1,IELEC)：温度応力を与える(=1.0)、与えない(=0.0)、CTEMP2(2,IELEC)：温度ひずみ
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IFAIL=0   ! IFAIL：ゆるみが生じているケーブルの本数
!
! *** ＜ケーブル部材（断面積の符号がプラスの部材）に圧縮力が生じている場合には、ゼロにする＞
      DO NE=1,NELEC
        IV=NODEC(NE,3)
        EA=PROPC(IV,1)*ABS(PROPC(IV,2))   !EA[N]=ヤング率[N/mm2]*断面積[mm2]
        IF(MCTEMP.EQ.0)THEN     ! MCTEMP=0：ケーブルの温度応力解析を行わない場合
          STRC1(NE)=(BL(NE)-CAL(NE))*EA/CAL(NE)
        ELSE                    ! MCTEMP=1：ケーブルの温度応力解析を行う場合
          STRC1(NE)=((BL(NE)-CAL(NE))*EA/CAL(NE))-CTEMP2(2,NE)*EA
        ENDIF
        STRSC(NE)=STRC1(NE)
        IF((STRSC(NE).LT.0.0).AND.(PROPC(IV,2).GT.0.0))THEN        !PROPC(IV,2)断面積
          STRSC(NE)=0.0
          IFAIL=IFAIL+1
        ENDIF
      ENDDO
!
allocate (NFAIL(IFAIL))           ! ゆるみを生じているケーブルの要素番号
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
! *** ＜ケーブル等の応力の最大値VALMM(ISTRE,1)と最小値VALMM(ISTRE,2)をファイル＃８に保存する（INDEX=1の場合）＞
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
! *** ＜ゆるみが生じたケーブルの個数IFAILとその部材番号NFAIL(IFAIL)をファイルに書き込む＞
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
! **** ＜膜要素の剛性マトリクスを作成する＞ ***
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
      DIMENSION NODEM(NELEM,4)          ! 膜要素のデータ（1～3：節点番号、4：部材番号）
      DIMENSION NORDR(NPOIN*NDOFN)
      DIMENSION NWSUM(0:NPOIN*NDOFN)
      DIMENSION PROPM(NVARM,7,3)        !膜要素の材料特性（種類番号、特性項目(1～7)、段階Tri-Linear(1～3)）
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
IVARM=NODEM(NE,4)   !NODEM(NE,4)：部材番号
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
! *** ＜Tri-Linearモデルによる材料剛性の設定＞ ***
!
  CNV=1.0E+00               ! CNV：単位の換算係数（長さがmm単位の場合1.0、長さがm単位の場合1000.0）
  TH=PROPM(IVARM,6,1)       ! PROPM(i,6,1)：膜厚[mm]
  F=SQRT(PX**2+PY**2-PX*PY+3*PXY**2)    ! F：相当応力[N/m]
  F=F/CNV/TH             ! F：単位の変換（座標がmm単位の場合N/mm2、座標がm単位の場合N/m→N/mm2）
  if(F.lt.PROPM(IVARM,5,1))then              ! PROPM(IVRAM,5,1)：第１降伏点応力σy
    j=1
  elseif(F.lt.PROPM(IVARM,5,2))then          ! PROPM(IVRAM,5,2)：第２降伏点応力σy
    j=2
  else
    j=3                                 ! j=1,2,3（第１～第３勾配）
  endif
  ET1=PROPM(IVARM,1,j)*TH*CNV     !引張剛性ET=ヤング率E*厚さt
  ET2=PROPM(IVARM,2,j)*TH*CNV
  POA1=PROPM(IVARM,3,j)
  POA2=POA1*ET2/ET1           ! 相反の定理を満足するようにPOA2を設定する
  GXY=PROPM(IVARM,4,j)*TH*CNV

!  ET1=PROPM(IVARM,1,1)
!  ET2=PROPM(IVARM,2,1)
!  POA1=PROPM(IVARM,3,1)
!  POA2=POA1*ET2/ET1           ! 相反の定理を満足するようにPOA2を設定する
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
! 注意 EEをMAIN PROGRAM側に移動したプログラムがあるが、材料非線形で解析を行う場合には、ここで更新する必要がある。
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
! *** ＜膜要素IJKのIJ方向と材料の主軸方向が傾いている場合の弾性剛性マトリクスＤの変換＞
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
! *** ＜要素剛性マトリクスの作成＞ ***
!
        CALL NEWTON(SM,PX,PY,PXY,P,XXI,YYI,ZZI,XXJ,YYJ,ZZJ,XXK,YYK,ZZK,U,NUNBAL,E,PPX,PPY,PPXY,EE,NFINAL,NE,SURF)
!
! *** ＜リンクリングの判定と処理＞ ***
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
! *** ＜膜張力の等価節点力を全体の節点力ベクトルＱに加える＞ ***
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
! *** ＜要素剛性マトリックス→全体剛性マトリックス＞ ***
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
