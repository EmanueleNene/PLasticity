c --- definitions for buffer status
      integer   W_ReadAndWrite,    W_Read,     W_Write
      parameter(W_ReadAndWrite=0,  W_Read=1,   W_Write=2)

c --- definitions of integer values of characters
      integer   W_ChBlank,          W_ChAster,       W_ChSlash,
     &          W_ChZero,		    W_ChBackSlash,
     &          W_ChPeriod,         W_ChTilda
      character*32 W_StringBlank
      parameter(W_ChBlank=32,       W_ChAster=42,    W_ChSlash=47,
     &          W_ChPeriod=46,      W_ChTilda=126,
     &          W_ChZero=48,        W_ChBackSlash=92,
     &          W_StringBlank='                                ')

c --- useful constant :)
      double precision  W_Pi
      parameter  (W_Pi=3.14159265358979323846d0)

c --- maximum number of nodes
      integer    W_NdMax
      parameter (W_NdMax=28)

c --- maximum row size
      integer    W_NrMax
      parameter (W_NrMax=W_NdMax*6)

      integer          W_FirstMatIndx, W_LastWrkIndx, W_FirstWrkIndx
      parameter       (W_FirstWrkIndx = 1,
     &                 W_LastWrkIndx  = 5*6 + 36,
     &                 W_FirstMatIndx = W_FirstWrkIndx+W_LastWrkIndx)
c
      integer    W_LenWrkIndx,
     &           W_WrkIndx1,  W_WrkIndx2,  W_WrkIndx3,
     &           W_WrkIndx4,  W_WrkIndx5,  W_WrkIndx6,
     &           W_WrkIndx7,  W_WrkIndx8,  W_WrkIndx9
      parameter (W_LenWrkIndx = 9,
     &           W_WrkIndx1  =  1,
     &           W_WrkIndx2  =  2,
     &           W_WrkIndx3  =  3,
     &           W_WrkIndx4  =  4,
     &           W_WrkIndx5  =  5,
     &           W_WrkIndx6  =  6,
     &           W_WrkIndx7  =  7,
     &           W_WrkIndx8  =  8,
     &           W_WrkIndx9  =  9)

      integer    W_HedLen
      parameter (W_HedLen = 103)
c  W_HedLen  starting location after file header read/write
c            W_HedLen =   1    +    1    +   100   +   1
c                        byte      kbfint   record    byte
c                       count               length   count


      integer    W_LenMatRecIndx, W_LenMatWrkIndx,
     &           W_iStrain, W_iStress, W_iStrainHyg, W_iStrainPl,
     &           W_iEqvPl, W_iStrainCr, W_iEqvCr, W_iStrainIni,
     &           W_iMisc, W_iMatMatrix, W_iStrainInc, W_iStateVar
c
      parameter (W_LenMatRecIndx = 10,
     &           W_LenMatWrkIndx = 12,
     &           W_iStrain      =  1,
     &           W_iStress      =  2,
     &           W_iStrainHyg   =  3,
     &           W_iStrainPl    =  4,
     &           W_iEqvPl       =  5,
     &           W_iStrainCr    =  6,
     &           W_iEqvCr       =  7,
     &           W_iStrainIni   =  8,
     &           W_iMisc        =  9,
     &           W_iStateVar    = 10,
     &           W_iMatMatrix   = 11, ! not saved
     &           W_iStrainInc   = 12) ! not saved

c	W_iStrain     -    elastic strain
c	W_iStress     -    stress
c	W_iStrainHyg  -    hygro strtain
c	W_iStrainPl   -    plastic strain
c	W_iEqvPl      -    equivalent plastic strain and plastic work
c	W_iStrainCr   -    creep strain
c	W_iEqvCr      -    equivalent creep strain and creep work
c	W_iStrainIni  -    initial strtain
c	W_iMisc       -    miscellaneous records
c	W_iMatMatrix  -    elastic matrix
c	W_iStrainInc  -    incremental strain


    
c --- define the element work index for the common features of elements
c     it is only stabilization variables now
      integer    W_LenStabIndx,
     &           W_iStabMat,     W_iStabForc,     W_iStabVar
      parameter (W_LenStabIndx = 3,
     &           W_iStabMat = 1, W_iStabForc = 2, W_iStabVar = 3)
     
c --- all system dependent data lengths
      integer  W_AddrStrideDp,W_AddrStrideI,W_IntPerDp,W_IntPerLong

	  parameter (W_AddrStrideDp=8)

#if defined(FULL_64)
	  parameter (W_IntPerDp=1,W_IntPerLong=1,W_AddrStrideI=8)
#else
      parameter (W_IntPerDp=2,W_IntPerLong=2,W_AddrStrideI=4)
#endif

c	W_AddrStrideI - address stride in integer vector
c	W_AddrStrideDp - address stride in double precision vector

c --- definitions for calls to memmory allocator
      integer   W_LenByte, W_LenInt, W_LenDouble, W_LenLongInt

#if defined(ANSYS_MEM)
      parameter(W_LenInt=1,    W_LenByte=0,
     &          W_LenDouble=3, W_LenLongInt=6)
#else
      parameter(W_LenInt=8/W_IntPerLong, W_LenByte=1,
     &          W_LenDouble=8, W_LenLongInt=8)
#endif
c --- number of contact params
      integer    W_NumContParam
      parameter (W_NumContParam=27)

c --- values for A_MatrSolvMeth
      integer     W_SvKcl,     W_SvSparse, W_SvWatson
      parameter  (W_SvKcl = 1, W_SvSparse = 2, W_SvWatson = 3)
      
c --- values for element
      integer     W_eSlave,         W_eMaster, 
     &            W_eBeam2nd,       W_eBeam3nd,
     &            W_eLink,          W_eBeamEB
      parameter  (W_eSlave  = 101,  W_eMaster = 102, 
     &            W_eBeam2nd = 103, W_eBeam3nd = 104,
     &            W_eLink = 100,    W_eBeamEB = 105)     

c --- values for M_DofBits
      integer     W_ForcBit,      W_DispBit,     
     &            W_CpBit,        W_TranBit
      parameter  (W_ForcBit  = 1, W_DispBit = 2, 
     &            W_CpBit = 3,    W_TranBit = 4)

c --- output control information
c --- define the default bit settings for nodal and element results
      integer    W_ElResMask,       W_NdResMask,    W_AllResMask
      parameter (W_ElResMask=1016,  W_NdResMask=6,  W_AllResMask=1022)


c --- define the item bit pointers for output:
      integer    W_OutNdSol,    W_OutReact,     W_OutElIndx,    
     &			 W_OutStress,   W_OutElaStra,   W_OutPlStra,   
     &			 W_OutCrStra,   W_OutGyStra,    W_OutMisc
      
      parameter (W_OutNdSol=1,  W_OutReact=2,   W_OutElIndx=3,  
     &			 W_OutStress=4, W_OutElaStra=5, W_OutPlStra=6,
     &			 W_OutCrStra=7, W_OutGyStra=8,  W_OutMisc=9)

c --- SMP locks (if you add locks, change P_MaxLock in smp_info.inc)
      integer    
     & W_LockAutoStep, W_LockSvr,     W_LockElemInq, W_LockText,
     & W_LockMaxStra,  W_LockErr,     W_LockRst,     W_LockBin,
     & W_LockCut,      W_LockConvVar, W_LockLoadVect,W_LockSparse,
     & W_LockTrans,    W_LockAnsys,   W_LockTab
c --- in case of new locks, adjust P_MaxLock in smp_info.inc
	  parameter  (
     & W_LockElemInq=1,  W_LockLoadVect=2,   W_LockConvVar=3,
     & W_LockAutoStep=4, W_LockSvr=5,        W_LockText=6,
     & W_LockErr=7,      W_LockBin=8,        W_LockCut=9,
     & W_LockMaxStra=10, W_LockRst=11,       W_LockSparse=12,
     & W_LockTrans=13,   W_LockAnsys=14,     W_LockTab=15)

c	W_LockElemInq  - lock element inquiries
c	W_LockLoadVect - lock load vector
c	W_LockConvVar  - lock convergence data
c	W_LockAutoStep - lock auto load data
c	W_LockSvr      - lock SVR file
c	W_LockErr      - lock error handling
c	W_LockFull     - lock FULL file
c	W_LockCut      - lock for C_CutEle,C_CutMat
c	W_LockMaxStra  - lock for max. total strain
c	W_LockRst      - lock for RST file
c	W_LockBin      - lock for binary i/o
c	W_LockText     - lock for text i/o
c	W_LockSparse   - lock for global sparse matrix

c --- Sparse matrix lib, header index
      integer     W_SpLength,    W_SpNRows,     W_SpNTerms,
     &            W_SpNLoads,    W_SpNDofs,     W_SpStoreMeth,
     &            W_SpPtrMat,    W_SpPtrDiag,   W_SpPtrCol,    
     &            W_SpPtrLv,     W_SpPtrRowEnd
     &            
      parameter  (W_SpLength=10, W_SpNRows=1,    W_SpNTerms=2,
     &            W_SpNLoads=3,  W_SpNDofs=4,    W_SpStoreMeth=5,
     &            W_SpPtrMat=6,  W_SpPtrDiag=7,  W_SpPtrCol=8,
     &            W_SpPtrLv=9,   W_SpPtrRowEnd=10)

c	W_SpNRows      # of rows in matrix 
c	W_SpNTerms     # of terms in index array and double array
c	W_SpNLoads     # of solutions being processed
c	W_SpNDofs      maximum number of DOFs per node
c	W_SpStoreMeth  0 - store in rows, 1 - in columns
c	W_SpPtrMat     pointer to matrix(s) array
c	W_SpPtrDiag    pointer to iDiag position vector
c	W_SpPtrCol     pointer to iCol index vector
c	W_SpPtrLv      pointer to load vector
c	W_SpPtrRowEnd  pointer to end of each row (= W_SpPtrDiag)

c --- definition of sizes of nonlinear material tables
      integer    W_nMatCol, W_nMatRow
      parameter (W_nMatCol=2, W_nMatRow=4)

c**************************************************************
c --- PARAMETERS BELOW ARE NOT CHANGABLE IN THE LINKED VERSION
c --- THEY ARE RELATED TO DATA BASE QUERIES

c --- for data base inquire functions
      integer   DB_SELECTED, DB_MAXDEFINED, DB_MAXRECLENG
      parameter(DB_SELECTED=1, DB_MAXDEFINED=14, DB_MAXRECLENG=15)

c --- define symbolic names for element result file records
      integer    EDEMS,     EDENS,     EDEEL,     EDEPL,      
     &           EDECR,     EDETH,     EDENL
     
      parameter (EDEMS = 1, EDENS = 3, EDEEL = 6, EDEPL = 7,
     &           EDECR = 8, EDETH = 9, EDENL =15)


c --- properties definitions
      integer IELCSZ, IETYP, JETYP,
     & KYOP1, KYOP2, KYOP3, KYOP4, KYOP5, KYOP6 ,
     & KYOP7, KYOP8, KYOP9, KYO10, KYO11, KYO12,
     & KDIM,  ISHAP, MNODE, KDOFS, NEWGEN, NMDRLC, 
     & MATRQD, NCOMPN, NMNDMX, NMNDMN, NMNDST, NMDFPN, 
     & MATRXS, NMNDNE, JSECT, NMPRES, NMMOIS,
     & NMSFST, NMSMIS, NMSSVR, NMPSVR, NMCNSVR, KSURF,
     & JBEAM, NMSTVR

c --- size of characters array
      parameter (IELCSZ = 200)

c --- define element characteristic indices
      parameter  (
     & IETYP = 1, JETYP = 2,
     & KYOP1 = 3, KYOP2 = 4,KYOP3 = 5, KYOP4 = 6,  KYOP5 = 7,
     & KYOP6 = 8, KYOP7 = 9,KYOP8 = 10,KYOP9 = 11, KYO10 = 12, 
     & KYO11 = 13,KYO12 = 14)

      parameter (
     & KDIM   = 21, ISHAP  = 22, MNODE  = 25, KDOFS  = 34,
     & NEWGEN = 38, NMDRLC = 39, MATRQD = 47,
     & NCOMPN = 52, NMNDMX = 61, NMNDMN = 62, NMNDST = 63,
     & NMDFPN = 64, MATRXS = 66, NMNDNE = 67, JSECT  = 71, 
     & NMPRES = 75, NMMOIS = 83, NMSFST =108, NMSMIS =109, 
     & NMSSVR =122, NMPSVR =127, NMCNSVR=130, KSURF  =147,
     & JBEAM  =153, NMSTVR =160)

c --- degrees of freedom
      integer    UX,  UY,  UZ,  ROTX,  ROTY,   ROTZ
	  parameter (UX=1,UY=2,UZ=4,ROTX=8,ROTY=16,ROTZ=32)

c --- identify matrix capability with parameters
      integer    STIFM,   MASSM,   DAMPM,   SSTIFM
      parameter (STIFM=2, MASSM=4, DAMPM=8, SSTIFM=16)

c --- user TB properties
      integer    NLMP_CREEP, NLMP_USERMAT, NLMP_USTATE, 
     &           NLMP_PLAST, NLMP_CYCPLAST
      parameter (NLMP_CREEP = 2, NLMP_USERMAT = 13, 
     &           NLMP_USTATE = 25, NLMP_PLAST = 1, 
     &           NLMP_CYCPLAST = 24) 

c --- element integer data definitions
      integer   EL_MAT,     EL_TYPE,    EL_REAL, 
     &          EL_SHAPE,   EL_DEAD,    EL_DIM
      parameter(EL_MAT=1,   EL_TYPE=2,  EL_REAL=3,
     &          EL_DEAD=6,  EL_SHAPE=8, EL_DIM=10)