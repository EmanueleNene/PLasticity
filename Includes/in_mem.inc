      integer     M_Limit
#if defined(MOIST)
      parameter  (M_Limit=28)
#else
      parameter  (M_Limit=26)
#endif

      PTRFTN   M_PtrComm

      pointer (M_PtrDisp,M_Disp)
      pointer (M_PtrIncI,M_IncI)
      pointer (M_PtrIncS,M_IncS)
      pointer (M_PtrVelo,M_Velo)
      pointer (M_PtrAcce,M_Acce)
      pointer (M_PtrERot,M_ERot)
      pointer (M_PtrSvrIndx,M_SvrIndx)
      pointer (M_PtrElemOutMask,M_ElemOutMask)
      pointer (M_PtrMatRecI,M_MatRecI)
      pointer (M_PtrMatRecDp,M_MatRecDp)
      pointer (M_PtrMasterCp,M_MasterCp)
      pointer (M_PtrOvr,M_Ovr)
      pointer (M_PtrDofBits,M_DofBits)
      pointer (M_PtrElemOrd,M_ElemOrd)
      pointer (M_PtrFrwdList,M_FrwdList)
      pointer (M_PtrBackList,M_BackList)
      pointer (M_PtrSoluVects,M_SoluVects)
      pointer (M_PtrForcTot,M_ForcTot)
      pointer (M_PtrForcNR,M_ForcNR)
      pointer (M_PtrElemWrkI,M_ElemWrkI) 
      pointer (M_PtrElemWrkDp,M_ElemWrkDp)
      pointer (M_PtrForcStab,M_ForcStab)
      pointer (M_PtrElSet,M_ElSet)
      pointer (M_PtrSlaveCpDofs,M_SlaveCpDofs)
      pointer (M_PtrSparseMtrx,SparseMtrx)
#if defined(MOIST)
	  pointer (M_PtrBeams,M_Beams)
      pointer (M_PtrHygStrains,M_HygStrains)
#endif
      pointer (M_PtrSvrD,M_SvrD)

      integer 
     & M_ElemOutMask(*), M_MasterCp(*), M_SvrIndx(*), 
     & M_MatRecI(*), M_ElemOrd(*), M_FrwdList(*), M_BackList(*), 
     & M_ElemWrkI(*), M_DofBits(*), M_ElSet(*), M_SlaveCpDofs(*),
     & SparseMtrx(*), M_SvrD(*)
#if defined(MOIST)
     & ,M_Beams(*)
#endif

      double precision  
     & M_Disp(*), M_IncI(*), M_IncS(*), M_Velo(*),
     & M_Acce(*), M_ERot(3,*), M_MatRecDp(*),
     & M_Ovr(*), M_SoluVects(*), M_ForcTot(*), 
     & M_ForcNR(*), M_ElemWrkDp(*), M_ForcStab(*) 
#if defined(MOIST)
     & ,M_HygStrains(*)
#endif

      common /in_mem/ M_PtrComm(M_Limit)
	
      equivalence
     &  (M_PtrDisp,          M_PtrComm(1)),
     &  (M_PtrIncI,          M_PtrComm(2)),
     &  (M_PtrIncS,          M_PtrComm(3)),
     &  (M_PtrVelo,          M_PtrComm(4)),
     &  (M_PtrAcce,          M_PtrComm(5)), 
     &  (M_PtrERot,          M_PtrComm(6)) 


      equivalence
     &  (M_PtrMatRecI,       M_PtrComm(7)),
     &  (M_PtrMatRecDp,      M_PtrComm(8))

      equivalence
     &  (M_PtrSvrIndx,       M_PtrComm(9))

      equivalence
     &  (M_PtrElemOutMask,   M_PtrComm(10))


      equivalence
     &  (M_PtrMasterCp,      M_PtrComm(11)),
     &  (M_PtrOvr,           M_PtrComm(12))

      equivalence
     &  (M_PtrDofBits,       M_PtrComm(13)),
     &  (M_PtrElemOrd,       M_PtrComm(14)),
     &  (M_PtrFrwdList,      M_PtrComm(15)),
     &  (M_PtrBackList,      M_PtrComm(16)),
     &  (M_PtrSoluVects,     M_PtrComm(17)),
     &  (M_PtrForcTot,       M_PtrComm(18)),
     &  (M_PtrForcNR,        M_PtrComm(19))

      equivalence
     &  (M_PtrElemWrkI,      M_PtrComm(20)),
     &  (M_PtrElemWrkDp,     M_PtrComm(21)),
     &  (M_PtrForcStab,      M_PtrComm(22))

      equivalence
     &  (M_PtrElSet,	     M_PtrComm(23)),
     &  (M_PtrSlaveCpDofs,   M_PtrComm(24)),
     &  (M_PtrSparseMtrx,    M_PtrComm(25)),
#if defined(MOIST)
     &  (M_PtrBeams,         M_PtrComm(26)),
     &  (M_PtrHygStrains,    M_PtrComm(27)),
#endif
     &  (M_PtrSvrD,          M_PtrComm(M_Limit)	 )


c	M_Limit        - limit number of memory pointers

c	M_PtrComm()    - vector of memory pointers

c --- Memory pointers for global data

c	M_Disp         - total displacement location
c	M_IncI         - displacement increment, over the equilibrium iteration
c	M_IncS         - displacement increment, over the substep
c	M_Velo         - velocity
c	M_Acce         - acceleration (reserved)
c	M_ERot         - rotational DOFs Euler parameters


c	M_MatRecI      - element integer work space
c	M_MatRecDp     - element DP work space
c	M_Ovr          - in-core solver saved data

c	M_ElemOutMask  - output mask
c	M_MasterCp     - slave-to-master correspondence
c   M_SlaveCpDofs  - slave CP Dofs

c	M_DofBits      - bits for each dof for disp, force, master, etc.
c	M_ElemOrd      - element element order vector
c	M_FrwdList     - external to internal node vector
c	M_BackList     - internal to external node vector
c	M_ForcTot      - total applied force vector
c	M_ForcNR       - Newton-Raphson force vector
c	M_ForcStab     - force vector due to nonlinear stabilization

c	M_ElemWrkI     - element inter working array for the beam elements
c	M_ElemWrkDp    - element double precision working array for the beam elements
c	M_Beams        - array storing the numbers of beam elements
c	M_HygStrains   - nodal hygro-strains
c   M_PtrSvrD      - element incore SVR data