      subroutine AllocSolMem(kerr)
c --- primary function: allocate memory needed during solution
c
c     NOTE: Any memory added here MUST BE DEALLOCATED inside
c           the companion routine DeallocSolMem.F
c
c  output arguments:
c     kerr     (int,sc,out)      - error flag   0: OK, 1: error
c
#include "impcom.inc"
#include "constants.inc"
#include "dof_info.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "solver_info.inc"
#if defined(MOIST)
#include "cont_info.inc"
#endif
c
      PTRFTN    MemAlloc, MemAllocL

      integer*8 iL
      integer   kerr,ndinqr,i


#if defined(DEBUG)
      call SubBeg('AllocSolMem')
#endif

      D_NumNd = ndinqr(0,DB_MAXDEFINED)

c --- external to internal nodal pointer
      M_PtrFrwdList = MemAlloc(D_NumNd,W_LenInt,'M_PtrFrwdList   ')
      D_FwdLen = D_NumNd

c --- internal to external nodal pointer
      M_PtrBackList = MemAlloc(D_NumNd,W_LenInt,'M_PtrBackList   ')
      call DefNodeVects (D_NumEl,M_ElemOrd(1),D_NumNd,D_BacLen)

      if (kerr.ne.0) goto 999

      D_NumDof = 0
      do i = 1,6
         D_NumDof = D_NumDof + 1
         D_DofList(D_NumDof) = i
      enddo
      D_DofList(7) = 0

c --- length of DOF vector
      D_ULen = D_BacLen*D_NumDof

c --- DOF bits for constraints, forces, masters, etc
      M_PtrDofBits = MemAlloc(D_ULen,W_LenInt,'M_PtrDofBits    ')
      call vi_zero (M_DofBits(1),D_ULen)

c --- get space for u vectors
      iL = 4
      if (A_Trans.gt.0) iL = iL + 1
      iL = iL*D_ULen
      M_PtrSoluVects = MemAllocL(iL,W_LenDouble,'M_SoluVects     ')
      M_PtrDisp = M_PtrSoluVects
      M_PtrIncI = M_PtrDisp + D_ULen*W_AddrStrideDp
      M_PtrIncS = M_PtrIncI + D_ULen*W_AddrStrideDp
      M_PtrVelo = M_PtrIncS + D_ULen*W_AddrStrideDp
      call v_zero (M_Disp(1),D_ULen)
      call v_zero (M_IncI(1),D_ULen)
      call v_zero (M_IncS(1),D_ULen)
      call v_zero (M_Velo(1),D_ULen)
      
      if (A_Trans.gt.0) then
          M_PtrAcce = M_PtrVelo + D_ULen*W_AddrStrideDp
          call v_zero (M_Acce(1),D_ULen)
      else
          M_PtrAcce = PTRFTNNULL
      endif

c --- get space for global rotations
      M_PtrERot = MemAlloc(3*D_BacLen,W_LenDouble,'M_PtrERot       ')
      call v_zero (M_ERot(1,1),3*D_BacLen)

c --- get space for total force vector
      M_PtrForcTot = MemAlloc(D_ULen,W_LenDouble,'M_PtrForcTot    ')
      call v_zero (M_ForcTot(1),D_ULen)

c --- get space for n-r force vector
      M_PtrForcNR = MemAlloc(D_ULen,W_LenDouble,'M_PtrForcNR     ')
      call v_zero (M_ForcNR(1),D_ULen)

c --- get space for nonlinear stabilization
      if (A_KeyStab.gt.0) then
          M_PtrForcStab = MemAlloc(D_ULen,W_LenDouble,'M_PtrForcStab  ')
          call v_zero (M_ForcStab(1),D_ULen)
      else
          M_PtrForcStab = PTRFTNNULL
      endif

c --- get space for contact
      call ContactPreset (M_ElemOrd(1),D_BacLen)

c --- create working memory for new generation of elements
      call MakeElemWrkSpace
      M_PtrSvrD = PTRFTNNULL
      M_PtrSvrIndx = PTRFTNNULL

#if defined (PARDISO)
      if (abs(A_MatrSolvMeth) .eq. W_SvSparse) then      
          do i = 1,64
             P_pt(i) = 0
          enddo
      endif
#endif

#if defined(MOIST)
c --- get space for contact
      if (K_nMast.gt.0) then
          call BeamPreset (D_NumEl)
          call InitHygVars()
      endif
#endif
      
      call SectPreset(D_NumEl)
      
  999 continue
#if defined(DEBUG)
      call SubEnd('AllocSolMem')
#endif

      return
      end
