      subroutine DeallocSolMem
c --- primary function: deallocate memory used during solution
c
c     NOTE: Any memory added in AllocSolMem.F MUST BE DEALLOCATED inside this routine
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "cont_info.inc"
#include "dof_info.inc"
#include "solu_info.inc"
#include "solver_info.inc"
#include "sec_info.inc"
#if defined(MOIST)
#include "beam_hyg.inc"
#endif
c
#if defined (PARDISO)
      integer i, phase, error, idummy, iparm(64), GetWriteUnit, msglvl
      SP_INT  ldummy
      double precision dummy
#endif      
      
#if defined(DEBUG)
      call SubBeg('DeallocSolMem')
#endif
      
c --- free space for n-r force vector
      call MemFree (M_PtrForcNR)

c --- free space for force vector due to stabilization
      call MemFree (M_PtrForcStab)

c --- free space for total force vector
      call MemFree (M_PtrForcTot)

c --- free space for u vectors
      if (M_PtrSoluVects.ne.PTRFTNNULL) then
          call MemFree (M_PtrSoluVects)
          M_PtrDisp = PTRFTNNULL
          M_PtrIncI = PTRFTNNULL
          M_PtrIncS = PTRFTNNULL
          M_PtrVelo = PTRFTNNULL
          M_PtrAcce = PTRFTNNULL
      endif

c --- free space for global rotations
      call MemFree (M_PtrERot)

c --- free DOF Bits for constraints, forces, masters, etc
      call MemFree (M_PtrDofBits)

c --- free internal to external nodal pointer
      call MemFree (M_PtrBackList)

c --- free external to internal nodal pointer
      call MemFree (M_PtrFrwdList)

c --- free element order vectors
      call MemFree (M_PtrElemOrd)

c --- free space for contact database
      if (K_nMast.gt.0 .or. K_nFoun.gt.0) then
          call MemFree (M_PtrElSet)
          call MemFree (K_PtrFounDp)
          call MemFree (K_PtrSeg)
          call MemFree (K_PtrContInfo)
          call MemFree (K_PtrPiloNd)
          call MemFree (K_PtrConSrch)
      endif

c --- deallocate coupling slave-to-master correspondence
      if (D_NumCp .gt. 0) then
          call MemFree (M_PtrMasterCp)
          call MemFree (M_PtrSlaveCpDofs)
      endif
      
c --- free SVR in-memory space
      call MemFree(M_PtrOvr)

c --- free working memory for elements
      call MemFree(M_PtrElemWrkI)
      call MemFree(M_PtrElemWrkDp)

c --- deallocate elements workspace
      call MemFree (M_PtrMatRecDp)
      call MemFree (M_PtrMatRecI)

c --- deallocate space with SVRs
      call MemFree (M_PtrSvrD)
      call MemFree (M_PtrSvrIndx)

c --- cleanup solution database objects
      call iniAnsObj(1)

c --- deallocate sparse object
      call SpDealloc (M_PtrSparseMtrx)
#if defined(WATSON)
      if (abs(A_MatrSolvMeth) .eq. W_SvWatson) then
          call MemFree (P_iperm)
          call MemFree (P_invp)
          call wsmp_clear
      endif
#endif

#if defined (PARDISO)
      call mkl_free_buffers()
c --- this functionality may cause crash when the called for the second time
c      if (abs(A_MatrSolvMeth) .eq. W_SvSparse) then      
c          phase = -1 ! release internal memory
c          dummy = 0.0d0
c          idummy = 0
c          ldummy = 0
c          call vi_zero(iparm(1),64)
c          iparm(1) = 1
c          iparm(7) = GetWriteUnit()
c          msglvl = 1 ! print statistical information          
c          call pardiso (P_pt(1),idummy,idummy,idummy,phase,idummy,
c     &                  dummy,ldummy,idummy,
c     &                  idummy,idummy,iparm(1),msglvl,
c     &                  dummy,dummy,error)
c      endif
#endif

#if defined(MOIST)
      if (K_nMast.gt.0) then
          call MemFree (M_PtrBeams)
          call MemFree (M_PtrHygStrains)
          call MemFree (H_pBCData)
          call MemFree (H_pDpVArs)
          call MemFree (H_pIntVArs)
      endif
#endif
      call MemFree(B_PtrRealN)
      call MemFree(B_PtrSectData)
      
#if defined(DEBUG)
      call SubEnd('DeallocSolMem')
#endif

      return
      end