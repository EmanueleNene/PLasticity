      subroutine StepIni
c --- primary function: setup for upcoming load step
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "dof_info.inc"
#include "solu_info.inc"
c
      integer  cpinqr
      PTRFTN   MemAlloc

      integer  leng

#if defined(DEBUG)
      call SubBeg('StepIni')
#endif

c --- create workspace for beams
      call WorkSpacesIni

      D_NumCp=cpinqr(0,DB_MAXDEFINED)
      D_LenCp=cpinqr(0,DB_MAXRECLENG)
c --- coupling slave-to-master correspondence
      if (D_NumCp.gt.0) then
          leng = D_NumNd*D_NumDof
          M_PtrMasterCp = MemAlloc(leng,W_LenInt,'M_PtrMasterCp   ')
          call vi_zero (M_MasterCp(1),leng)
      endif

      call MakeCpVect

#if defined(DEBUG)
      call SubEnd('StepIni')
#endif

      return
      end
