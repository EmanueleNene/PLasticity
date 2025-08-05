      subroutine SolverEnd
c --- subroutine to finish solution routine
c
#include "impcom.inc"
#include "constants.inc"
#include "in_mem.inc"
#include "solu_info.inc"
#include "SVR_file.inc"
#include "OVR_file.inc"
#include "RST_file.inc"
#include "dof_info.inc"
c
      integer InqSMP
c
#if defined(DEBUG)
      call SubBeg('SolverEnd')
#endif

      call Bin_Clo (S_BlockNum,'DELE')
      
c --- close results file buffer
      if (R_Unit.gt.0 .and. R_BlockNum.gt.0) then
          call Bin_Clo(R_BlockNum,'KEEP')
      endif
      R_Unit  = 0
      R_FilePos = 0

c --- deallocate memory used during solution
      call DeallocSolMem()

c --- zero corresponding lengths
      D_ULen = 0
      D_BacLen = 0
      D_FwdLen = 0
      D_NumEl = 0

c --- clean material memory
      call MemFree (M_PtrMatRecDp)
      call MemFree (M_PtrMatRecI)

c --- terminate parallel
      if (InqSMP(7).gt.0) call FiniSMP

#if defined(DEBUG)
      call SubEnd('SolverEnd')
#endif

      return
      end
