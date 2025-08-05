      subroutine LineSearchIni
c --- primary function: initializations for the start of a line search pass
c
#include "impcom.inc"
#include "solu_info.inc"
c
#if defined(DEBUG)
      call SubBeg('LineSearchIni')
#endif

c --- initialize the start of the iteration
      call SetFormKeys

c --- only form stiffness
      if (A_Trans .eq. 0) A_MatrReq(3) = 0
      A_MatrReq(4) = 0
c
c --- never update history variables or saved variables
      A_KeyMatUpd = 0
      A_KeySvrUpd = 0

#if defined(DEBUG)
      call SubEnd('LineSearchIni')
#endif

      return
      end
