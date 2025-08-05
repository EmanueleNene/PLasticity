      subroutine OutIni
c --- primary function: initialize the start of an output pass
c
#include "impcom.inc"
#include "solu_info.inc"
#include "dof_info.inc"
c
      logical  IfEqvDP

c --- unset the first pass key
      A_KeyFirstPass = 0
      A_KeyFirstIter = 0

c --- initialize the element keys
      call vi_zero (A_MatrReq(1),10)
      A_KeySvrUpd = 1

c --- set the output key
      A_KeyOut = 1

c --- set newton-raphson updating keys
      A_KeyMatUpd = 0
      A_KeyBCsUpd = 0
      A_KeyMoistUpd = 0
c --- only do if converged or not stopping the analysis
      if (A_KeyConv.eq.1 .or. A_IterMax.eq.1) then
          A_KeyMatUpd = 1
c --- update BCs at end of this load step to be BCs at start of next
          if (IfEqvDP(A_TimeCur,A_TimeEnd)) A_KeyBCsUpd = 1
c --- update the previous moisture info at end of this substep
          A_KeyMoistUpd = 1
      endif
c

      return
      end
