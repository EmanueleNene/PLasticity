      subroutine MtxReformCheck
c --- primary function: make ultimate determination whether to reform matrix
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
#if defined(DEBUG)
      call SubBeg('MtxReformCheck')
#endif

c --- force reform if first step of new analysis
      if (A_KeyFirstPass.eq.1) goto 200

c --- reform if there are any changes in the constraints
      if (A_SubStepNum.eq.1 .and. A_IterNum.eq.1) goto 200

c --- reform if any element matrices have changed
      if (A_KeyReform.eq.1) goto 200

c --- reuse the old matrix
 100  A_KeyReform = 0
      goto 999

c --- reform the matrix
 200  A_KeyReform = 1

  999 continue

#if defined(DEBUG)
      call SubEnd('MtxReformCheck')
#endif

      return
      end
