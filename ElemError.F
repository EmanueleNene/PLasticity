      function ElemError (AllowBisec,kcut)
c --- primary function: process error from element formation
c
c  output arguments:
c     ElemError (int,sc,out)     - return statement number
c     kcut      (int,sc,out)     - kcut for bisection
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
c
      integer   ElemError, AllowBisec, kcut, GetErr

#if defined(DEBUG)
      call SubBeg('ElemError')
#endif

c --- error in element formulation: bisection or exit
      if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
c --- try bisection
          call IfBisect (5,kcut)
          call PutErr (0)
          ElemError = 150
      else
c --- get out of here
          if (A_KeyFirstPass.eq.0) then
c --- don't print out more than 1 error
              if (GetErr().eq.0 .and. C_CutMat.eq.0) then
	            call DistortMessage(3)
              endif
              call StopAll (8)
          endif
          ElemError = 999
      endif

#if defined(DEBUG)
      call SubEnd('ElemError')
#endif

      return
      end
