      function SolvError (kcut)
c --- primary function: process error output from equation solution
c
c  output arguments:
c     SolvError (int,sc,out)    - the return statement number
c     kcut      (int,sc,out)    - kcut for bisection option
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer SolvError, kcut
      logical IfEqvDP,IfMinTime  

#if defined(DEBUG)
      call SubBeg('SolvError')
#endif

      SolvError = 0
      IfMinTime = IfEqvDP(A_TimeInc,A_MinTimeStep) .or. 
     &                   (A_TimeInc.lt.A_MinTimeStep)
      
      if (A_KeyFirstPass.eq.1 .or.
     &   (A_BisecNum.gt.0 .and. IfMinTime) ) then
c --- error on the first pass or on minimum step: abort
          SolvError = 999
      else
c --- try bisection second
          call IfBisect (2,kcut)
          call PutErr(0)
          SolvError = 150
      endif

#if defined(DEBUG)
      call SubEnd('SolvError')
#endif

      return
      end
