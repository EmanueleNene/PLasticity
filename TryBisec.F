      function TryBisec (AllowBisec,kcut)
c --- primary function: routine to try bisection if not converged
c
c  input arguments:
c     AllowBisec    (int,sc,in)
c
c  output arguments:
c     TryBisec  (int,sc,out)  - return statement
c     kcut      (int,sc,out)  - kcut for bisection
c
#include "impcom.inc"
#include "constants.inc"
#include "solu_info.inc"
c
      integer            TryBisec, AllowBisec, kcut
      double precision   dperr(2)



#if defined(DEBUG)
      call SubBeg('TryBisec')
#endif

      TryBisec = 0

      if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
          call IfBisect (3,kcut)
          TryBisec = 150
          goto 999
      endif
      dperr(1) = A_TimeCur
      dperr(2) = dble(A_SubStepNum)
      call erhandler('TryBisec',5000,2,
     &    'Solution failed to converge at time %g (substep %i).',
     &     dperr(1),' ')

      call StopAll (1)
      TryBisec = 999

  999 continue

#if defined(DEBUG)
      call SubEnd('TryBisec')
#endif

      return
      end
