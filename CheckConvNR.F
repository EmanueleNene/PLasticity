      function CheckConvNR (AllowBisec,kcut)
c --- primary functions: check for force convergence
c
c  input arguments:
c     AllowBisec     (int,sc,in)   - key indicates bisection is possible
c
c  output arguments:
c     CheckConvNR  (int,sc,out)    - branch return path
c     kcut         (int,sc,out)    - kcut value
c
#include "impcom.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "in_mem.inc"
c
      integer   GetErr

      integer   CheckConvNR,AllowBisec,kcut,nondefnorm
      logical   CheckIfConvOK


#if defined(DEBUG)
      call SubBeg('CheckConvNR')
#endif

      CheckConvNR = 0
      nondefnorm = 0

c --- check for force convergence
      call CheckConvForc (M_IncI(1),M_ForcTot(1),M_ForcNR(1),
     &                    C_Forc,nondefnorm)

      
c --- check convergence pattern: A_KeyBisect the load if needed
      if (AllowBisec.gt.0 .and. C_Forc.eq.0
     &                    .and. A_KeyFirstPass.eq.0) then
          if (AllowBisec.eq.2 .and. A_IterNum.ge.3) then
c --- let it go through
              C_Forc = 1
              CheckConvNR = 700
              goto 999
          else
              call IfBisect (1,kcut)
              if (kcut.eq.1) then
                  CheckConvNR = 150
                  goto 999
              endif
          endif
      endif
c --- dirty fix to skip the contact removal in UpdContConn      
      if (nondefnorm.ne.0) C_nDeBo = -1
c --- it is very wrong or force norm is not a number (NAN) which is normally caused by predictor
      if (GetErr().eq.1 .or. nondefnorm.eq.1) then
c --- dof limit exceeded: A_KeyBisect or exit
          if (AllowBisec.eq.1 .and. A_KeyFirstPass.eq.0) then
c --- try bisection
              call IfBisect (4,kcut)
              call PutErr (0)
              CheckConvNR = 150
              goto 999
          else
c --- get out of here
              call StopAll (2)
              CheckConvNR = 999
              goto 999
          endif
      endif

c --- NOTE: we used A_ConvLineSrch before to pass over line search status
c           if solution is converged but line search is not,
c           this means there could be some problems on the element level
c           we do not process this right now, but we probably will in future

      if (CheckIfConvOK()) then
          CheckConvNR = 700
          goto 999
      else
c --- if the only thing that did not converge is the linesearch, than send it back
          if (A_IterNum.eq.(A_IterMax+1) .and. 
     &        C_Forc.eq.1 .and. C_Disp.eq.1) then
              CheckConvNR = 300
              goto 999
          endif
      endif

      
      if (A_IterNum.ge.(A_IterMax+1) .and. A_IterNum.gt.1) then
          A_IterNum = A_IterNum - 1
          A_CumIterNum = A_CumIterNum - 1
          CheckConvNR = 650
          goto 999
      endif

  999 continue

#if defined(DEBUG)
      call SubEnd('CheckConvNR')
#endif

      return
      end
