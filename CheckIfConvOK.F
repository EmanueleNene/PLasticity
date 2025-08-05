      function CheckIfConvOK
c --- primary functions: check for force convergence
c
c  output arguments:
c     CheckIfConvOK (int,sc,out)    - if convergence really achived
c
#include "impcom.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "dof_info.inc"
#include "in_mem.inc"

      integer GetWriteUnit
      integer iott
      logical CheckIfConvOK, IfMinTime, IfEqvDP

#if defined(DEBUG)
      call SubBeg('CheckIfConvOK')
#endif

      CheckIfConvOK = .False.
      IfMinTime = IfEqvDP(A_TimeInc,A_MinTimeStep) .or. 
     &                   (A_TimeInc.lt.A_MinTimeStep)
      
      if (.not.(A_BisecNum.gt.0 .or. IfMinTime)) then
          if (A_KeyPlast.gt.0 .and. A_StraMax.gt.A_StraLimit) goto 999
      endif
      
      if (C_Forc.eq.1 .and. (C_Disp.eq.1 .and. A_ConvLineSrch)) then
c --- solution is converged!
          A_KeyConv = 1
          CheckIfConvOK = .True.
c --- send a null displacement increment to the output
c     we've already done the stress calculations
          call v_zero (M_IncI(1),D_ULen)
          A_CumIterNum = A_CumIterNum - 1
          A_IterNum = A_IterNum - 1

          iott = GetWriteUnit()
          write (iott,2000) A_IterNum, A_CumIterNum
 2000     format ('--> SOLUTION CONVERGED AFTER ITERATION',i4,
     &            ' (TOTAL',i6,')'/)
      endif

999   continue

#if defined(DEBUG)
      call SubEnd('CheckIfConvOK')
#endif

      return
      end
