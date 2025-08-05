      subroutine SubStepIni
c --- primary function: initializations for the start of a new substep
c
#include "impcom.inc"
#include "solu_info.inc"
#include "conv_stat.inc"
#include "OVR_file.inc"
c
#if defined(DEBUG)
      call SubBeg('SubStepIni')
#endif

c --- update time
      A_TimeCur = A_TimeCur + A_TimeInc
c --- truncate roundoff beyond end time (possible with auto time step)
      if (A_TimeCur.gt.A_TimeEnd) A_TimeCur = A_TimeEnd

c --- set key for nonlinear stabilization
      if (A_KeyStab.gt.0) then
          if (A_SubStepNum.eq.1) then
              if (A_StabOpt.eq.2) then
                  A_StabFirstStep = 3
              else
                  A_StabFirstStep = 1
              endif
          else
              A_StabFirstStep = 0
          endif
      endif

c --- protect against zero divide
      if (A_TimeEnd .gt. A_TimeBeg) then
          A_CurLoadFact = (A_TimeCur-A_TimeBeg)/(A_TimeEnd-A_TimeBeg)
      else
          A_CurLoadFact = 1.0d0
      endif


c --- initialize auto load gap time and equilibrium convergence estimates
      if (A_KeyAutoStep.eq.1) then
          C_tGap = -1.0d0
          C_nIterEsti1b = 0
          C_nIterEsti2b = 0
      endif

c --- copy SVR to OVR for possible restart upon bisection
      if ((O_SaveKey.eq.1 .and. A_BisecNum.eq.0) .and.
     &    (A_LoadStepNum.gt.1 .or. A_SubStepNum.gt.1)) call SvrToOvr

      A_TotReactCur = 0.0d0

c --- update transient integration constants
      if (A_Trans.gt.0) call UpdTransParams
c --- update gravity
      call UpdGrav
      
#if defined(DEBUG)
      call SubEnd('SubStepIni')
#endif

      return
      end
