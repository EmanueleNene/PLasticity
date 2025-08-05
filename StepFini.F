      subroutine StepFini
c --- primary function: finish up a load step
c
#include "impcom.inc"
#include "solu_info.inc"
c
#if defined(DEBUG)
      call SubBeg('StepFini')
#endif

c --- set the beginning time of the next step to the end time of this step
      A_TimeBeg = A_TimeEnd

#if defined(DEBUG)
      call SubEnd('StepFini')
#endif

      return
      end
