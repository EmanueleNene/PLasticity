      subroutine SetCutBackFact (ishape)
c --- primary function: update cutback and error keys from elements
c
#include "impcom.inc"
#include "conv_stat.inc"
#include "constants.inc"
c
      integer ishape
      double precision cutfactor


      if (ishape.ne.0) then
          call lockSMP (W_LockCut)
          C_CutEle = ishape
          call unlockSMP (W_LockCut)
      endif

      return
      end
