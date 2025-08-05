      subroutine InitHygMem()
#include "impcom.inc"      
#include "beam_hyg.inc"
#include "constants.inc"
      
      PTRFTN  MemAlloc
      integer InqSMP
      integer nProcs, i
      
      nProcs = InqSMP(0)

      i = H_DpSize*nProcs
      H_pDpVArs = MemAlloc(i,W_LenDouble,'H_pDpVArs       ')

      i = H_IntSize*nProcs
      H_pIntVArs = MemAlloc(i,W_LenInt,'H_pIntVArs      ')
      
      return
      end