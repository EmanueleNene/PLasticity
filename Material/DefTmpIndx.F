      subroutine DefTmpIndx (nTens, firstIndex, matTmpIndex)
c --- primary function: define work array index
c
#include "impcom.inc"
#include "constants.inc"
c
      integer          firstIndex, nTens, matTmpIndex (W_LenWrkIndx)
      integer          nTensSq

      nTensSq    = nTens * nTens

      firstIndex = 1
      matTmpIndex (W_WrkIndx1) = firstIndex
      matTmpIndex (W_WrkIndx2) = firstIndex + nTens
      matTmpIndex (W_WrkIndx3) = firstIndex + nTens*2
      matTmpIndex (W_WrkIndx4) = firstIndex + nTens*3
      matTmpIndex (W_WrkIndx5) = firstIndex + nTens*4
      matTmpIndex (W_WrkIndx6) = firstIndex + nTens*7
      matTmpIndex (W_WrkIndx7) = firstIndex + nTens*8
      matTmpIndex (W_WrkIndx8) = firstIndex + nTens*9
      matTmpIndex (W_WrkIndx9) = firstIndex + nTens*10
      firstIndex               = firstIndex + nTens*10 + nTensSq

      return
      end
