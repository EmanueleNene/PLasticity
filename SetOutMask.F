      subroutine SetOutMask (OutCode,kk,nnodes,nelem,entmsk)
c --- primary function: set or clear bit for output control mask
c
c  input arguments:
c     OutCode   (int,sc,in)      - output item number
c     kk       (int,sc,in)       - loop counter
c     nnodes   (int,sc,in)       - number of nodes
c     nelem    (int,sc,in)       - number of elements
c     entmsk   (int,sc,inout)    - output mask
c
c  output arguments:
c     entmsk   (int,sc,inout)    - output mask
c
#include "impcom.inc"
#include "constants.inc"
c
      integer nnodes,nelem,kk,entmsk,OutCode
c
#if defined(DEBUG)
      call SubBeg('SetOutMask')
#endif

      if (OutCode .eq. 101) then
c --- write all
          entmsk = W_AllResMask
      elseif (OutCode .eq. 102) then
c --- write basic
          if (kk .le. nnodes) then
              entmsk = ibset(entmsk,W_OutNdSol)
              entmsk = ibset(entmsk,W_OutReact)
          endif
c --- no longer included in basic output
c          if (kk .le. nelem) then
c              entmsk = ibset (entmsk,W_OutElIndx)
c              entmsk = ibset (entmsk,W_OutStress)
c          endif
      else
c --- print or write a specific item
          entmsk = ibset (entmsk,OutCode - 1)
      endif

#if defined(DEBUG)
      call SubEnd('SetOutMask')
#endif
      return
      end
