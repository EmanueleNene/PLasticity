      subroutine SetOutDefault (nelem,nnodes,entmsk,leng)
c --- primary function: set default bit settings for output controls
c
c  input arguments:
c     nelem    (int,sc,in)          - number of elements in ordered set
c     nnodes   (int,sc,in)          - number of nodes in compressed set
c
c  output arguments:
c     entmsk   (int,ar(leng),inout) - output control mask (leng is max(nelem,nnodes))
c
#include "impcom.inc"
#include "constants.inc"
#include "RST_file.inc"
#include "in_mem.inc"
c
      integer n,leng,nelem,nnodes,entmsk(leng)

#if defined(DEBUG)
      call SubBeg('SetOutDefault')
#endif

c --- initialize the masking per node and element entity if space allocated
      if (M_PtrElemOutMask .ne. PTRFTNNULL) then
          n = min (nelem,nnodes)
          call vi_init (entmsk(1),n,W_AllResMask)
          if (nelem .eq. nnodes) goto 999
          if (nelem .gt. nnodes) then
              call vi_init (entmsk(n+1),nelem-n,W_ElResMask)
          else
              call vi_init (entmsk(n+1),nnodes-n,W_NdResMask)
          endif
      endif

c --- initialize the global output mask
 999  R_OutMask = W_AllResMask


#if defined(DEBUG)
      call SubEnd('SetOutDefault')
#endif

      return
      end
