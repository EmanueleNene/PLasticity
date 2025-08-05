      subroutine In_ExDP (ivin,ivext,n)
c --- primary function: internal format to exter format for double precision
c
#include "impcom.inc"
c
#if defined (LITTLE_ENDIAN)
      integer   i, n

      integer*4 ivin(*), ivext(*), iwork, iwork2

      do i = 1,2*n,2
         iwork = ivin(i)
         call SwapBytes (iwork)
         iwork2 = ivin(i+1)
         call SwapBytes (iwork2)
         ivext(i) = iwork2
         ivext(i+1) = iwork
      enddo
#else
      integer  i, n

      integer*4 ivin(*), ivext(*)

      do i = 1,2*n
         ivext(i) = ivin(i)
      enddo

#endif
      return
      end
