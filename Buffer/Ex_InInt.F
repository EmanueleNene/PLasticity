      subroutine Ex_InInt (intext,intin,n)
c --- primary function: exter format to internal format for integers
c
#include "impcom.inc"
c
#if defined (LITTLE_ENDIAN)
      integer  i,n
#if defined (FULL_64)
      integer*4  intext(*), iwork
      integer*8  intin(*)

      do i = n, 1, -1
         intin(i) = intext(i)
         iwork = intin(i)
         call SwapBytes (iwork)
         intin(i) = iwork
      enddo
#else
      integer  intext(*),intin(*),iwork
      do i = n, 1, -1
         iwork = intext(i)
         call SwapBytes (iwork)
         intin(i) = iwork
      enddo
#endif
#else

      integer  n, intext(*), intin(*)
      call vi_move(intin(1),intext(1),n)

#endif
      return
      end
