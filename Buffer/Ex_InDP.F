      subroutine Ex_InDP (ivext,ivin,n)
c --- primary function: exter format to internal format for double precision
c
#include "impcom.inc"
c
#if defined (LITTLE_ENDIAN)      
      integer  n, i

      integer*4  ivin(*), ivext(*), iwork, iwork2

      do i = 1,2*n,2
         iwork = ivext(i)
         iwork2 = ivext(i+1)
c --- swap the bytes
         call SwapBytes (iwork)
         call SwapBytes (iwork2)
c --- swap the words
         ivin(i) = iwork2
         ivin(i+1) = iwork
      enddo
#else
      integer  i,n

      integer*4  ivin(*), ivext(* )

      do i = 1,2*n
         ivin(i) = ivext(i)
      enddo

#endif
      return
      end
