      subroutine In_ExInt (intin,intext,n)
c --- primary function: internal format to exter format for integers
c
#include "impcom.inc"
c
#if defined (LITTLE_ENDIAN)
      integer  n, i

#if defined (FULL_64)
      integer*4  intext(*), iwork
      integer*8  intin(*), j

c --- convert at the end so as to not overwrite good data before conversion
      j = 2*n
      do i = n,1,-1
         intext(j) = intin(i)
         iwork = intext(j)
         call SwapBytes (iwork)
         intext(j) = iwork
         j = j - 1
      enddo

c --- slide data to the beginning
      do i = 1,n
         intext(i) = intext(i+n)
      enddo

#else
      integer  intext(*),intin(*),iwork

      do i = n,1,-1
         iwork = intin(i)
         call SwapBytes (iwork)
         intext(i) = iwork
      enddo

#endif
#else
      integer  n,intin(n),intext(n)
      
      call vi_move(intin(1),intext(1),n)

#endif

      return
      end
