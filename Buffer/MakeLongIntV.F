      subroutine MakeLongIntV (n,iVectL)
c primary function: convert internal long integer to external
c
c  input arguments:
c     n         (int,sc,in)         - number of long integers
c
c  output arguments:
c     iVectL    (int*8,ar(n),inout) - long integer vector
c
#include "impcom.inc"
c
      integer n
      integer*8 iVectL(n)

      integer    i, j, iwork(2)
      integer*8  iworkL
      equivalence (iworkL,iwork(1))

#if defined (IEEE_SYS)
      do i = 1,n
         iworkL = iVectL(i)
         j = iwork(1)
         iwork(1) = iwork(2)
         iwork(2) = j
         iVectL(i) = iworkL
      enddo
#endif

      return
      end
