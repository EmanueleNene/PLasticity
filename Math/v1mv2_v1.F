      subroutine v1mv2_v1 (v1,v2,n)
c --- primary function: subtract vector v2 from vector v1 and save in vector v1
c
c  input arguments:
c     v2       (dp,ar(n),in)     - v2 vector
c     n        (int,sc,in)       - length of vectors v1,v2
c
c  output arguments:
c     v1       (dp,ar(n),inout)  - v1 vector
c
#include "impcom.inc"
c
      integer  n, i
      double precision   v1(n),v2(n)
 
      do i = 1,n
         v1(i) = v1(i) - v2(i)
      enddo

      return
      end
