      subroutine v1pcv2_v1 (v1,v2,n,const)
c --- primary function: multiply vector v2 to constant, add to vector v1,
c                       and store in vector v1
c
c  input arguments:
c     v2       (dp,ar(n),in)     - v2 vector
c     const    (int,sc,in)       - constant
c     n        (int,sc,in)       - length of vectors v1,v2
c
c  output arguments:
c     v1       (dp,ar(n),inout)  - v1 vector  
c
#include "impcom.inc"

      integer  n, i
      double precision  v1(n),v2(n),const
 
      do i = 1,n
         v1(i) = v1(i) + v2(i)*const
      enddo

      return
      end
