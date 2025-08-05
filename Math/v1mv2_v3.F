      subroutine v1mv2_v3 (v1,v2,v3,n)
c --- primary function: subtract vector v2 from vector v1 to get vector v3
c
c  input arguments:
c     v1       (dp,ar(n),in)     - vector v1
c     v2       (dp,ar(n),in)     - vector v2
c     n        (int,sc,in)       - length of vectors v1,v2,c
c
c  output arguments:
c     v3       (dp,ar(n),out)    - vector v3
c
c
#include "impcom.inc"
c
      integer  n, i
      double precision  v1(n),v2(n),v3(n)
 
      do i = 1,n
         v3(i) = v1(i) - v2(i)
      enddo

      return
      end
