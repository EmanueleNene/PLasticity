      subroutine v_movei (v1,inc1,v2,inc2,n)
c --- primary function: move one vector into another (with increment)
c
c  input arguments:
c     v1       (dp,ar(inc1,n),in)  - vector v1
c     inc1     (int,sc,in)         - increment on vector v1
c     inc2     (int,sc,in)         - increment on vector v2
c     n        (int,sc,in)         - number of items to be moved
c
c  output arguments:
c     v2       (dp,ar(inc2,n),in)  - vector v2
c
#include "impcom.inc"
c
      integer  inc1,inc2,n, i,i1,i2
      double precision  v1(*),v2(*)
 
      i1 = 1
      i2 = 1
      do i = 1,n
         v2(i2) = v1(i1)
         i1 = i1 + inc1
         i2 = i2 + inc2
      enddo
c
      return
      end
