      subroutine v_zeroi (v,inc,n)
c --- primary function: set a vector to zero (with increments)
c
c  input arguments:
c     inc      (int,sc,in)         - increment (first dimension) of vector v
c     n        (int,sc,out)        - number of times to be repeated
c
c  output arguments:
c     v        (dp,ar(inc,n),out)  - vector v
c
#include "impcom.inc"
c
      integer  inc,n, i,ia
      double precision  v(*)
c
      ia = 1
      do i = 1,n
         v(ia) = 0.0d0
         ia = ia + inc
      enddo

      return
      end
