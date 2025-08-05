      subroutine cv1_v1 (v,n,const)
c --- primary function: multiply a vector by a constant
c
c  input arguments:
c     n        (int,sc,in)       - length of vector n
c     const    (dp,sc,in)        - constant to multiply v
c
c  output arguments:
c     v        (dp,ar(n),inout)  - vector v
c
#include "impcom.inc"
c
      integer  n, i
      double precision   v(n),const
 
      do i = 1,n
         v(i) = const*v(i)
      enddo

      return
      end
