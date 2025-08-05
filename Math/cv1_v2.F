      subroutine cv1_v2 (v1,v2,n,const)
c --- primary function: multiply a vector by a constant
c
c  input arguments:
c     v1       (dp,ar(n),in)     - vector v1
c     n        (int,sc,in)       - length of vectors v1, v2
c     const    (dp,sc,in)        - constant to multiply v1
c
c  output arguments:
c     v2       (dp,ar(n),out)    - vector v2
c
#include "impcom.inc"
c
      integer  n, i
      double precision   v1(n),v2(n),const
 
      do i = 1,n
         v2(i) = const*v1(i)
      enddo

      return
      end
