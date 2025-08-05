      subroutine cv1_v2i (v1,inc1,v2,inc2,n,factor)
c --- primary function: multiply a vector (with increments) by a factor
c
c  input arguments:
c     v1       (dp,ar(inc1,n),in)  - vector v1
c     inc1     (int,sc,in)         - increment (first dimension) of vector v1
c     inc2     (int,sc,in)         - increment (first dimension) of vector v2
c     n        (int,sc,in)         - length (second dimension) of vector v1,v2
c     factor   (dp,sc,in)          - factor to multiply vector v1
c
c  output arguments:
c     v2       (dp,ar(inc2,n),out) - vector v2
c
#include "impcom.inc"
c
      integer  inc1, inc2, n, i, i1, i2
      double precision  v1(*),v2(*),factor
c
      i1 = 1
      i2 = 1
      do  i = 1,n
         v2(i2) = v1(i1) * factor
         i1 = i1 + inc1
         i2 = i2 + inc2
      enddo
c
      return
      end
