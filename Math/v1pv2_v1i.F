      subroutine v1pv2_v1i (v1,inca,v2,incb,n)
c --- primary function: add vector v2 to vector v1 and store in vector v1
c
c  input arguments:
c     v1       (dp,ar(inca,n),inout) - vector v1
c     inca     (int,sc,in)           - increment (first dimension) of vector v1
c     v2       (dp,ar(incb,n),in)    - vector v2
c     incb     (int,sc,in)           - increment (first dimension) of vector v2
c     n        (int,sc,in)           - length (second dimension) of vectors v1, v2
c
c  output arguments:
c     v1       (dp,ar(inca,n),inout) - vector v1
c
#include "impcom.inc"
c
      integer  inca,incb,n, i,i1,i2
      double precision  v1(*),v2(*)
c
      i1 = 1
      i2 = 1
      do i = 1,n
         v1(i1) = v1(i1) + v2(i2)
         i1 = i1 + inca
         i2 = i2 + incb
      enddo
c
      return
      end
