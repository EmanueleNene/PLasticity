      subroutine v_zero (v,n)
c --- primary function: initialize a vector to zero
c
c  input arguments:
c     n     (int,sc,in)       - number of words to zero out
c
c  output arguments:
c     v     (dp,ar(n),inout)  - zeroed vector
c
#include "impcom.inc"
c
      integer  n, i
      double precision  v(n)

      do i = 1,n
         v(i) = 0.0d0
      enddo

      return
      end
