      function v_sum (v,n)
c --- primary function: sum the components of a vector
c
c  input arguments:
c     v        (dp,ar(n),in)     - vector v
c     n        (int,sc,in)       - length of vector va
c
c  output arguments:
c     v_sum    (dp,sc,out)       - vector sum
c
#include "impcom.inc"
c
      integer  n, i
      double precision   v(n),v_sum, con
 
      con = 0.0d0
      do i = 1,n
         con = con + v(i)
      enddo
      v_sum = con

      return
      end
