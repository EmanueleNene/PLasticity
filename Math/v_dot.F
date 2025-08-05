      function  v_dot (v1,v2,n)
c --- primary function: compute dot product of vectors v1 and v2
c
c  input arguments:
c     v1       (dp,ar(n),in)     - vector v1
c     v2       (dp,ar(n),in)     - vector v2
c     n        (int,sc,in)       - length of vectors v1 and v2
c
c  output arguments:
c     v_dot    (dp,sc,out)       - dot product of v1 and v2
c
#include "impcom.inc"
c
      integer  n, i
      double precision  v1(n),v2(n),v_dot,con
 
      con = 0.0d0

      do i = 1,n
         con = con + v1(i)*v2(i)
      enddo

      v_dot = con

      return
      end
