      subroutine v_norm (v,n)
c --- primary function: normalize a vector to unit length
c
c  input arguments:
c     v        (dp,ar(n),inout)  - vector v
c     n        (int,sc,inout)    - dimension length of vector v
c
c  output arguments:
c     v        (dp,ar(n),inout)  - normalized vector v
c     n        (int,sc,inout)    - n = 0 if error in operation
c
#include "impcom.inc"
c
      integer  n, i
      double precision   v(n), con

      con = 0.0d0
      do i = 1,n
         con = con + v(i)*v(i)
      enddo

      con = sqrt(con)

      if (con .lt. TINY) then
          n = 0
      else
          con = 1.0d0/con
          do i = 1,n
             v(i) = v(i)*con
          enddo
      endif

      return
      end
