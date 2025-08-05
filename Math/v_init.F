      subroutine v_init (v,n,const)
c --- primary function: initialize a vector to a constant
c
c  input arguments:
c     n        (int,sc,in)       - length of vector v
c     const    (dp,sc,in)        - constant to set vector v to
c
c  output arguments:
c     v        (dp,ar(n),out)    - vector v
c
#include "impcom.inc"
c
      integer   n,i
      double precision   v(n),const
 
      do i = 1,n
         v(i) = const
      enddo

      return
      end
