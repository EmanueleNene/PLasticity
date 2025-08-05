      subroutine v_cross (a,b,c)
c --- primary function: calculate c = a x b
c
c  input arguments:
c     a        (dp,ar(3),in)     - first vector to be cross-multiplied
c     b        (dp,ar(3),in)     - second vector to be cross-multiplied
c                                             
c  output arguments:
c     c        (dp,ar(3),out)    - resulting vector
c
c
#include "impcom.inc"
c
      double precision  a(3),b(3),c(3)
 
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      return
      end
