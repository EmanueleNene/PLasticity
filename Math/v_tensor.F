      subroutine v_tensor (a,b,C)
c --- primary function: calculate c = a (x) b
c
c  input arguments:
c     a        (dp,ar(3),in)     - first vector to be cross-multiplied
c     b        (dp,ar(3),in)     - second vector to be cross-multiplied
c                                             
c  output arguments:
c     C        (dp,ar(3,3),out)  - resulting tensor
c
#include "impcom.inc"
c
      double precision  a(3),b(3),c(3,3)
      integer i,j
      
      do i = 1,3
          do j = 1,3
             c(j,i) = a(j)*b(i)
          enddo
      enddo

      return
      end
