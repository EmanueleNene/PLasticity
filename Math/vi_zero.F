      subroutine vi_zero (ivect,n)
c --- primary function: set an integer vector ivect to zero
c
c  input arguments:
c     n           (int,sc,in)        - vector lenght
c
c  output arguments:
c     ivect       (int,ar(n),inout)  - vector
c
#include "impcom.inc"
c
      integer  n,ivect(n), i

      do i = 1,n
         ivect(i) = 0
      enddo

      return
      end
