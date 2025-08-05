      subroutine vi_init (ivect,n,icon)
c --- primary function: set an integer vector ivect to a constant
c
c  input arguments:
c     n           (int,sc,in)        - vector lenght
c     icon        (int,sc,in)        - constant
c
c  output arguments:
c     ivect       (int,ar(n),inout)  - vector
c
#include "impcom.inc"
c
      integer   n,ivect(n),icon, i
 
      do i = 1,n
         ivect(i) = icon
      enddo

      return
      end
