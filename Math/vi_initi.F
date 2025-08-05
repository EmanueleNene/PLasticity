      subroutine vi_initi (ivect,inc,n,icon)
c --- primary function: set an integer vector ivect to a constant
c
c  input arguments:
c     inc      (int,sc,in)         - increment (first dimension) of vector ivect
c     n        (int,sc,in)         - number of times to be repeated
c     const    (int,sc,in)         - constant to set components of vector ivect to
c
c  output arguments:
c     ivect    (int,ar(n),inout)  - vector
c
#include "impcom.inc"
c
      integer   n,ivect(*),icon,i,ia,inc
      
      ia = 1
      do i = 1,n
         ivect(ia) = icon
         ia = ia + inc
      enddo

      return
      end
